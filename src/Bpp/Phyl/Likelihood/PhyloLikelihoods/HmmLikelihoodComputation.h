// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMLIKELIHOODCOMPUTATION_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMLIKELIHOODCOMPUTATION_H

#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/NumTools.h>

#include "../DataFlow/TransitionMatrix.h"
#include "HmmPhyloEmissionProbabilities.h"

// From the STL:
#include <vector>
#include <memory>

namespace bpp
{
class CondLikelihood : public Value<Eigen::MatrixXd>
{
private:
  /**
   * @brief Dimension of the data : states X sites
   *
   */

  Dimension<Eigen::MatrixXd> targetDimension_;

public:
  static ValueRef<Eigen::MatrixXd> create (Context& c, NodeRefVec&& deps, const Dimension<Eigen::MatrixXd>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (CondLikelihood), deps);
    checkDependencyVectorSize (typeid (CondLikelihood), deps, 2);

    // dependency on the name, to make objects different
    checkNthDependencyIsValue<std::string>(typeid(CondLikelihood), deps, 1);

    return cachedAs<Value<Eigen::MatrixXd>>(c, std::make_shared<CondLikelihood>(std::move (deps), dim));
  }

  CondLikelihood (NodeRefVec&& deps, const Dimension<Eigen::MatrixXd>& dim)
    : Value<Eigen::MatrixXd>(std::move (deps)), targetDimension_ (dim)
  {
    this->accessValueMutable().resize(dim.rows, dim.cols);
  }

  std::string debugInfo () const override
  {
    using namespace numeric;
    const auto name = accessValueConstCast<std::string>(*this->dependency (1));
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_) + ":name= " + name;
  }

  // ForwardHmmLikelihood_DF additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const CondLikelihood*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    // NodeRef derivNode = this->dependency (3);
    // const auto nDeriv = accessValueConstCast<size_t> (*derivNode);

    throw Exception("CondLikelihood::derive is done in dependency class.");
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return CondLikelihood::create (c, std::move (deps), targetDimension_);
  }

  Eigen::MatrixXd& getCondLikelihood()
  {
    return this->accessValueMutable();
  }

  const Eigen::MatrixXd& getCondLikelihood() const
  {
    return this->accessValueConst();
  }

private:
  // Nothing happens here, computation is done in ForwardHmmLikelihood_DF class
  void compute() override {}

  friend class ForwardHmmLikelihood_DF;
  friend class ForwardHmmDLikelihood_DF;
  friend class ForwardHmmD2Likelihood_DF;
};

/*
 * Computation of Forward Likelihood Arrays
 *
 *
 * Dependencies are:
 *  Value<VectorXd> : Starting vector of states probabililies
 *  Value<MatrixXd> : TransitionMatrix
 *  Value<MatrixLik> : Matrix of Emission likelihoods states X sites
 *
 * After computation, its value stores the conditional forward
 * likelihoods of the sites, P(x_j|x_1,...,x_{j-1}), where the x are
 * the observed states.
 *
 * The conditional matrix of the likelihoods per hidden state
 * Pr(x_1...x_j, y_j=i)/Pr(x_1...x_j) (with y the hidden states), is
 * stored, available through getForwardCondLikelihood.
 *
 */

class ForwardHmmLikelihood_DF : public Value<RowLik>
{
private:
  /*
   * @brief conditional forward likelihoods : Will be used by
   * backward likelihoods computation.
   *
   * condLik_(i,j) corresponds to @f$Pr(x_1...x_j, y_j=i)/Pr(x_1...x_j)@f$,
   * where the x are the observed states, and y the hidden states.
   *
   * @f$ \sum_i \text{condLik\_(i,j)} = 1 @f$
   */

  ValueRef<Eigen::MatrixXd> condLik_;

  /*
   * @brief Conditional partial likelihood, used for computation.
   *
   * parCondLik_(i,j) corresponds to Pr(x_1...x_j, y_{j+1}=i)/Pr(x_1...x_j),
   * where the x are the observed states, and y the hidden states.
   *
   * @f$ \sum_i \text{parCondLik\_(i,j)} = 1 @f$
   */

  std::vector<Eigen::VectorXd> parCondLik_;

  /*
   * @brief Dimension of the data : states X sites
   *
   */

  Dimension<Eigen::MatrixXd> targetDimension_;

public:
  using Self = ForwardHmmLikelihood_DF;

  static ValueRef<RowLik> create (Context& c, NodeRefVec&& deps, const Dimension<Eigen::MatrixXd>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 3);

    checkNthDependencyIsValue<Eigen::VectorXd>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<Eigen::MatrixXd>(typeid (Self), deps, 1);
    checkNthDependencyIsValue<MatrixLik>(typeid (Self), deps, 2);

    auto sself = std::make_shared<Self>(std::move (deps), dim);
    sself->build(c);

    return cachedAs<Value<RowLik>>(c, sself);
  }

  ForwardHmmLikelihood_DF (NodeRefVec&& deps, const Dimension<Eigen::MatrixXd>& dim)
    : Value<RowLik>(std::move (deps)), condLik_(), parCondLik_((size_t)dim.cols), targetDimension_ (dim)
  {
    for (auto& v:parCondLik_)
    {
      v.resize(dim.rows);
    }

    this->accessValueMutable().resize(targetDimension_.cols);
  }

  void build(Context& c)
  {
    auto fname = NumericConstant<std::string>::create(c, "forwardCondLik");

    condLik_ = CondLikelihood::create(c, {this->shared_from_this(), fname}, targetDimension_);

    const auto& hmmEq = dynamic_pointer_cast<Value<Eigen::VectorXd>>(this->dependency(0))->targetValue();

    if (hmmEq.rows() != targetDimension_.rows)
      throw BadSizeException("ForwardHmmLikelihood_DF: bad dimension for starting vector", size_t(hmmEq.rows()), size_t(targetDimension_.rows));


    const auto& hmmTrans = dynamic_pointer_cast<Value<Eigen::MatrixXd>>(this->dependency(1))->targetValue();

    if (hmmTrans.cols() != hmmTrans.rows())
      throw BadSizeException("ForwardHmmLikelihood_DF: Transition matrix should be square", size_t(hmmTrans.cols()), size_t(hmmTrans.rows()));
    if (hmmTrans.rows() != targetDimension_.rows)
      throw BadSizeException("ForwardHmmLikelihood_DF: bad number of rows for transition matrix", size_t(hmmTrans.rows()), size_t(targetDimension_.rows));

    const auto& hmmEmis = dynamic_pointer_cast<Value<MatrixLik>>(this->dependency(2))->targetValue();

    if (hmmEmis.rows() != targetDimension_.rows)
      throw BadSizeException("ForwardHmmLikelihood_DF: bad number of states for emission matrix", size_t(hmmEmis.rows()), size_t(targetDimension_.rows));
    if (hmmEmis.cols() != targetDimension_.cols)
      throw BadSizeException("ForwardHmmLikelihood_DF: bad number of sites for emission matrix", size_t(hmmEmis.cols()), size_t(targetDimension_.cols));
  }

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // ForwardHmmLikelihood_DF additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final;

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

  ValueRef<Eigen::MatrixXd> getForwardCondLikelihood() const
  {
    return condLik_;
  }

  const std::vector<Eigen::VectorXd>& getParCondLik() const
  {
    return parCondLik_;
  }

private:
  void compute() override;
};


/*
 * Computation of 1st order Derived Forward Likelihood Arrays
 *
 * Dependencies are:
 *  Value<VectorXd> : Starting vector of states probabililies
 *  Value<MatrixXd> : TransitionMatrix
 *  Value<MatrixLik> : Matrix of Emission likelihoods states X sites
 *
 *  ForwardHmmLikelihood_DF : Forward Computations
 *
 *  Value<VectorXd> : Derivatives of starting vector of states probabililies
 *  Value<MatrixXd> : Derivatives of TransitionMatrix
 *  Value<MatrixLik> : Derivatives Matrix of Emission likelihoods states X sites
 *
 * After computation, its value stores the derivates of the
 * conditional forward likelihoods of the sites,
 * dP(x_j|x_1,...,x_{j-1}), where the x are the observed states.
 *
 * The derivates of the conditional matrix of the likelihoods per hidden state
 * d(Pr(x_1...x_j, y_j=i)/Pr(x_1...x_j)) (with y the hidden states), is
 * stored, available through getForwardDCondLikelihood.
 *
 */


class ForwardHmmDLikelihood_DF : public Value<RowLik>
{
private:
  /**
   * @brief derivatives of the conditional forward likelihoods :
   * Will be used by likelihoods computation of 2nd order
   *
   * dcondLik_(i,j) corresponds to d(Pr(x_1...x_j, y_j=i)/Pr(x_1...x_j)),
   * where the x are the observed states, and y the hidden states.
   *
   */

  ValueRef<Eigen::MatrixXd> dCondLik_;

  /*
   * @brief Conditional partial likelihood derivatives, used for
   * computation.
   *
   */

  std::vector<Eigen::VectorXd> dParCondLik_;

  /**
   * @brief Dimension of the data : states X sites
   *
   */

  Dimension<Eigen::MatrixXd> targetDimension_;

public:
  using Self = ForwardHmmDLikelihood_DF;

  static ValueRef<RowLik> create (Context& c, NodeRefVec&& deps, const Dimension<Eigen::MatrixXd>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 7);

    checkNthDependencyIsValue<Eigen::VectorXd>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<Eigen::MatrixXd>(typeid (Self), deps, 1);
    checkNthDependencyIsValue<MatrixLik>(typeid (Self), deps, 2);

    checkNthDependencyIs<ForwardHmmLikelihood_DF>(typeid (Self), deps, 3);

    checkNthDependencyIsValue<Eigen::VectorXd>(typeid (Self), deps, 4);
    checkNthDependencyIsValue<Eigen::MatrixXd>(typeid (Self), deps, 5);
    checkNthDependencyIsValue<MatrixLik>(typeid (Self), deps, 6);

    auto sself = std::make_shared<Self>(std::move (deps), dim);
    sself->build(c);

    return cachedAs<Value<RowLik>>(c, sself);
  }

  ForwardHmmDLikelihood_DF (NodeRefVec&& deps, const Dimension<Eigen::MatrixXd>& dim)
    : Value<RowLik>(std::move (deps)), dCondLik_(), dParCondLik_((size_t)dim.cols), targetDimension_ (dim)
  {
    for (auto& v:dParCondLik_)
    {
      v.resize(dim.rows);
    }
    this->accessValueMutable().resize(targetDimension_.cols);
  }

  void build(Context& c)
  {
    auto fname = NumericConstant<std::string>::create(c, "forwardDcondLik");

    dCondLik_ = CondLikelihood::create(c, {this->shared_from_this(), fname}, targetDimension_);
  }

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // ForwardHmmDLikelihood_DF additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final;

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

  ValueRef<Eigen::MatrixXd> getForwardDCondLikelihood() const
  {
    return dCondLik_;
  }

  const std::vector<Eigen::VectorXd>& getParDCondLik() const
  {
    return dParCondLik_;
  }

private:
  void compute() override;
};

/*
 * Computation of 2nd order Derived Forward Likelihood Arrays
 *
 * Dependencies are:
 *  Value<VectorXd> : Starting vector of states probabililies
 *  Value<MatrixXd> : TransitionMatrix
 *  Value<MatrixLik> : Matrix of Emission likelihoods states X sites
 *
 *  ForwardHmmLikelihood_DF : Forward Computations
 *
 *  Value<VectorXd> : 1st Derivatives of starting vector of states probabililies
 *  Value<MatrixXd> : 1st Derivatives of TransitionMatrix
 *  Value<MatrixLik> : 1st Derivatives Matrix of Emission likelihoods states X sites
 *
 *  ForwardHmmDLikelihood_DF : 1st order derivatives Forward Computations
 *
 *  Value<VectorXd> : 2nd Derivatives of starting vector of states probabililies
 *  Value<MatrixXd> : 2nd Derivatives of TransitionMatrix
 *  Value<MatrixLik> : 2nd Derivatives Matrix of Emission likelihoods states X sites
 *
 * After computation, its value stores the 2nd derivates of the
 * conditional forward likelihoods of the sites,
 * d2P(x_j|x_1,...,x_{j-1}), where the x are the observed states.
 *
 * The derivates of the conditional matrix of the likelihoods per hidden state
 * d(Pr(x_1...x_j, y_j=i)/Pr(x_1...x_j)) (with y the hidden states), is
 * stored, available through getForwardCondLikelihood.
 *
 */


class ForwardHmmD2Likelihood_DF : public Value<RowLik>
{
private:
  /**
   * @brief derivatives of the conditional forward likelihoods :
   * Will be used by backward likelihoods computation.
   *
   * d2condLik_(i,j) corresponds to d2(Pr(x_1...x_j,
   * y_j=i)/Pr(x_1...x_j)), where the x are the observed states, and
   * y the hidden states.
   */

  ValueRef<Eigen::MatrixXd> d2CondLik_;

  /**
   * @brief Dimension of the data : states X sites
   *
   */

  Dimension<Eigen::MatrixXd> targetDimension_;

public:
  using Self = ForwardHmmD2Likelihood_DF;

  static ValueRef<RowLik> create (Context& c, NodeRefVec&& deps, const Dimension<Eigen::MatrixXd>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 11);

    checkNthDependencyIsValue<Eigen::VectorXd>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<Eigen::MatrixXd>(typeid (Self), deps, 1);
    checkNthDependencyIsValue<MatrixLik>(typeid (Self), deps, 2);

    checkNthDependencyIs<ForwardHmmLikelihood_DF>(typeid (Self), deps, 3);

    checkNthDependencyIsValue<Eigen::VectorXd>(typeid (Self), deps, 4);
    checkNthDependencyIsValue<Eigen::MatrixXd>(typeid (Self), deps, 5);
    checkNthDependencyIsValue<MatrixLik>(typeid (Self), deps, 6);

    checkNthDependencyIs<ForwardHmmDLikelihood_DF>(typeid (Self), deps, 7);

    checkNthDependencyIsValue<Eigen::VectorXd>(typeid (Self), deps, 8);
    checkNthDependencyIsValue<Eigen::MatrixXd>(typeid (Self), deps, 9);
    checkNthDependencyIsValue<MatrixLik>(typeid (Self), deps, 10);

    auto sself = std::make_shared<Self>(std::move (deps), dim);
    sself->build(c);

    return cachedAs<Value<RowLik>>(c, sself);
  }

  ForwardHmmD2Likelihood_DF (NodeRefVec&& deps, const Dimension<Eigen::MatrixXd>& dim)
    : Value<RowLik>(std::move (deps)), d2CondLik_(), targetDimension_ (dim)
  {
    this->accessValueMutable().resize(targetDimension_.cols);
  }

  void build(Context& c)
  {
    auto fname = NumericConstant<std::string>::create(c, "forwardD2condLik");

    d2CondLik_ = CondLikelihood::create(c, {this->shared_from_this(), fname}, targetDimension_);
  }

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // ForwardHmmD2Likelihood_DF additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    throw Exception("ForwardHmmD2Likelihood_DF::derive not implemented.");
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute() override;
};


/////////////////////////////////////////////////////////////////////////

/*
 * Computation of Backward Likelihood Arrays
 *
 *
 * Dependencies are:
 *  Value<RowLik> : Vector of conditional Forward Likelihoods
 *  Value<MatrixXd> : TransitionMatrix
 *  Value<MatrixLik> : Matrix of Emission likelihoods states X sites
 *  Value<size_t> : level of derivation
 *
 * After computation, stores the conditional likelihoods of the
 * sites for all states.
 */

class BackwardHmmLikelihood_DF : public Value<Eigen::MatrixXd>
{
private:
  /**
   * @brief backward likelihood
   *
   * Its value stores the ratio @f$ \frac{Pr(x_{j+1}...x_n |
   * y_j=i)}{P(x_{j+1}...x_n |x_1,...,x_j)}@f$ of the backward
   * conditional likelihoods per hidden state divided per conditional
   * state likelihood, where the x are the observed states, and y the
   * hidden states,
   *
   */

  /**
   * @brief Dimension of the data : states X sites
   *
   */

  Dimension<Eigen::MatrixXd> targetDimension_;

public:
  using Self = BackwardHmmLikelihood_DF;

  static ValueRef<Eigen::MatrixXd> create (Context& c, NodeRefVec&& deps, const Dimension<Eigen::MatrixXd>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 3);

    checkNthDependencyIsValue<RowLik>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<Eigen::MatrixXd>(typeid (Self), deps, 1);
    checkNthDependencyIsValue<MatrixLik>(typeid (Self), deps, 2);

    return cachedAs<Value<Eigen::MatrixXd>>(c, std::make_shared<Self>(std::move (deps), dim));
  }

  BackwardHmmLikelihood_DF (NodeRefVec&& deps, const Dimension<Eigen::MatrixXd>& dim)
    : Value<Eigen::MatrixXd>(std::move (deps)), targetDimension_ (dim)
  {
    this->accessValueMutable().resize(dim.rows, dim.cols);

    const auto& hmmScale = accessValueConstCast<RowLik>(*this->dependency(0));
    if (hmmScale.cols() != dim.cols)
      throw BadSizeException("BackwardHmmLikelihood_DF: bad dimension for forward likelihoods vector", size_t(hmmScale.cols()), size_t(dim.cols));


    const auto& hmmTrans = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(1));
    if (hmmTrans.cols() != dim.rows)
      throw BadSizeException("BackwardHmmLikelihood_DF: bad size for transition matrix", size_t(hmmTrans.cols()), size_t(dim.rows));

    const auto& hmmEmis = accessValueConstCast<MatrixLik>(*this->dependency(2));
    if (hmmEmis.rows() != dim.rows)
      throw BadSizeException("BackwardHmmLikelihood_DF: bad number of states for emission matrix", size_t(hmmEmis.rows()), size_t(dim.rows));
    if (hmmEmis.cols() != dim.cols)
      throw BadSizeException("BackwardHmmLikelihood_DF: bad number of sites for emission matrix", size_t(hmmEmis.cols()), size_t(dim.cols));
  }

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // BackwardHmmLikelihood_DF additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    // NodeRef derivNode = this->dependency (3);
    // const auto nDeriv = accessValueConstCast<size_t> (*derivNode);

    throw Exception("BackwardHmmLikelihood_DF::derive To be finished.");

    return Self::create (c, {this->dependency(0)->derive (c, node)}, targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute() override;
};
}
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMLIKELIHOODCOMPUTATION_H
