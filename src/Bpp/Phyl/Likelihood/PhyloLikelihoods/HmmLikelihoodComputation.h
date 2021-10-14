//
// File: HmmLikelihoodComputation.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: jeudi 13 aoÃÂ»t 2020, ÃÂ  17h 46
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

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
class CondLikelihood : public Value<MatrixLik>
{
private:
  /**
   * @brief Dimension of the data : states X sites
   *
   */

  Dimension<MatrixLik> targetDimension_;

public:
  static ValueRef<MatrixLik> create (Context& c, NodeRefVec&& deps, const Dimension<MatrixLik>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (CondLikelihood), deps);
    checkDependencyVectorSize (typeid (CondLikelihood), deps, 2);

    // dependency on the name, to make objects different
    checkNthDependencyIsValue<std::string>(typeid(CondLikelihood), deps, 1);

    return cachedAs<Value<MatrixLik> >(c, std::make_shared<CondLikelihood>(std::move (deps), dim));
  }

  CondLikelihood (NodeRefVec&& deps, const Dimension<MatrixLik>& dim)
    : Value<MatrixLik>(std::move (deps)), targetDimension_ (dim)
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

  MatrixLik& getCondLikelihood()
  {
    return this->accessValueMutable();
  }

  const MatrixLik& getCondLikelihood() const
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
  /**
   * @brief conditional forward likelihoods : Will be used by
   * backward likelihoods computation.
   *
   * condLik_(i,j) corresponds to Pr(x_1...x_j, y_j=i)/Pr(x_1...x_j),
   * where the x are the observed states, and y the hidden states.
   */

  typedef std::vector<VectorLik> vVectorLik;

  ValueRef<MatrixLik> condLik_;

  /*
   * @brief Conditional partial likelihood, used for computation.
   *
   * parCondLik_(i,j) corresponds to Pr(x_1...x_j, y_{j+1}=i)/Pr(x_1...x_j),
   * where the x are the observed states, and y the hidden states.
   */

  std::vector<VectorLik> parCondLik_;

  /*
   * @brief Dimension of the data : states X sites
   *
   */

  Dimension<MatrixLik> targetDimension_;

public:
  using Self = ForwardHmmLikelihood_DF;

  static ValueRef<RowLik> create (Context& c, NodeRefVec&& deps, const Dimension<MatrixLik>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 3);

    checkNthDependencyIsValue<Eigen::VectorXd>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<Eigen::MatrixXd>(typeid (Self), deps, 1);
    checkNthDependencyIsValue<MatrixLik>(typeid (Self), deps, 2);

    auto sself = std::make_shared<Self>(std::move (deps), dim);
    sself->build(c);

    return cachedAs<Value<RowLik> >(c, sself);
  }

  ForwardHmmLikelihood_DF (NodeRefVec&& deps, const Dimension<MatrixLik>& dim)
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

    const auto& hmmEq = dynamic_pointer_cast<Value<Eigen::VectorXd> >(this->dependency(0))->getTargetValue();

    if (hmmEq.rows() != targetDimension_.rows)
      throw BadSizeException("ForwardHmmLikelihood_DF: bad dimension for starting vector", size_t(hmmEq.rows()), size_t(targetDimension_.rows));


    const auto& hmmTrans = dynamic_pointer_cast<Value<Eigen::MatrixXd> >(this->dependency(1))->getTargetValue();

    if (hmmTrans.cols() != hmmTrans.rows())
      throw BadSizeException("ForwardHmmLikelihood_DF: Transition matrix should be square", size_t(hmmTrans.cols()), size_t(hmmTrans.rows()));
    if (hmmTrans.rows() != targetDimension_.rows)
      throw BadSizeException("ForwardHmmLikelihood_DF: bad number of rows for transition matrix", size_t(hmmTrans.rows()), size_t(targetDimension_.rows));

    const auto& hmmEmis = dynamic_pointer_cast<Value<MatrixLik> >(this->dependency(2))->getTargetValue();

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

  ValueRef<MatrixLik> getForwardCondLikelihood() const
  {
    return condLik_;
  }

  const std::vector<VectorLik>& getParCondLik() const
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

  ValueRef<MatrixLik> dCondLik_;

  /*
   * @brief Conditional partial likelihood derivatives, used for
   * computation.
   *
   */

  std::vector<VectorLik> dParCondLik_;

  /**
   * @brief Dimension of the data : states X sites
   *
   */

  Dimension<MatrixLik> targetDimension_;

public:
  using Self = ForwardHmmDLikelihood_DF;

  static ValueRef<RowLik> create (Context& c, NodeRefVec&& deps, const Dimension<MatrixLik>& dim)
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

    return cachedAs<Value<RowLik> >(c, sself);
  }

  ForwardHmmDLikelihood_DF (NodeRefVec&& deps, const Dimension<MatrixLik>& dim)
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

  ValueRef<MatrixLik> getForwardDCondLikelihood() const
  {
    return dCondLik_;
  }

  const std::vector<VectorLik>& getParDCondLik() const
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

  ValueRef<MatrixLik> d2CondLik_;

  /**
   * @brief Dimension of the data : states X sites
   *
   */

  Dimension<MatrixLik> targetDimension_;

public:
  using Self = ForwardHmmD2Likelihood_DF;

  static ValueRef<RowLik> create (Context& c, NodeRefVec&& deps, const Dimension<MatrixLik>& dim)
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

    return cachedAs<Value<RowLik> >(c, sself);
  }

  ForwardHmmD2Likelihood_DF (NodeRefVec&& deps, const Dimension<MatrixLik>& dim)
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

class BackwardHmmLikelihood_DF : public Value<MatrixLik>
{
private:
  /**
   * @brief backward likelihood
   *
   * Its value stores the backward conditional likelihoods per
   * hidden state Pr(x_{j+1}...x_n | y_j=i), where the x are the
   * observed states, and y the hidden states, divided per
   * conditional state likelihood: P(x_{j+1}...x_n |x_1,...,x_j).
   *
   */

  /**
   * @brief Dimension of the data : states X sites
   *
   */

  Dimension<MatrixLik> targetDimension_;

public:
  using Self = BackwardHmmLikelihood_DF;

  static ValueRef<MatrixLik> create (Context& c, NodeRefVec&& deps, const Dimension<MatrixLik>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 3);

    checkNthDependencyIsValue<RowLik>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<Eigen::MatrixXd>(typeid (Self), deps, 1);
    checkNthDependencyIsValue<MatrixLik>(typeid (Self), deps, 2);

    return cachedAs<Value<MatrixLik> >(c, std::make_shared<Self>(std::move (deps), dim));
  }

  BackwardHmmLikelihood_DF (NodeRefVec&& deps, const Dimension<MatrixLik>& dim)
    : Value<MatrixLik>(std::move (deps)), targetDimension_ (dim)
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
