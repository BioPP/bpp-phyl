// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMLIKELIHOOD_DF_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMLIKELIHOOD_DF_H

#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/NumTools.h>

#include "../DataFlow/LikelihoodCalculation.h"
#include "../DataFlow/TransitionMatrix.h"
#include "HmmLikelihoodComputation.h"
#include "HmmPhyloEmissionProbabilities.h"

// From the STL:
#include <vector>
#include <memory>

namespace bpp
{
/**
 * @brief A simple implementation of hidden Markov models recursion,
 * in DataFlow implementation.
 *
 * This class builds DF objects linked with HMM computation, and
 * owns a HmmLikelihoodComputation object.
 *
 */

class HmmLikelihood_DF :
  public AlignedLikelihoodCalculation
{
private:
  Context& context_;

protected:
  /**
   * @brief The alphabet describing the hidden states.
   */

  std::shared_ptr<HmmStateAlphabet> hiddenAlphabet_;
  std::shared_ptr<HmmPhyloEmissionProbabilities> emissionProbabilities_;

  /************************
   * DF objects & their targets
   *
   */

  /**
   * DF TransitionMatrix
   *
   */

  std::shared_ptr<ConfiguredTransitionMatrix> matrix_;

  /**
   * DF equilibrium vector from transitionmatrix for computation
   *
   */

  ValueRef<Eigen::VectorXd> hmmEq_;

  /**
   * DF Matrix from transitionmatrix for computation
   *
   */

  ValueRef<Eigen::MatrixXd> hmmTrans_;

  /**
   * DF Matrix from emission likelihoods for computation
   *
   */

  ValueRef<MatrixLik> hmmEmis_;

  /**
   * DF Conditional Likelihoods for sites:
   *
   */

  ValueRef<RowLik> forwardLik_;


  /**
   * DF Backward Likelihoods for sites per state
   *
   * backwardLik_(i,j) corresponds to Pr(x_{j+1}...x_n | y_j=i)/Pr(x_{j+1}|x_1...x_j)
   * where the x are the observed states, and y the hidden states.
   *
   */

  ValueRef<Eigen::MatrixXd> backwardLik_;

  /**
   * Hidden Posterior Probabilities
   *
   * hiddenPostProb_(i,j) corresponds to Pr(y_j=i | x_1...x_n)
   * where the x are the observed states, and y the hidden states.
   *
   */

  ValueRef<Eigen::MatrixXd> hiddenPostProb_;

  Eigen::Index nbStates_, nbSites_;

public:
  /**
   * @brief Build a new HmmLikelihood_DF object.
   *
   * @warning the HmmTransitionMatrix and HmmEmissionProbabilities
   * object passed as argument must be non-null and point toward the
   * same HmmStateAlphabet instance.
   *
   */

  HmmLikelihood_DF(
      Context& context,
      std::shared_ptr<HmmStateAlphabet> hiddenAlphabet,
      std::shared_ptr<HmmTransitionMatrix> transitionMatrix,
      std::shared_ptr<HmmPhyloEmissionProbabilities> emissionProbabilities,
      const std::string& prefix = "");

  HmmLikelihood_DF(const HmmLikelihood_DF& lik) :
    AlignedLikelihoodCalculation(lik),
    context_(lik.context_),
    hiddenAlphabet_(lik.hiddenAlphabet_),
    matrix_(lik.matrix_),
    hmmEq_(lik.hmmEq_),
    hmmTrans_(lik.hmmTrans_),
    hmmEmis_(lik.hmmEmis_),
    forwardLik_(lik.forwardLik_),
    backwardLik_(lik.backwardLik_),
    hiddenPostProb_(lik.hiddenPostProb_),
    nbStates_(lik.nbStates_),
    nbSites_(lik.nbSites_)
  {}

  virtual ~HmmLikelihood_DF() {}

  HmmLikelihood_DF* clone() const { return new HmmLikelihood_DF(*this); }

  void makeLikelihoods() {}

public:
  const HmmStateAlphabet& hmmStateAlphabet() const { return *hiddenAlphabet_; }

  HmmStateAlphabet& hmmStateAlphabet() { return *hiddenAlphabet_; }

  /*
   *@ brief Access to the Transition Matrix
   *
   * !! No check if DF up to date
   *
   */
  const Eigen::MatrixXd& getHmmTransitionMatrix() const { return hmmTrans_->accessValueConst(); }

  const HmmPhyloEmissionProbabilities& getHmmEmissionProbabilities() const { return *emissionProbabilities_; }

  HmmPhyloEmissionProbabilities& getHmmEmissionProbabilities() { return *emissionProbabilities_; }

  void setParameters(const ParameterList& pl)
  {
    setParametersValues(pl);
  }

  void setNamespace(const std::string& nameSpace);

  const Eigen::MatrixXd& getHiddenStatesPosteriorProbabilities() const
  {
    return hiddenPostProb_->targetValue();
  }

  Eigen::VectorXd getHiddenStatesPosteriorProbabilitiesForASite(size_t site) const
  {
    auto& mat = hiddenPostProb_->targetValue();
    return mat.col(Eigen::Index(site));
  }

  DataLik getLikelihoodForASite(size_t site) const
  {
    auto vec = getHiddenStatesPosteriorProbabilitiesForASite(site);

    return hmmEmis_->accessValueConst().col(int(site)).dot(vec);
  }
};
}
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMLIKELIHOOD_DF_H
