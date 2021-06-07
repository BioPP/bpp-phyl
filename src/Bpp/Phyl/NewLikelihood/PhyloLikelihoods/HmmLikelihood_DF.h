//
// File: HmmLikelihood_DF.h
// Created by: Laurent Guéguen
// Created on: jeudi 13 août 2020, à 17h 46
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use, 
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info". 

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability. 

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or 
  data to be ensured and,  more generally, to use and operate it in the 
  same conditions as regards security. 

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _HMMLIKELIHOOD_DF_H_
#define _HMMLIKELIHOOD_DF_H_

#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/NumTools.h>
#include <Bpp/Numeric/Matrix/Matrix.h>

#include "../DataFlow/TransitionMatrix.h"
#include "../DataFlow/LikelihoodCalculation.h"

#include "HmmPhyloEmissionProbabilities.h"

#include "HmmLikelihoodComputation.h"


//From the STL:
#include <vector>
#include <memory>

namespace bpp {

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
    
    ValueRef<Eigen::MatrixXd> hmmEmis_;

    /**
     * DF Conditional Likelihoods for sites:
     *
     */
    
    ValueRef<Eigen::RowVectorXd> forwardLik_;


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

    HmmLikelihood_DF(const HmmLikelihood_DF& lik):
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
    {
    }

    virtual ~HmmLikelihood_DF() {}

    HmmLikelihood_DF* clone() const { return new HmmLikelihood_DF(*this); }

    void makeLikelihoods() {};

  public:
    const HmmStateAlphabet& getHmmStateAlphabet() const { return *hiddenAlphabet_; }

    HmmStateAlphabet& getHmmStateAlphabet() { return *hiddenAlphabet_; }

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
      return hiddenPostProb_->getTargetValue();
    }

    Eigen::Ref<const Eigen::VectorXd> getHiddenStatesPosteriorProbabilitiesForASite(size_t site) const
    {
      auto& mat = hiddenPostProb_->getTargetValue();
      return mat.col(Eigen::Index(site));
    }

    double getLikelihoodForASite(size_t site) const 
    {
      auto vec = getHiddenStatesPosteriorProbabilitiesForASite(site);

      return vec.dot(hmmEmis_->accessValueConst().col(Eigen::Index(site)));
    }

  };

}

#endif //_HMMLIKELIHOOD_DF_H_

