//
// File: HmmProcessPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: lundi 23 septembre 2013, à 22h 56
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

#ifndef _HMM_PROCESS_PHYLOLIKELIHOOD_H_
#define _HMM_PROCESS_PHYLOLIKELIHOOD_H_


#include "MultiProcessSequencePhyloLikelihood.h"

#include "../HmmSequenceEvolution.h"

#include "HmmPhyloEmissionProbabilities.h"

// From SeqLib:
#include <Bpp/Seq/Container/AlignedValuesContainer.h>

#include "HmmLikelihood_DF.h"
#include <Bpp/Numeric/Hmm/HmmTransitionMatrix.h>


namespace bpp
{

/**
 * @brief Likelihood framework based on a hmm of simple likelihoods
 *
 * The resulting likelihood is the likelihood of the given Hmm with
 * the site emission probabilities proportional to the computed
 * likelihoods of the process.
 */
  class HmmProcessPhyloLikelihood :
    public MultiProcessSequencePhyloLikelihood
  {
  private:
    std::shared_ptr<HmmPhyloEmissionProbabilities> Hpep_;

    /*
     * @brief LikelihoodCalculation in context of HMM.
     *
     */
  
    mutable std::shared_ptr<HmmLikelihood_DF> hmm_;

  public:
    HmmProcessPhyloLikelihood(
      const AlignedValuesContainer& data,
      HmmSequenceEvolution& processSeqEvol,
      CollectionNodes& collNodes,
      size_t nSeqEvol = 0,
      size_t nData = 0,
      bool verbose = true,
      bool patterns = true);

    HmmProcessPhyloLikelihood(const HmmProcessPhyloLikelihood& mlc) :
      AbstractPhyloLikelihood(mlc),
      AbstractAlignedPhyloLikelihood(mlc),
      MultiProcessSequencePhyloLikelihood(mlc),
      Hpep_(std::shared_ptr<HmmPhyloEmissionProbabilities>(mlc.Hpep_->clone())),
      hmm_(std::shared_ptr<HmmLikelihood_DF>(mlc.hmm_->clone())) {}

    virtual ~HmmProcessPhyloLikelihood() {}

    HmmProcessPhyloLikelihood* clone() const { return new HmmProcessPhyloLikelihood(*this); }

  public:
    void setNamespace(const std::string& nameSpace);

    void fireParameterChanged(const ParameterList& parameters);

    /**
     * @name The likelihood functions.
     *
     * @{
     */

    std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation () const
    {
      return hmm_;
    }

    std::shared_ptr<AlignedLikelihoodCalculation> getAlignedLikelihoodCalculation () const
    {
      return hmm_;
    }

    /*
     *@brief return the posterior probabilities of subprocess on each site.
     *
     *@return MatrixXd sites x states
     */

    VVdouble getPosteriorProbabilitiesPerSitePerProcess() const
    {
      VVdouble pp;
      auto mat = hmm_->getHiddenStatesPosteriorProbabilities().transpose();
      copyEigenToBpp(mat, pp);
      return pp;
    }

    const Eigen::MatrixXd& getHmmTransitionMatrix() const
    {
      return hmm_->getHmmTransitionMatrix();
    }

  };
} // end of namespace bpp.

#endif  // _HMM_PROCESS_PHYLOLIKELIHOOD_H_
