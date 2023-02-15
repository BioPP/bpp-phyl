//
// File: HmmPhyloAlphabet.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: mardi 27 octobre 2015, ÃÂ  19h 09
//

/*
  Copyright or ÃÂ© or Copr. CNRS, (November 16, 2004)
  
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

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMPHYLOALPHABET_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMPHYLOALPHABET_H



// From Numeric
#include <Bpp/Numeric/Hmm/HmmStateAlphabet.h>

#include "SetOfAlignedPhyloLikelihood.h"
#include "MultiProcessSequencePhyloLikelihood.h"
#include "SingleProcessPhyloLikelihood.h"

namespace bpp
{
/**
 * @brief Hidden states alphabet.
 *
 * Implementation of HmmStateAlphabet where Alphabet States are all
 * the PhyloLikelihoods belonging to a SetOfAlignedPhyloLikelihood.
 */
class HmmPhyloAlphabet :
  public virtual HmmStateAlphabet,
  public AbstractParametrizable
{
private:
  Context& context_;

  /**
   * @brief vector of aligned phylolikelihoods. These
   * phylolikelihoods are not owned by the alphabet.
   */
  std::vector<std::shared_ptr<AlignedPhyloLikelihoodInterface> > vAP_;

  size_t nbSites_;

public:
  HmmPhyloAlphabet(SetOfAlignedPhyloLikelihoodInterface& soap) :
    AbstractParametrizable(""),
    context_(soap.context()),
    vAP_(),
    nbSites_(0)
  {
    const std::vector<size_t>& nphyl = soap.getNumbersOfPhyloLikelihoods();

    for (size_t i = 0; i < nphyl.size(); ++i)
    {
      auto ap = soap.getAlignedPhyloLikelihood(nphyl[i]);
      vAP_.push_back(ap);
      includeParameters_(ap->getParameters());
    }

    nbSites_ = soap.getNumberOfSites();
  }

  HmmPhyloAlphabet(MultiProcessSequencePhyloLikelihood& mpsp) :
    AbstractParametrizable(""),
    context_(mpsp.context()),
    vAP_(),
    nbSites_(0)
  {
    auto nb = mpsp.getNumberOfSubstitutionProcess();
    ParameterList pl;
    for (size_t i = 0; i < nb; i++)
    {
      // no parameters for each SingleProcessPhyloLikelihood
      auto spl = make_shared<SingleProcessPhyloLikelihood>(getContext(), mpsp.getLikelihoodCalculationForAProcess(i), pl);

      if (!nbSites_)
        nbSites_ = spl->getNumberOfSites();
      else if (nbSites_ != spl->getNumberOfSites())
        throw BadSizeException("HmmPhyloAlphabet::HmmPhyloAlphabet : different numbers of sites", spl->getNumberOfSites(), nbSites_);

      vAP_.push_back(spl);
    }

    includeParameters_(mpsp.getParameters());
  }


  HmmPhyloAlphabet(const HmmPhyloAlphabet& hpa) :
    AbstractParametrizable(hpa),
    context_(hpa.context_),
    vAP_(hpa.vAP_),
    nbSites_(hpa.nbSites_)
  {}

  HmmPhyloAlphabet* clone() const override { return new HmmPhyloAlphabet(*this);}

  virtual ~HmmPhyloAlphabet() {}

  Context& getContext() { return context_;}

  /**
   * @param stateIndex The index of a hidden state.
   * @return The corresponding hidden state.
   * @see getNumberOfStates
   */
  const Clonable& getState(size_t stateIndex) const override
  {
    return *vAP_[stateIndex];
  }

  const AlignedPhyloLikelihoodInterface& alignedPhyloLikelihood(size_t stateIndex) const
  {
    return *vAP_[stateIndex];
  }

  std::shared_ptr<const AlignedPhyloLikelihoodInterface> getAlignedPhyloLikelihood(size_t stateIndex) const
  {
    return vAP_[stateIndex];
  }

  size_t getNumberOfStates() const override
  {
    return vAP_.size();
  }

  size_t getNumberOfSites() const
  {
    return nbSites_;
  }

  /**
   * @brief Tell if this instance can work with the instance of alphabet given as input.
   *
   * In many case, this will return true if the pointer provided as argument refers to this object.
   *
   * @param stateAlphabet The alphabet to check.
   * @return true if the matrix is compatible with the given alphabet.
   */
  bool worksWith(const HmmStateAlphabet& stateAlphabet) const override
  {
    return &stateAlphabet == this;
  }
};
}
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMPHYLOALPHABET_H
