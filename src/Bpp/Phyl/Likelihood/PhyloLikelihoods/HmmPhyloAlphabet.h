// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMPHYLOALPHABET_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMPHYLOALPHABET_H


// From Numeric
#include <Bpp/Numeric/Hmm/HmmStateAlphabet.h>

#include "AlignedPhyloLikelihoodSet.h"
#include "MultiProcessSequencePhyloLikelihood.h"
#include "SingleProcessPhyloLikelihood.h"

namespace bpp
{
/**
 * @brief Hidden states alphabet.
 *
 * Implementation of HmmStateAlphabet where Alphabet States are all
 * the PhyloLikelihoods belonging to a AlignedPhyloLikelihoodSet.
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
  std::vector<std::shared_ptr<AlignedPhyloLikelihoodInterface>> vAP_;

  size_t nbSites_;

public:
  HmmPhyloAlphabet(AlignedPhyloLikelihoodSetInterface& soap) :
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
