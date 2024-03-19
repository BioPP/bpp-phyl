// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SEQUENCEPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SEQUENCEPHYLOLIKELIHOOD_H


// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignmentData.h>

#include "PhyloLikelihood.h"
#include "AbstractPhyloLikelihood.h"
#include "SingleDataPhyloLikelihood.h"
#include "../SequenceEvolution.h"

namespace bpp
{
/**
 * @brief PhyloLikelihoods based on Sequence Evolution class, ie for
 * which there is a (set of) processes able to build a sequence.
 *
 */

/**
 * @brief Interface
 */
class SequencePhyloLikelihoodInterface :
  public virtual SingleDataPhyloLikelihoodInterface
{
protected:
  virtual ~SequencePhyloLikelihoodInterface() {}

  SequencePhyloLikelihoodInterface* clone() const override = 0;

public:
  /**
   * @brief The Sequence Evolution
   *
   */
  virtual const SequenceEvolution& sequenceEvolution() const = 0;

  virtual std::shared_ptr<const SequenceEvolution> getSequenceEvolution() const = 0;

  virtual const size_t getSequenceEvolutionNumber() const = 0;
};

class AbstractSequencePhyloLikelihood :
  public virtual SequencePhyloLikelihoodInterface,
  public virtual AbstractSingleDataPhyloLikelihood
{
protected:
  /**
   * @brief the Sequence Evolution
   */
  std::shared_ptr<SequenceEvolution> seqEvol_;

  /**
   * @brief the Sequence Evolution number
   */
  size_t nSeqEvol_;

public:
  AbstractSequencePhyloLikelihood(
      Context& context,
      std::shared_ptr<SequenceEvolution> se,
      size_t nSE = 0,
      size_t nData = 0) :
    AbstractPhyloLikelihood(context),
    AbstractSingleDataPhyloLikelihood(
      context, 0,
      (se->getSubstitutionProcessNumbers().size() != 0) ? se->substitutionProcess(se->getSubstitutionProcessNumbers()[0]).getNumberOfStates() : 0, nData),
    seqEvol_(se),
    nSeqEvol_(nSE)
  {}

  AbstractSequencePhyloLikelihood(const AbstractSequencePhyloLikelihood& aspl) :
    AbstractSingleDataPhyloLikelihood(aspl),
    seqEvol_(aspl.seqEvol_),
    nSeqEvol_(aspl.nSeqEvol_)
  {}

  AbstractSequencePhyloLikelihood& operator=(const AbstractSequencePhyloLikelihood& aspl)
  {
    AbstractSingleDataPhyloLikelihood::operator=(aspl);
    seqEvol_ = aspl.seqEvol_;
    nSeqEvol_ = aspl.nSeqEvol_;
    return *this;
  }

public:
  /**
   * @brief The Sequence Evolution
   *
   */
  const SequenceEvolution& sequenceEvolution() const override
  {
    return *seqEvol_;
  }

  std::shared_ptr<const SequenceEvolution> getSequenceEvolution() const override
  {
    return seqEvol_;
  }

  const size_t getSequenceEvolutionNumber() const override
  {
    return nSeqEvol_;
  }

  /**
   * @name the Likelihood interface
   *
   */
  std::shared_ptr<const Alphabet> getAlphabet() const override
  {
    if (seqEvol_->getSubstitutionProcessNumbers().size() == 0)
      return NULL;
    else
      return seqEvol_->substitutionProcess(seqEvol_->getSubstitutionProcessNumbers()[0]).getModel(0, 0)->getAlphabet();
  }
};


class AbstractParametrizableSequencePhyloLikelihood :
  public virtual AbstractSequencePhyloLikelihood,
  public virtual AbstractParametrizable
{
public:
  AbstractParametrizableSequencePhyloLikelihood(
      Context& context,
      std::shared_ptr<SequenceEvolution> se,
      size_t nSE = 0,
      size_t nData = 0) :
    AbstractPhyloLikelihood(context),
    AbstractSequencePhyloLikelihood(context, se, nSE, nData),
    AbstractParametrizable("")
  {
    // initialize INDEPENDENT parameters:
    addParameters_(se->getIndependentParameters());
  }

  AbstractParametrizableSequencePhyloLikelihood(const AbstractParametrizableSequencePhyloLikelihood& apspl) :
    AbstractSequencePhyloLikelihood(apspl),
    AbstractParametrizable(apspl)
  {}

  AbstractParametrizableSequencePhyloLikelihood& operator=(const AbstractParametrizableSequencePhyloLikelihood& apspl)
  {
    AbstractSequencePhyloLikelihood::operator=(apspl);
    AbstractParametrizable::operator=(apspl);
    return *this;
  }

protected:
  virtual void fireParameterChanged(const ParameterList& parameters) override
  {
    seqEvol_->matchParametersValues(parameters);
  }

public:
  void setData(std::shared_ptr<const AlignmentDataInterface> sites, size_t nData = 0) override
  {
    AbstractSingleDataPhyloLikelihood::setData(sites, nData);
  }

  void setNamespace(const std::string& nameSpace) override
  {
    seqEvol_->setNamespace(nameSpace);
  }


  ParameterList getNonDerivableParameters() const override
  {
    return seqEvol_->getNonDerivableParameters();
  }

  ParameterList getDerivableParameters() const override
  {
    return seqEvol_->getBranchLengthParameters(true);
  }

  ParameterList getBranchLengthParameters() const override
  {
    return seqEvol_->getBranchLengthParameters(true);
  }

  ParameterList getRootFrequenciesParameters() const override
  {
    return seqEvol_->getRootFrequenciesParameters(true);
  }

  ParameterList getRateDistributionParameters() const override
  {
    return seqEvol_->getRateDistributionParameters(true);
  }

  ParameterList getSubstitutionModelParameters() const override
  {
    return seqEvol_->getSubstitutionModelParameters(true);
  }

  virtual ParameterList getSubstitutionProcessParameters() const
  {
    return seqEvol_->getSubstitutionProcessParameters(true);
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SEQUENCEPHYLOLIKELIHOOD_H
