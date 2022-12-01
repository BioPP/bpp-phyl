//
// File: SequencePhyloLikelihood.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: mardi 28 avril 2015, ÃÂ  11h 41
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
 *
 */

class SequencePhyloLikelihood :
  public AbstractSingleDataPhyloLikelihood
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
  SequencePhyloLikelihood(Context& context, SequenceEvolution& se, size_t nSE = 0, size_t nData = 0) :
    AbstractPhyloLikelihood(context),
    //AbstractAlignedPhyloLikelihood(context, 0),
    AbstractSingleDataPhyloLikelihood(context, 0, (se.getSubstitutionProcessNumbers().size() != 0) ? se.substitutionProcess(se.getSubstitutionProcessNumbers()[0]).getNumberOfStates() : 0, nData),
    seqEvol_(&se),
    nSeqEvol_(nSE)
  {}

  SequencePhyloLikelihood(const SequencePhyloLikelihood& asd) :
    AbstractPhyloLikelihood(asd),
    //AbstractAlignedPhyloLikelihood(asd),
    AbstractSingleDataPhyloLikelihood(asd),
    seqEvol_(asd.seqEvol_),
    nSeqEvol_(asd.nSeqEvol_)
  {}

  virtual ~SequencePhyloLikelihood() {}

  SequencePhyloLikelihood* clone() const override = 0;

public:
  /**
   * @brief The Sequence Evolution
   *
   */
  const SequenceEvolution& sequenceEvolution() const
  {
    return *seqEvol_;
  }

  std::shared_ptr<const SequenceEvolution> getSequenceEvolution() const
  {
    return seqEvol_;
  }

  const size_t getSequenceEvolutionNumber() const
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

class AbstractSequencePhyloLikelihood :
  public SequencePhyloLikelihood,
  public AbstractParametrizable
{
public:
  AbstractSequencePhyloLikelihood(Context& context, SequenceEvolution& se, size_t nSE = 0, size_t nData = 0) :
    AbstractPhyloLikelihood(context),
    //AbstractAlignedPhyloLikelihood(context, 0),
    SequencePhyloLikelihood(context, se, nSE, nData),
    AbstractParametrizable("")
  {
    // initialize INDEPENDENT parameters:
    addParameters_(se.getIndependentParameters());
  }

  AbstractSequencePhyloLikelihood(const AbstractSequencePhyloLikelihood& asd) :
    AbstractPhyloLikelihood(asd),
    //AbstractAlignedPhyloLikelihood(asd),
    SequencePhyloLikelihood(asd),
    AbstractParametrizable(asd)
  {}

  virtual ~AbstractSequencePhyloLikelihood() {}

  AbstractSequencePhyloLikelihood* clone() const override = 0;

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
