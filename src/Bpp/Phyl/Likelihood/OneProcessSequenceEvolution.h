// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_ONEPROCESSSEQUENCEEVOLUTION_H
#define BPP_PHYL_LIKELIHOOD_ONEPROCESSSEQUENCEEVOLUTION_H


// From bpp-core:
#include <Bpp/Numeric/ParameterAliasable.h>

// From the STL:
#include <memory>

#include "SequenceEvolution.h"
#include "SubstitutionProcess.h"

namespace bpp
{
/**
 * @brief Evolution of a sequence performed by a unique
 * SubstitutionProcess all along the sequence.
 *
 */

class OneProcessSequenceEvolution :
  public virtual SequenceEvolution,
  public AbstractParameterAliasable
{
protected:
  std::shared_ptr<SubstitutionProcessInterface> subsProc_;

  /**
   * @brief the substitution process number.
   */
  size_t nProc_;

  /**
   * @brief not nice, for inheritance compatibility.
   */
  std::vector<size_t> vProc_;

public:
  OneProcessSequenceEvolution(std::shared_ptr<SubstitutionProcessInterface> process, size_t nProc = 0);

  OneProcessSequenceEvolution(const OneProcessSequenceEvolution& evol);

  OneProcessSequenceEvolution& operator=(const OneProcessSequenceEvolution& evol);

  OneProcessSequenceEvolution* clone() const override
  {
    return new OneProcessSequenceEvolution(*this);
  }

  bool isCompatibleWith(const AlignmentDataInterface& data) const override
  {
    return subsProc_->isCompatibleWith(data);
  }

  const size_t getSubstitutionProcessNumber() const
  {
    return nProc_;
  }

  const SubstitutionProcessInterface& substitutionProcess() const
  {
    return *subsProc_;
  }

  std::shared_ptr<const SubstitutionProcessInterface> getSubstitutionProcess() const
  {
    return subsProc_;
  }

  SubstitutionProcessInterface& substitutionProcess()
  {
    return *subsProc_;
  }

  std::shared_ptr<SubstitutionProcessInterface> getSubstitutionProcess()
  {
    return subsProc_;
  }

  const std::vector<size_t>& getSubstitutionProcessNumbers() const override
  {
    return vProc_;
  }

  const SubstitutionProcessInterface& substitutionProcess(size_t number) const override
  {
    return *subsProc_;
  }

  std::shared_ptr<const SubstitutionProcessInterface> getSubstitutionProcess(size_t number) const override
  {
    return subsProc_;
  }

  /**
   * @brief Get the tree (topology and branch lengths).
   *
   * @return The tree of this OneProcessSequenceEvolution object.
   */
  std::shared_ptr<const ParametrizablePhyloTree> tree() const
  {
    return subsProc_->getParametrizablePhyloTree();
  }

  void fireParameterChanged(const ParameterList& pl) override
  {
    // Updates substitution process:
    subsProc_->matchParametersValues(pl);
  }


  /**
   * @brief Get several categories of parameters
   *
   **/
  ParameterList getBranchLengthParameters(bool independent) const override
  {
    return subsProc_->getBranchLengthParameters(independent);
  }

  ParameterList getRootFrequenciesParameters(bool independent) const override
  {
    return subsProc_->getRootFrequenciesParameters(independent);
  }

  ParameterList getRateDistributionParameters(bool independent) const override
  {
    return subsProc_->getRateDistributionParameters(independent);
  }

  ParameterList getSubstitutionModelParameters(bool independent) const override
  {
    return subsProc_->getSubstitutionModelParameters(independent);
  }

  ParameterList getSubstitutionProcessParameters(bool independent) const override
  {
    if (independent)
      return subsProc_->getIndependentParameters();
    else
      return subsProc_->getParameters();
  }

  const ParameterList& getParameters() const override
  {
    return subsProc_->getParameters();
  }

  const ParameterList& getIndependentParameters() const override
  {
    return subsProc_->getIndependentParameters();
  }

  /**
   * @brief get (Non)Derivable INDEPENDENT parameters
   *
   **/
  ParameterList getNonDerivableParameters() const override
  {
    return subsProc_->getNonDerivableParameters();
  }
};
} // end namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_ONEPROCESSSEQUENCEEVOLUTION_H
