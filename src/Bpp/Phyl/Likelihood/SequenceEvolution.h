// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_SEQUENCEEVOLUTION_H
#define BPP_PHYL_LIKELIHOOD_SEQUENCEEVOLUTION_H


#include "SubstitutionProcess.h"

// From bpp-core:
#include <Bpp/Numeric/ParameterAliasable.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>

// From bpp-seq:

#include <Bpp/Seq/Container/AlignmentData.h>

// From the STL:
#include <memory>

namespace bpp
{
/**
 * @brief This interface describes the evolution process of a sequence.
 *
 * It main purpose is to provide the necessary calculus for each
 * management of site specific substitution process, such as partition
 * models and HMM.
 *
 * This object has the INDEPENDENT parameters of the processes.
 */
class SequenceEvolution :
  public virtual ParameterAliasable
{
public:
  virtual SequenceEvolution* clone() const override = 0;

public:
  virtual bool isCompatibleWith(const AlignmentDataInterface& data) const = 0;

  virtual const std::vector<size_t>& getSubstitutionProcessNumbers() const = 0;

  virtual const SubstitutionProcessInterface& substitutionProcess(size_t number) const = 0;

  virtual std::shared_ptr<const SubstitutionProcessInterface> getSubstitutionProcess(size_t number) const = 0;

  const StateMapInterface& stateMap() const
  {
    return substitutionProcess(getSubstitutionProcessNumbers()[0]).stateMap();
  }

  std::shared_ptr<const StateMapInterface> getStateMap() const
  {
    return substitutionProcess(getSubstitutionProcessNumbers()[0]).getStateMap();
  }

  /**
   * @brief Get the branch lengths parameters.
   *
   * @return A ParameterList with all branch lengths.
   */
  virtual ParameterList getBranchLengthParameters(bool independent) const = 0;

  /**
   * @brief Get the parameters associated to substitution model(s).
   *
   * @return A ParameterList.
   */
  virtual ParameterList getSubstitutionModelParameters(bool independent) const = 0;

  /**
   * @brief Get the parameters associated to substitution processes(s).
   *
   * @return A ParameterList.
   */
  virtual ParameterList getSubstitutionProcessParameters(bool independent) const = 0;

  /**
   * @brief Get the parameters associated to the rate distribution(s).
   *
   * @return A ParameterList.
   */
  virtual ParameterList getRateDistributionParameters(bool independent) const = 0;

  /**
   * @brief Get the parameters associated to the root frequencies(s).
   *
   * @return A ParameterList.
   */
  virtual ParameterList getRootFrequenciesParameters(bool independent) const = 0;


  /**
   * @brief All non derivable parameters.
   *
   * Usually, this contains all substitution model parameters and
   * rate distribution, and alias.
   *
   * @return A ParameterList.
   */
  virtual ParameterList getNonDerivableParameters() const = 0;
};
} // end namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_SEQUENCEEVOLUTION_H
