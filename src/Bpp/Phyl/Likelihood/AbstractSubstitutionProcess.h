// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_ABSTRACTSUBSTITUTIONPROCESS_H
#define BPP_PHYL_LIKELIHOOD_ABSTRACTSUBSTITUTIONPROCESS_H


#include "SubstitutionProcess.h"

// From the STL:
#include <memory>

// From bpp-core:
#include <Bpp/Numeric/AbstractParameterAliasable.h>

namespace bpp
{
/**
 * @brief A partial implementation of the SubstitutionProcess interface.
 *
 * This class handles a pointer toward a ParametrizableTree object, as well
 * as convenient arrays for storing previously computed probabilities.
 */
class AbstractSubstitutionProcess :
  public virtual SubstitutionProcessInterface,
  public virtual AbstractParameterAliasable
{
public:
  AbstractSubstitutionProcess() :
    AbstractParameterAliasable("")
  {}
  
  size_t getNumberOfClasses() const
  {
    auto dist=getRateDistribution();
    return dist?dist->getNumberOfCategories():1;
  }

  size_t getNumberOfStates() const
  {
    return stateMap().getNumberOfModelStates();
  }

  std::shared_ptr<const Alphabet> getAlphabet() const
  {
    return stateMap().getAlphabet();
  }

  bool isCompatibleWith(const AlignmentDataInterface& data) const
  {
    return data.getAlphabet()->getAlphabetType() == getAlphabet()->getAlphabetType();
  }

  /**
   * @brief get NonDerivable parameters
   *
   **/

  ParameterList getNonDerivableParameters() const;

};
} // end namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_ABSTRACTSUBSTITUTIONPROCESS_H
