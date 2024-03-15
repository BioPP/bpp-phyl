// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "StateMap.h"

using namespace bpp;
using namespace std;

CanonicalStateMap::CanonicalStateMap(std::shared_ptr<const Alphabet> alphabet, bool includeGaps) :
  AbstractStateMap(alphabet)
{
  size_t i = 0;
  while (states_.size() < alphabet->getSize())
  {
    if (!alphabet->isGap(alphabet->getIntCodeAt(i)))
      states_.push_back(alphabet->getIntCodeAt(i));
    i++;
  }

  if (includeGaps)
    states_.push_back(alphabet->getGapCharacterCode());
}

CanonicalStateMap::CanonicalStateMap(const StateMapInterface& sm, bool includeGaps) :
  AbstractStateMap(sm.getAlphabet())
{
  for (size_t i = 0; i < sm.getNumberOfModelStates(); ++i)
  {
    states_.push_back(sm.getAlphabetStateAsInt(i));
  }
  if (includeGaps)
    states_.push_back(sm.getAlphabet()->getGapCharacterCode());
}

MarkovModulatedStateMap::MarkovModulatedStateMap(const StateMapInterface& unitMap, unsigned int nbClasses) :
  AbstractStateMap(unitMap.getAlphabet()),
  nbClasses_(nbClasses)
{
  for (unsigned int j = 0; j < nbClasses; ++j)
  {
    for (size_t i = 0; i < unitMap.getNumberOfModelStates(); ++i)
    {
      states_.push_back(unitMap.getAlphabetStateAsInt(i));
    }
  }
}
