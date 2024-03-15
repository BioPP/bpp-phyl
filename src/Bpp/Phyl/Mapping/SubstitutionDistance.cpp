// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SubstitutionDistance.h"

using namespace bpp;

void AbstractSubstitutionDistance::setDistances(std::shared_ptr<const AlphabetIndex2> distances)
{
  distances_ = distances;
  distancesHaveChanged();
}
