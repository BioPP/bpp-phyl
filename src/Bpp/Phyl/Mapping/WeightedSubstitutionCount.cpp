// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "WeightedSubstitutionCount.h"

using namespace bpp;

void AbstractWeightedSubstitutionCount::setWeights(std::shared_ptr<const AlphabetIndex2> weights)
{
  weights_ = weights;
  weightsHaveChanged();
}
