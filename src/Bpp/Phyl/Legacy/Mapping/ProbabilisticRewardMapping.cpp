// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "ProbabilisticRewardMapping.h"

using namespace bpp;

void LegacyProbabilisticRewardMapping::setTree(const Tree& tree)
{
  LegacyAbstractRewardMapping::setTree(tree);
  for (size_t i = 0; i < getNumberOfSites(); i++) 
    mapping_[i].resize(getNumberOfBranches());
}

void LegacyProbabilisticRewardMapping::setNumberOfSites(size_t numberOfSites)
{
  LegacyAbstractRewardMapping::setNumberOfSites(numberOfSites);
  mapping_.resize(numberOfSites);
  for (size_t i = 0; i < numberOfSites; i++) {
    mapping_[i].resize(getNumberOfBranches());
  }
}
 
