// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "ProbabilisticSubstitutionMapping.h"

using namespace bpp;

void LegacyProbabilisticSubstitutionMapping::setTree(const Tree& tree)
{
  LegacyAbstractSubstitutionMapping::setTree(tree);
  for (size_t i = 0; i < getNumberOfSites(); i++) {
    mapping_[i].resize(getNumberOfBranches());
    for (size_t j = 0; j < getNumberOfBranches(); j++) {
      mapping_[i][j].resize(1);
    }
  }
}

void LegacyProbabilisticSubstitutionMapping::setNumberOfSites(size_t numberOfSites)
{
  LegacyAbstractSubstitutionMapping::setNumberOfSites(numberOfSites);
  mapping_.resize(numberOfSites);
  for (size_t i = 0; i < numberOfSites; i++) {
    mapping_[i].resize(getNumberOfBranches());
    for (size_t j = 0; j < getNumberOfBranches(); j++) {
      mapping_[i][j].resize(getNumberOfSubstitutionTypes());
    }
  }
}
 
