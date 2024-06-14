// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "ProbabilisticRewardMapping.h"

using namespace bpp;
using namespace std;

void ProbabilisticRewardMapping::setNumberOfSites(size_t numberOfSites)
{
  if (numberOfSites != getNumberOfSites() || (usePatterns_ && numberOfSites != numberOfDistinctSites_))
  {
    AbstractRewardMapping::setNumberOfSites(numberOfSites);

    numberOfDistinctSites_ = numberOfSites;
    usePatterns_ = false;

    unique_ptr<mapTree::EdgeIterator> nIT = allEdgesIterator();

    for ( ; !nIT->end(); nIT->next())
    {
      (**nIT)->setNumberOfSites(numberOfSites);
    }
  }
}
