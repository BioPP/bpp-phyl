// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "ProbabilisticSubstitutionMapping.h"

using namespace bpp;
using namespace std;

void ProbabilisticSubstitutionMapping::setNumberOfSites(size_t numberOfSites)
{
  if (numberOfSites != getNumberOfSites() || (usePatterns_ && numberOfSites != numberOfDistinctSites_))
  {
    AbstractSubstitutionMapping::setNumberOfSites(numberOfSites);

    numberOfDistinctSites_ = numberOfSites;
    usePatterns_ = false;

    unique_ptr<mapTree::EdgeIterator> nIT = allEdgesIterator();

    for ( ; !nIT->end(); nIT->next())
    {
      (**nIT)->setNumberOfSites(numberOfSites);
    }
  }
}

void ProbabilisticSubstitutionMapping::setNumberOfSitesAndTypes(size_t numberOfSites, size_t numberOfTypes)
{
  if (numberOfSites != getNumberOfSites() || (usePatterns_ && numberOfSites != numberOfDistinctSites_))
  {
    numberOfDistinctSites_ = numberOfSites;
    usePatterns_ = false;

    AbstractSubstitutionMapping::setNumberOfSites(numberOfSites);
    AbstractSubstitutionMapping::setNumberOfSubstitutionTypes(numberOfTypes);

    unique_ptr<mapTree::EdgeIterator> nIT = allEdgesIterator();
    for ( ; !nIT->end(); nIT->next())
    {
      (**nIT)->setNumberOfSitesAndTypes(numberOfDistinctSites_, numberOfTypes);
    }
  }
}


void ProbabilisticSubstitutionMapping::setNumberOfSubstitutionTypes(size_t numberOfTypes)
{
  AbstractSubstitutionMapping::setNumberOfSubstitutionTypes(numberOfTypes);

  unique_ptr<mapTree::EdgeIterator> nIT = allEdgesIterator();
  for ( ; !nIT->end(); nIT->next())
  {
    (**nIT)->setNumberOfTypes(numberOfTypes);
  }
}


void ProbabilisticSubstitutionMapping::fillMappingVectorForSite(size_t siteIndex, VVdouble& counts) const
{
  unique_ptr<mapTree::EdgeIterator> nIT = allEdgesIterator();
  size_t nT = (**nIT)->getNumberOfTypes();
  VectorTools::resize2(counts, getNumberOfBranches(), nT);

  auto count_i = counts.begin();
  for ( ; !nIT->end(); nIT->next())
  {
    for (size_t t = 0; t < nT; t++)
    {
      (*count_i)[t] = (***nIT)(siteIndex, t);
    }
    count_i++;
  }
}
