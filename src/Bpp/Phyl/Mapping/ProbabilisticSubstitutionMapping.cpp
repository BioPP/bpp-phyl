//
// File: ProbabilisticSubstitutionMapping.cpp
// Authors:
//   Julien Dutheil
// Created: 2006-04-05 10:47:00
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/


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
