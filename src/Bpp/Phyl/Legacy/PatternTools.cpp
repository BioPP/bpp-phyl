//
// File: PatternTools.cpp
// Authors:
//   Julien Dutheil
// Created: 2003-03-20 13:36:54
//

/*
  Copyright or ÃÂ© or Copr. CNRS, (November 16, 2004)
  
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


#include "../Tree/TreeTemplateTools.h"
#include "PatternTools.h"

using namespace bpp;

// From the STL:
#include <iostream>
#include <algorithm>

using namespace std;


/******************************************************************************/

AlignedValuesContainer* PatternToolsOld::getSequenceSubset(const AlignedValuesContainer& sequenceSet, const vector<string>& names)
{
  AlignedValuesContainer* result;

  if (dynamic_cast<const SiteContainer*>(&sequenceSet))
  {
    const SiteContainer& sitecontainer = dynamic_cast<const SiteContainer&>(sequenceSet);

    VectorSiteContainer* sequenceSubset = new VectorSiteContainer(sequenceSet.getAlphabet());
    result = sequenceSubset;

    for (unsigned int i = 0; i < names.size(); i++)
    {
      const Sequence* newSeq = &sitecontainer.getSequence(names[i]);
      if (!newSeq)
        throw SequenceNotFoundException("PatternToolsOld ERROR: name not found in sequence file: ", names[i]);
      sequenceSubset->addSequence(*newSeq);
    }
  }
  else if (dynamic_cast<const VectorProbabilisticSiteContainer*>(&sequenceSet))
  {
    const VectorProbabilisticSiteContainer& sitecontainer = dynamic_cast<const VectorProbabilisticSiteContainer&>(sequenceSet);

    VectorProbabilisticSiteContainer* sequenceSubset = new VectorProbabilisticSiteContainer(sequenceSet.getAlphabet());
    result = sequenceSubset;

    for (unsigned int i = 0; i < names.size(); i++)
    {
      shared_ptr<BasicProbabilisticSequence> newSeq = sitecontainer.getSequence(names[i]);
      if (!newSeq)
        throw SequenceNotFoundException("PatternToolsOld ERROR: name not found in sequence file: ", names[i]);
      sequenceSubset->addSequence(newSeq);
    }
  }
  else
    throw Exception("PatternToolsOld::getSequenceSubset : this should not happen.");

  result->setSitePositions(sequenceSet.getSitePositions());
  return result;
}

/******************************************************************************/

AlignedValuesContainer* PatternToolsOld::shrinkSiteSet(const AlignedValuesContainer& siteSet)
{
  if (siteSet.getNumberOfSites() == 0)
    throw Exception("PatternToolsOld::shrinkSiteSet siteSet is void.");
  AlignedValuesContainer* result;

  const SiteContainer* sc = dynamic_cast<const SiteContainer*>(&siteSet);
  const VectorProbabilisticSiteContainer* psc = dynamic_cast<const VectorProbabilisticSiteContainer*>(&siteSet);

  if (!sc && !psc)
    throw Exception("PatternToolsOld::shrinkSiteSet : this should not happen.");

  vector<const CruxSymbolListSite*> sites;
  for (unsigned int i = 0; i < siteSet.getNumberOfSites(); i++)
  {
    const CruxSymbolListSite* currentSite = sc ?
                                            dynamic_cast<const CruxSymbolListSite*>(&sc->getSite(i)) :
                                            dynamic_cast<const CruxSymbolListSite*>(psc->getSite(i).get());

    bool siteExists = false;
    for (unsigned int j = 0; !siteExists && j < sites.size(); j++)
    {
      if (SymbolListTools::areSymbolListsIdentical(*currentSite, *sites[j]))
        siteExists = true;
    }
    if (!siteExists)
      sites.push_back(currentSite);
  }
  result = sc ? dynamic_cast<AlignedValuesContainer*>(new VectorSiteContainer(sites, siteSet.getAlphabet(), false)) :
           dynamic_cast<AlignedValuesContainer*>(new VectorProbabilisticSiteContainer(sites, siteSet.getAlphabet(), false)); // We do not check positions here.
  result->setSequencesNames(siteSet.getSequencesNames(), false);
  return result;
}

/******************************************************************************/

Vint PatternToolsOld::getIndexes(const AlignedValuesContainer& sequences1, const AlignedValuesContainer& sequences2)
{
  size_t nbSites = sequences1.getNumberOfSites();
  Vint indexes(nbSites);

  const SiteContainer* sc1 = dynamic_cast<const SiteContainer*>(&sequences1);
  const VectorProbabilisticSiteContainer* psc1 = dynamic_cast<const VectorProbabilisticSiteContainer*>(&sequences1);
  const SiteContainer* sc2 = dynamic_cast<const SiteContainer*>(&sequences2);
  const VectorProbabilisticSiteContainer* psc2 = dynamic_cast<const VectorProbabilisticSiteContainer*>(&sequences2);

  if (!(sc1 && sc2) && !(psc1 && psc2))
    throw Exception("PatternToolsOld::getIndexes : this should not happen.");

  for (size_t i = 0; i < nbSites; i++)
  {
    // For each site in sequences1,
    indexes[i] = -1;
    const CruxSymbolListSite* site1 = sc1 ?
                                      dynamic_cast<const CruxSymbolListSite*>(&sc1->getSite(i)) :
                                      dynamic_cast<const CruxSymbolListSite*>(psc1->getSite(i).get());

    for (size_t j = 0; j < sequences2.getNumberOfSites(); j++)
    {
      const CruxSymbolListSite* site2 = sc2 ?
                                        dynamic_cast<const CruxSymbolListSite*>(&sc2->getSite(i)) :
                                        dynamic_cast<const CruxSymbolListSite*>(psc2->getSite(i).get());

      if (SymbolListTools::areSymbolListsIdentical(*site1, *site2))
      {
        indexes[i] = static_cast<int>(j);
        break;
      }
    }
  }
  return indexes;
}

/******************************************************************************/
