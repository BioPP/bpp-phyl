//
// File: PatternTools.cpp
// Created by: Julien Dutheil
// Created on: Thu Mar 20 13:36:54 2003
//

/*
  Copyright or © or Copr. CNRS, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use, 
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info". 

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability. 

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or 
  data to be ensured and,  more generally, to use and operate it in the 
  same conditions as regards security. 

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#include "PatternTools.h"

#include "Tree/TreeTemplateTools.h"

// From SeqLib:
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <map>
#include <algorithm>

using namespace std;

/******************************************************************************/

SiteContainer* PatternTools::getSequenceSubset(const SiteContainer& sequenceSet, const Node& node) throw (Exception)
{
  size_t nbSites=sequenceSet.getNumberOfSites();
  
  VectorSiteContainer* sequenceSubset = new VectorSiteContainer(sequenceSet.getAlphabet());
  vector<const Node *> leaves = TreeTemplateTools::getLeaves(node);
  for (vector<const Node *>::iterator i = leaves.begin(); i < leaves.end(); i++)
  {
    const Sequence* newSeq = 0;

    try
    {
      newSeq = &sequenceSet.getSequence((*i)->getName());
      sequenceSubset->addSequence(*newSeq);
    }
    catch (std::exception const& e)
    {
      ApplicationTools::displayWarning("Leaf name not found in sequence file: " + (*i)->getName() + " : Replaced with unknown sequence");
      
      BasicSequence seq((*i)->getName(),"",sequenceSet.getAlphabet());
      seq.setToSizeR(nbSites);
      SymbolListTools::changeGapsToUnknownCharacters(seq);
      sequenceSubset->addSequence(seq);
    }
  }

  sequenceSubset->setSitePositions(sequenceSet.getSitePositions());
  
  return sequenceSubset;
}

/******************************************************************************/

SiteContainer* PatternTools::getSequenceSubset(const SiteContainer& sequenceSet, const vector<string>& names) throw (Exception)
{
  VectorSiteContainer* sequenceSubset = new VectorSiteContainer(sequenceSet.getAlphabet());
  for (unsigned int i = 0; i < names.size(); i++)
  {
    const Sequence* newSeq = &sequenceSet.getSequence(names[i]);
    if (!newSeq) throw SequenceNotFoundException("PatternTools ERROR: name not found in sequence file: ", names[i]);
    sequenceSubset->addSequence(*newSeq);
  }
  sequenceSubset->setSitePositions(sequenceSet.getSitePositions());

  return sequenceSubset;
}

/******************************************************************************/

SiteContainer* PatternTools::shrinkSiteSet(const SiteContainer& siteSet) throw (Exception)
{
  if (siteSet.getNumberOfSites() == 0) throw Exception("PatternTools::shrinkSiteSet siteSet is void.");
  vector<const Site *> sites;
  for (unsigned int i = 0; i < siteSet.getNumberOfSites(); i++)
  {
    const Site* currentSite = &siteSet.getSite(i);
    bool siteExists = false;
    for (unsigned int j = 0; !siteExists && j < sites.size(); j++)
    {
      if (SiteTools::areSitesIdentical(* currentSite, * sites[j])) siteExists = true;
    }
    if (!siteExists)	sites.push_back(currentSite);
  }
  SiteContainer* result = new VectorSiteContainer(sites, siteSet.getAlphabet(), false); //We do not check positions here.
  result->setSequencesNames(siteSet.getSequencesNames(), false);
  return result;
}

/******************************************************************************/

Vint PatternTools::getIndexes(const SiteContainer& sequences1, const SiteContainer& sequences2)
{
  size_t nbSites = sequences1.getNumberOfSites(); 
  Vint indexes(nbSites);
  for (size_t i = 0; i < nbSites; i++) {
    //For each site in sequences1,
    indexes[i] = -1;
    const Site* site = &sequences1.getSite(i);
    for (size_t j = 0; j < sequences2.getNumberOfSites(); j++)
    {
      if (SiteTools::areSitesIdentical(*site, sequences2.getSite(j)))
      {
        indexes[i] = static_cast<int>(j);
        break;
      }
    }
  }
  return indexes;
}

/******************************************************************************/

