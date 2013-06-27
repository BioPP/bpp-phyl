//
// File: SitePatterns.cpp
// Created by: Julien Dutheil
// Created on: Tue Nov 29 15:37 2005
//  from file PatternTools.cpp
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "SitePatterns.h"

// From the SeqLib library:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

SitePatterns::SitePatterns(const SiteContainer* sequences, bool own) :
  names_(sequences->getSequencesNames()),
  sites_(),
  weights_(),
  indices_(),
  sequences_(sequences),
  alpha_(sequences->getAlphabet()),
  own_(own)
{
  size_t nbSites = sequences->getNumberOfSites();
  vector<SortableSite> ss(nbSites);
  for (size_t i = 0; i < nbSites; i++)
  {
    const Site* currentSite = &sequences->getSite(i);
    SortableSite* ssi = &ss[i];
    ssi->siteS = currentSite->toString();
    ssi->siteP = currentSite;
    ssi->originalPosition = i;
  }

  if (nbSites > 0)
  {
    // Quick sort according to site contents:
    sort(ss.begin(), ss.end());

    // Now build patterns:

    SortableSite* ss0 = &ss[0];
    const Site* previousSite = ss0->siteP;
    indices_.resize(nbSites);
    indices_[ss0->originalPosition] = 0;
    sites_.push_back(previousSite);
    weights_.push_back(1);

    size_t currentPos = 0;
    for (size_t i = 1; i < nbSites; i++)
    {
      SortableSite* ssi = &ss[i];
      const Site* currentSite = ssi->siteP;
      bool siteExists = SiteTools::areSitesIdentical(*currentSite, *previousSite);
      if (siteExists)
      {
        weights_[currentPos]++;
      }
      else
      {
        sites_.push_back(currentSite);
        weights_.push_back(1);
        currentPos++;
      }
      indices_[ssi->originalPosition] = currentPos;
      previousSite = currentSite;
    }
  }
}

/******************************************************************************/

SiteContainer* SitePatterns::getSites() const
{
  SiteContainer* sites = new VectorSiteContainer(sites_, alpha_);
  sites->setSequencesNames(names_, false);
  return sites;
}

/******************************************************************************/

