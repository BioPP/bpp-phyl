// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "GivenDataSubstitutionProcessSequenceSimulator.h"

// From bpp-seq:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

unique_ptr<SiteContainerInterface> GivenDataSubstitutionProcessSequenceSimulator::simulate(
  size_t numberOfSites) const
{
  if (numberOfSites > calcul_->getNumberOfSites())
    throw BadIntegerException("GivenDataSubstitutionProcessSequenceSimulator::simulate. Too many sites to simulate.",(int)numberOfSites);

  auto seqNames = vSiteSim_[0]->getSequenceNames();
  auto sites = make_unique<VectorSiteContainer>(seqNames, getAlphabet());
  sites->setSequenceNames(seqNames, true);

  for (size_t j = 0; j < numberOfSites; j++)
  {
    auto site = vSiteSim_[calcul_->getRootArrayPosition(j)]->simulateSite();
    site->setCoordinate(static_cast<int>(j));
    sites->addSite(site);
  }
  return sites;
}

/******************************************************************************/

