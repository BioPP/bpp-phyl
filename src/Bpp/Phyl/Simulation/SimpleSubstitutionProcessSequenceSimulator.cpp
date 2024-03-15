// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SimpleSubstitutionProcessSequenceSimulator.h"

// From boo-seq:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

unique_ptr<SiteContainerInterface> SimpleSubstitutionProcessSequenceSimulator::simulate(
    size_t numberOfSites) const
{
  auto seqNames = siteSim_->getSequenceNames();
  auto sites = make_unique<VectorSiteContainer>(seqNames, getAlphabet());
  sites->setSequenceNames(seqNames, true);
  for (size_t j = 0; j < numberOfSites; ++j)
  {
    auto site = siteSim_->simulateSite();
    site->setCoordinate(static_cast<int>(j));
    sites->addSite(site);
  }
  return sites;
}

/******************************************************************************/
