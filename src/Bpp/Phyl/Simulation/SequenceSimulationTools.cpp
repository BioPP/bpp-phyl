// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "GivenDataSubstitutionProcessSiteSimulator.h"
#include "SequenceSimulationTools.h"
#include "SequenceSimulator.h"
#include "SimpleSubstitutionProcessSiteSimulator.h"
#include "SiteSimulator.h"

// From bpp-seq:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;
using namespace std;

unique_ptr<SiteContainerInterface> SequenceSimulationTools::simulateSites(
    const SequenceSimulatorInterface& simulator,
    const vector<double>& rates)
{
  size_t numberOfSites = rates.size();
  vector<unique_ptr<Site>> vs(numberOfSites);
  for (size_t i = 0; i < numberOfSites; ++i)
  {
    auto s = simulator.siteSimulator(i).simulateSite(rates[i]);
    s->setCoordinate(static_cast<int>(i));
    vs[i] = std::move(s);
  }
  auto sites = make_unique<VectorSiteContainer>(vs, simulator.getAlphabet());
  sites->setSequenceNames(simulator.getSequenceNames(), true);

  return sites;
}

unique_ptr<SiteContainerInterface> SequenceSimulationTools::simulateSites(
    const SequenceSimulatorInterface& simulator,
    const vector<double>& rates,
    const vector<size_t>& states)
{
  size_t numberOfSites = rates.size();
  if (states.size() != numberOfSites)
    throw Exception("SequenceSimulationTools::simulateSites., 'rates' and 'states' must have the same length.");
  vector<unique_ptr<Site>> vs(numberOfSites);
  for (size_t i = 0; i < numberOfSites; ++i)
  {
    auto s = simulator.siteSimulator(i).simulateSite(states[i], rates[i]);
    s->setCoordinate(static_cast<int>(i));
    vs[i] = std::move(s);
  }

  auto sites = make_unique<VectorSiteContainer>(vs, simulator.getAlphabet());
  sites->setSequenceNames(simulator.getSequenceNames(), true);

  return sites;
}

unique_ptr<SiteContainerInterface> SequenceSimulationTools::simulateSites(
    const SequenceSimulatorInterface& simulator,
    const vector<size_t>& states)
{
  size_t numberOfSites = states.size();
  vector<unique_ptr<Site>> vs(numberOfSites);
  for (size_t i = 0; i < numberOfSites; ++i)
  {
    auto s = simulator.siteSimulator(i).simulateSite(states[i], 1);
    s->setCoordinate(static_cast<int>(i));
    vs[i] = std::move(s);
  }

  auto sites = make_unique<VectorSiteContainer>(vs, simulator.getAlphabet());
  sites->setSequenceNames(simulator.getSequenceNames(), true);

  return sites;
}
