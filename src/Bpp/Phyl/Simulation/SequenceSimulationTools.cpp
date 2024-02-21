//
// File: SequenceSimulationTools.cpp
// Authors:
//   Julien Dutheil
// Created: 2005-08-24 16:25:00
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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
