//
// File: SequenceSimulationTools.cpp
// Created by: Julien Dutheil
// Created on: Wed Aug  24 16:25 2005
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

#include "SequenceSimulationTools.h"

// From bpp-seq:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;
using namespace std;

SiteContainer* SequenceSimulationTools::simulateSites(const SiteSimulator& simulator, const vector<double>& rates)
{
  size_t numberOfSites = rates.size();
  vector<const Site*> vs(numberOfSites);
  for (size_t i = 0; i < numberOfSites; i++)
  {
    Site* s = simulator.simulateSite(rates[i]);
    s->setPosition(static_cast<int>(i));
    vs[i] = s;
  }
  SiteContainer* sites = new VectorSiteContainer(vs, simulator.getAlphabet());
  sites->setSequencesNames(simulator.getSequencesNames(), false);
  // Freeing memory:
  for (size_t i = 0; i < numberOfSites; i++)
  {
    delete vs[i];
  }

  return sites;
}

SiteContainer* SequenceSimulationTools::simulateSites(const SiteSimulator& simulator, const vector<double>& rates, const vector<size_t>& states)
throw (Exception)
{
  size_t numberOfSites = rates.size();
  if (states.size() != numberOfSites)
    throw Exception("SequenceSimulationTools::simulateSites., 'rates' and 'states' must have the same length.");
  vector<const Site*> vs(numberOfSites);
  for (size_t i = 0; i < numberOfSites; i++)
  {
    Site* s = simulator.simulateSite(states[i], rates[i]);
    s->setPosition(static_cast<int>(i));
    vs[i] = s;
  }
  SiteContainer* sites = new VectorSiteContainer(vs, simulator.getAlphabet());
  sites->setSequencesNames(simulator.getSequencesNames(), false);
  // Freeing memory:
  for (size_t i = 0; i < numberOfSites; i++)
  {
    delete vs[i];
  }

  return sites;
}

SiteContainer* SequenceSimulationTools::simulateSites(const SiteSimulator& simulator, const vector<size_t>& states)
throw (Exception)
{
  size_t numberOfSites = states.size();
  vector<const Site*> vs(numberOfSites);
  for (size_t i = 0; i < numberOfSites; i++)
  {
    Site* s = simulator.simulateSite(states[i]);
    s->setPosition(static_cast<int>(i));
    vs[i] = s;
  }
  SiteContainer* sites = new VectorSiteContainer(vs, simulator.getAlphabet());
  sites->setSequencesNames(simulator.getSequencesNames(), false);
  // Freeing memory:
  for (size_t i = 0; i < numberOfSites; i++)
  {
    delete vs[i];
  }

  return sites;
}

