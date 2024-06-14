// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_SIMULATION_SEQUENCESIMULATIONTOOLS_H
#define BPP_PHYL_SIMULATION_SEQUENCESIMULATIONTOOLS_H


#include "SequenceSimulator.h"

// From Seqlib:
#include <Bpp/Seq/Container/SiteContainer.h>

// From the STL:
#include <vector>

namespace bpp
{
/**
 * @brief Tools for sites and sequences simulation.
 */
class SequenceSimulationTools
{
public:
  SequenceSimulationTools() {}
  ~SequenceSimulationTools() {}

public:
  /**
   * @brief Simulate a set of sites knowing their rate.
   *
   * This method is rather slow.
   * consider using a discrete rate distribution and a SequenceSimulator,
   * which is really faster.
   * This method should be used only for continuous rate distribution, or
   * as estimated from posterior rates for instance.
   *
   * @see SequenceSimulator
   * @param simulator A SiteSimulator object to use to simulate sites.
   * @param rates     the rates to use, one for each site to simulate.
   * @return          A container with all simulated sites.
   */
  static std::unique_ptr<SiteContainerInterface> simulateSites(
      const SequenceSimulatorInterface& simulator,
      const std::vector<double>& rates);

  /**
   * @brief Simulate a set of sites knowing their rate and ancestral state.
   *
   * This method is rather slow.
   * consider using a discrete rate distribution and a SequenceSimulator,
   * which is really faster.
   * This method should be used only for continuous rate distribution, or
   * as estimated from posterior rates for instance.
   *
   * @see SequenceSimulator
   * @param simulator A SiteSimulator object to use to simulate sites.
   * @param rates     the rates to use, one for each site to simulate.
   * @param states    the ancestral states to use, one for each site to simulate.
   * @return          A container with all simulated sites.
   */
  static std::unique_ptr<SiteContainerInterface> simulateSites(
      const SequenceSimulatorInterface& simulator,
      const std::vector<double>& rates,
      const std::vector<size_t>& states);

  /**
   * @brief Simulate a set of sites knowing ancestral state.
   *
   * @see SequenceSimulator
   * @param simulator A SiteSimulator object to use to simulate sites.
   * @param states    the ancestral states to use, one for each site to simulate.
   * @return          A container with all simulated sites.
   */
  static std::unique_ptr<SiteContainerInterface> simulateSites(
      const SequenceSimulatorInterface& simulator,
      const std::vector<size_t>& states);
};
} // end of namespace bpp.
#endif // BPP_PHYL_SIMULATION_SEQUENCESIMULATIONTOOLS_H
