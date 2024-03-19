// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_SIMULATION_SEQUENCESIMULATOR_H
#define BPP_PHYL_SIMULATION_SEQUENCESIMULATOR_H


// From bpp-core:
#include <Bpp/Clonable.h>

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainer.h>

#include "SiteSimulator.h"

namespace bpp
{
/**
 * @brief The SequenceSimulator interface.
 * SequenceSimulator classes can simulate whole datasets.
 */
class SequenceSimulatorInterface :
  public virtual Clonable
{
public:
  SequenceSimulatorInterface() {}
  virtual ~SequenceSimulatorInterface() {}

  SequenceSimulatorInterface* clone() const override = 0;

public:
  virtual std::unique_ptr<SiteContainerInterface> simulate(size_t numberOfSites) const = 0;

  virtual const SiteSimulatorInterface& siteSimulator(size_t pos) const = 0;

  virtual std::vector<std::string> getSequenceNames() const = 0;

  virtual std::shared_ptr<const Alphabet> getAlphabet() const = 0;

  virtual const Alphabet& alphabet() const = 0;

  virtual void outputInternalSequences(bool inter) = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_SIMULATION_SEQUENCESIMULATOR_H
