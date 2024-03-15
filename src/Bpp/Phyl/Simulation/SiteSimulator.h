// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_SIMULATION_SITESIMULATOR_H
#define BPP_PHYL_SIMULATION_SITESIMULATOR_H

// From bpp-core
#include <Bpp/Clonable.h>

// From bpp-seq:
#include <Bpp/Seq/Site.h>

namespace bpp
{
/**
 * @brief The SiteSimulator interface.
 * SiteSimulator classes can simulate single sites.
 *
 * @see SequenceSimulator interface for simulating whole sequence sets.
 */
class SiteSimulatorInterface:
  public virtual Clonable
{
public:
  SiteSimulatorInterface() {}
  virtual ~SiteSimulatorInterface() {}

  SiteSimulatorInterface* clone() const override = 0;

public:
  virtual std::unique_ptr<Site> simulateSite() const = 0;

  virtual std::unique_ptr<Site> simulateSite(size_t rateClass) const = 0;

  virtual std::unique_ptr<Site> simulateSite(double rate) const = 0;

  virtual std::unique_ptr<Site> simulateSite(size_t ancestralStateIndex, double rate) const = 0;

  virtual std::vector<std::string> getSequenceNames() const = 0;
  
  virtual std::shared_ptr<const Alphabet> getAlphabet() const = 0;
  
  virtual const Alphabet& alphabet() const = 0;

  virtual void outputInternalSites(bool yn) = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_SIMULATION_SITESIMULATOR_H
