// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_SIMULATION_SIMPLESUBSTITUTIONPROCESSSEQUENCESIMULATOR_H
#define BPP_PHYL_SIMULATION_SIMPLESUBSTITUTIONPROCESSSEQUENCESIMULATOR_H


#include "../Likelihood/SequenceEvolution.h"
#include "GivenDataSubstitutionProcessSiteSimulator.h"
#include "SequenceSimulator.h"
#include "SimpleSubstitutionProcessSiteSimulator.h"
#include "SiteSimulator.h"

namespace bpp
{
/**
 * @brief Sequences simulation under a unique substitution process.
 */
class SimpleSubstitutionProcessSequenceSimulator :
  public virtual SequenceSimulatorInterface
{
private:
  std::shared_ptr<SiteSimulatorInterface> siteSim_;

public:
  SimpleSubstitutionProcessSequenceSimulator(std::shared_ptr<const SubstitutionProcessInterface> process) :
    siteSim_(std::make_shared<SimpleSubstitutionProcessSiteSimulator>(process))
  {}

  /**
   * @brief A posterior simulation, from a position in an alignment.
   */
  SimpleSubstitutionProcessSequenceSimulator(
      std::shared_ptr<LikelihoodCalculationSingleProcess> calcul,
      size_t pos) :
    siteSim_(std::make_shared<GivenDataSubstitutionProcessSiteSimulator>(calcul, pos))
  {}

  SimpleSubstitutionProcessSequenceSimulator(std::shared_ptr<SiteSimulatorInterface> simul) :
    siteSim_(simul) {}


  virtual ~SimpleSubstitutionProcessSequenceSimulator() {}

  SimpleSubstitutionProcessSequenceSimulator(const SimpleSubstitutionProcessSequenceSimulator& nhss) :
    siteSim_(nhss.siteSim_)
  {}

  SimpleSubstitutionProcessSequenceSimulator* clone() const override
  {
    return new SimpleSubstitutionProcessSequenceSimulator(*this);
  }

public:
  /**
   * @name The SequenceSimulator interface
   *
   *  @{
   */

  std::unique_ptr<SiteContainerInterface> simulate(size_t numberOfSites) const override;


  const SiteSimulatorInterface& siteSimulator(size_t pos) const override
  {
    return *siteSim_;
  }

  std::vector<std::string> getSequenceNames() const override
  {
    return siteSim_->getSequenceNames();
  }

  /** @} */

  /**
   * @name SiteSimulator and SequenceSimulator interface
   *
   * @{
   */
  std::shared_ptr<const Alphabet> getAlphabet() const override { return siteSim_->getAlphabet(); }

  const Alphabet& alphabet() const override { return siteSim_->alphabet(); }
  /** @} */

  /**
   * @brief Sets whether we will output the internal sequences or not.
   *
   *
   * @param yn Tell if we should output internal sequences.
   */
  void outputInternalSequences(bool yn) override
  {
    siteSim_->outputInternalSites(yn);
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_SIMULATION_SIMPLESUBSTITUTIONPROCESSSEQUENCESIMULATOR_H
