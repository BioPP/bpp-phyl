// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_SIMULATION_GIVENDATASUBSTITUTIONPROCESSSEQUENCESIMULATOR_H
#define BPP_PHYL_SIMULATION_GIVENDATASUBSTITUTIONPROCESSSEQUENCESIMULATOR_H


#include "GivenDataSubstitutionProcessSiteSimulator.h"
#include "SequenceSimulator.h"

namespace bpp
{
/**
 * @brief Sequences simulation under a unique substitution process, but with site specific
 * posterior probabilities.
 */
class GivenDataSubstitutionProcessSequenceSimulator :
  public virtual SequenceSimulatorInterface
{
private:
  std::shared_ptr<LikelihoodCalculationSingleProcess> calcul_;

  /*
   * @brief Vector of site specific site simulators on SHRUNKED data
   *
   * More efficient implementation is possible (all in ona).
   *
   */

  std::vector<std::shared_ptr<GivenDataSubstitutionProcessSiteSimulator>> vSiteSim_;

  /*
   * @brief Vector of branch indexes where the data is not used, and
   * substitution probabilities are computed from the prior process.
   *
   */

  std::vector<uint> vPriorBranch_;

public:
  GivenDataSubstitutionProcessSequenceSimulator(std::shared_ptr<LikelihoodCalculationSingleProcess> calcul, std::vector<uint> vPrior = std::vector<uint>()) :
    calcul_(calcul), vSiteSim_(), vPriorBranch_(vPrior)
  {
    for (size_t i = 0; i < calcul_->getNumberOfDistinctSites(); i++)
    {
      vSiteSim_.push_back(std::make_shared<GivenDataSubstitutionProcessSiteSimulator>(calcul_, i, true, vPriorBranch_));
    }
  }

  virtual ~GivenDataSubstitutionProcessSequenceSimulator()
  {}

  GivenDataSubstitutionProcessSequenceSimulator(const GivenDataSubstitutionProcessSequenceSimulator& nhss) :
    calcul_(nhss.calcul_), vSiteSim_(nhss.vSiteSim_)
  {}

  GivenDataSubstitutionProcessSequenceSimulator* clone() const override
  {
    return new GivenDataSubstitutionProcessSequenceSimulator(*this);
  }

public:
  /**
   * @name The SequenceSimulator interface
   * Here the numberOfSites is unused (awkward inheritance...)
   *
   */ 
  std::unique_ptr<SiteContainerInterface> simulate() const
  {
    return simulate(calcul_->getNumberOfSites());
  }

  std::unique_ptr<SiteContainerInterface> simulate(size_t numberOfSites) const override;

  const SiteSimulatorInterface& siteSimulator(size_t pos) const override
  {
    
    return *vSiteSim_[calcul_->getRootArrayPosition(pos)];
  }

  /**
   * @name SiteSimulator and SequenceSimulator interface
   *
   * @{
   */
  std::shared_ptr<const Alphabet> getAlphabet() const override
  {
    return vSiteSim_[0]->getAlphabet();
  }

  const Alphabet& alphabet() const override
  {
    return vSiteSim_[0]->alphabet();
  }

  std::vector<std::string> getSequenceNames() const override
  {
    return vSiteSim_[0]->getSequenceNames();
  }

  /**
   * @brief the number of mapped sites.
   */
  size_t getNumberOfSites() const
  {
    return calcul_->getNumberOfSites();
  }

  /** @} */

  /**
   * @brief Sets whether we will output the internal sequences or not.
   *
   *
   * @param yn Tell if we should output internal sequences.
   */
  void outputInternalSequences(bool yn) override
  {
    for (auto& siteSim : vSiteSim_)
    {
      siteSim->outputInternalSites(yn);
    }
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_SIMULATION_GIVENDATASUBSTITUTIONPROCESSSEQUENCESIMULATOR_H
