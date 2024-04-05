// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_SIMULATION_SUBSTITUTIONPROCESSSEQUENCESIMULATOR_H
#define BPP_PHYL_SIMULATION_SUBSTITUTIONPROCESSSEQUENCESIMULATOR_H

#include <Bpp/Phyl/Likelihood/SequenceEvolution.h>

#include "SequenceSimulator.h"
#include "SiteSimulator.h"

namespace bpp
{
/**
 * @brief Sequences simulation under position specific substitution process.
 *
 */

class SubstitutionProcessSequenceSimulator :
  public virtual SequenceSimulatorInterface
{
protected:
  /**
   * @brief the map of the process simulators.
   *
   */
  std::map<size_t, std::shared_ptr<SiteSimulatorInterface>> mProcess_;

  /**
   * @brief The vector of the site specific process in mProcess_;
   * is mutable because can be changed for each simulation (for ex
   * in case of HMM).
   */
  mutable std::vector<size_t> vMap_;

  /**
   * @brief all processes trees must have at least the same sequence
   * names as the first process of the map.
   */
  std::vector<std::string> seqNames_;

  /**
   * @brief correspondence map of seqNames positions of the several trees.
   * Reference is the tree of the first process of the map.
   *
   * mvPosNames[process id][i] is the position in the id_th tree
   * leaves names of the i_th name of seqName_.
   */
  std::map<size_t, std::vector<size_t>> mvPosNames_;

public:
  SubstitutionProcessSequenceSimulator(const SequenceEvolution& evol);

  SubstitutionProcessSequenceSimulator(const SubstitutionProcessSequenceSimulator&);

  SubstitutionProcessSequenceSimulator& operator=(const SubstitutionProcessSequenceSimulator&);

  SubstitutionProcessSequenceSimulator* clone() const override { return new SubstitutionProcessSequenceSimulator(*this); }

  virtual ~SubstitutionProcessSequenceSimulator() {}

  const SiteSimulatorInterface& siteSimulator(size_t pos) const override
  {
    if (pos > vMap_.size())
      throw BadIntegerException("Out of range position for SubstitutionProcessSequenceSimulator", (int)pos);
    return *mProcess_.at(vMap_[pos]);
  }

  std::vector<std::string> getSequenceNames() const override
  {
    return seqNames_;
  }

  /**
   * @brief Sets whether we will output the internal sequences or not.
   *
   *
   * @param yn Tell if we should output internal sequences.
   */
  void outputInternalSequences(bool yn) override;

  /**
   * @brief reset the set of processes.
   */
  void setMap(std::vector<size_t> vMap);

  /**
   * @brief the number of mapped sites.
   */
  size_t getNumberOfSites() const
  {
    return vMap_.size();
  }

  std::unique_ptr<SiteContainerInterface> simulate(size_t numberOfSites) const override;

  std::unique_ptr<SiteContainerInterface> simulate(const std::vector<double>& rates) const;

  std::unique_ptr<SiteContainerInterface> simulate(const std::vector<size_t>& states) const;

  std::unique_ptr<SiteContainerInterface> simulate(const std::vector<double>& rates, const std::vector<size_t>& states) const;

  std::shared_ptr<const Alphabet> getAlphabet() const override;

  const Alphabet& alphabet() const override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_SIMULATION_SUBSTITUTIONPROCESSSEQUENCESIMULATOR_H
