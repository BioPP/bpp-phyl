// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_SIMULATION_EVOLUTIONSEQUENCESIMULATOR_H
#define BPP_PHYL_SIMULATION_EVOLUTIONSEQUENCESIMULATOR_H


#include "../Likelihood/SequenceEvolution.h"
#include "SubstitutionProcessSequenceSimulator.h"

namespace bpp
{
class EvolutionSequenceSimulator :
  public SubstitutionProcessSequenceSimulator
{
private:
  const SequenceEvolution* seqEvol_;

public:
  EvolutionSequenceSimulator(const SequenceEvolution& evol);

  virtual ~EvolutionSequenceSimulator()
  {}

  EvolutionSequenceSimulator(const EvolutionSequenceSimulator& ess) :
    SubstitutionProcessSequenceSimulator(ess),
    seqEvol_(ess.seqEvol_)
  {}

  EvolutionSequenceSimulator& operator=(const EvolutionSequenceSimulator& ess)
  {
    SubstitutionProcessSequenceSimulator::operator=(ess);
    seqEvol_ = ess.seqEvol_;

    return *this;
  }

  EvolutionSequenceSimulator* clone() const { return new EvolutionSequenceSimulator(*this); }

public:
  /**
   * @brief Reset the succession of site simulators (useful in case
   * of HMM).
   *
   **/

  void resetSiteSimulators(size_t numberOfSites) const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_SIMULATION_EVOLUTIONSEQUENCESIMULATOR_H
