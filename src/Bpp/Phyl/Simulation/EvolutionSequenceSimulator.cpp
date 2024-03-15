// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Random/RandomTools.h>

#include "../Likelihood/AutoCorrelationSequenceEvolution.h"
#include "../Likelihood/HmmSequenceEvolution.h"
#include "../Likelihood/MixtureSequenceEvolution.h"
#include "../Likelihood/OneProcessSequenceEvolution.h"
#include "../Likelihood/PartitionSequenceEvolution.h"
#include "EvolutionSequenceSimulator.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

EvolutionSequenceSimulator::EvolutionSequenceSimulator(const SequenceEvolution& evol) :
  SubstitutionProcessSequenceSimulator(evol),
  seqEvol_(&evol)
{}

/******************************************************************************/

void EvolutionSequenceSimulator::resetSiteSimulators(size_t numberOfSites) const
{
  vMap_.resize(numberOfSites);

  const OneProcessSequenceEvolution* opse = dynamic_cast<const OneProcessSequenceEvolution*>(seqEvol_);

  if (opse)
  {
    size_t nProc = opse->getSubstitutionProcessNumbers()[0];
    for (size_t i = 0; i < numberOfSites; i++)
    {
      vMap_[i] = nProc;
    }
  }
  else
  {
    const PartitionSequenceEvolution* pse = dynamic_cast<const PartitionSequenceEvolution*>(seqEvol_);

    if (pse)
    {
      if (numberOfSites > pse->getNumberOfSites())
        throw IndexOutOfBoundsException("EvolutionSequenceSimulator::resetSiteSimulators", numberOfSites, 0, pse->getNumberOfSites());

      for (size_t i = 0; i < numberOfSites; i++)
      {
        vMap_[i] = pse->getSubstitutionProcessNumber(i);
      }
    }

    else
    {
      const vector<size_t> nProc = seqEvol_->getSubstitutionProcessNumbers();

      const MixtureSequenceEvolution* mse = dynamic_cast<const MixtureSequenceEvolution*>(seqEvol_);

      if (mse)
      {
        const vector<double>& vprob = mse->getSubProcessProbabilities();

        RandomTools::getSample(nProc, vprob, vMap_, true);
      }
      else
      {
        const AbstractHmmTransitionMatrix* htm = 0;

        const AutoCorrelationSequenceEvolution* ase = dynamic_cast<const AutoCorrelationSequenceEvolution*>(seqEvol_);

        if (ase)
          htm = &ase->hmmTransitionMatrix();

        const HmmSequenceEvolution* hse = dynamic_cast<const HmmSequenceEvolution*>(seqEvol_);

        if (hse)
          htm = &hse->hmmTransitionMatrix();
        if (htm)
        {
          vector<size_t> vInd = htm->sample(numberOfSites);
          for (size_t i = 0; i < numberOfSites; i++)
          {
            vMap_[i] = nProc[vInd[i]];
          }
        }
        else
          throw Exception("EvolutionSequenceSimulator::resetSiteSimulators : unknow Sequence Evolution.");
      }
    }
  }
}
