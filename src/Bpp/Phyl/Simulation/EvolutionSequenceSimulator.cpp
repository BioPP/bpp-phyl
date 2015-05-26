//
// File: EvolutionSequenceSimulator.cpp
// Created by: Laurent Guéguen
// Created on: vendredi 22 mai 2015, à 23h 20
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

#include "EvolutionSequenceSimulator.h"

#include "../NewLikelihood/OneProcessSequenceEvolution.h"
#include "../NewLikelihood/MixtureSequenceEvolution.h"
#include "../NewLikelihood/PartitionSequenceEvolution.h"
#include "../NewLikelihood/AutoCorrelationSequenceEvolution.h"
#include "../NewLikelihood/HmmSequenceEvolution.h"

#include <Bpp/Numeric/Random/RandomTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

EvolutionSequenceSimulator::EvolutionSequenceSimulator(const SequenceEvolution& evol) :
  SubstitutionProcessSequenceSimulator(evol),
  seqEvol_(&evol)
{
}

/******************************************************************************/

void EvolutionSequenceSimulator::resetSiteSimulators(size_t numberOfSites) const
{
  vMap_.resize(numberOfSites);

  const OneProcessSequenceEvolution* opse=dynamic_cast<const OneProcessSequenceEvolution*>(seqEvol_);

  if (opse)
  {
    size_t nProc=opse->getSubstitutionProcessNumbers()[0];
    for (size_t i=0; i<numberOfSites;i++)
      vMap_[i]=nProc;    
  }
  else
  {
    const PartitionSequenceEvolution* pse=dynamic_cast<const PartitionSequenceEvolution*>(seqEvol_);

    if (pse)
    {
      if (numberOfSites>pse->getNumberOfSites())
        throw IndexOutOfBoundsException("EvolutionSequenceSimulator::resetSiteSimulators",numberOfSites,0,pse->getNumberOfSites());

      for (size_t i=0; i<numberOfSites;i++)
        vMap_[i]=pse->getSubstitutionProcessNumber(i);
    }
    
    else
    {
      const vector<size_t> nProc=seqEvol_->getSubstitutionProcessNumbers();
      
      const MixtureSequenceEvolution* mse=dynamic_cast<const MixtureSequenceEvolution*>(seqEvol_);

      if (mse)
      {
        const vector<double>& vprob=mse->getSubProcessProbabilities();
        
        RandomTools::getSample(nProc, vprob, vMap_, true);
      }
      else
      {
        const AbstractHmmTransitionMatrix* htm=0;
        
        const AutoCorrelationSequenceEvolution* ase=dynamic_cast<const AutoCorrelationSequenceEvolution*>(seqEvol_);
        
        if (ase)
          htm=&ase->getHmmTransitionMatrix();

        const HmmSequenceEvolution* hse=dynamic_cast<const HmmSequenceEvolution*>(seqEvol_);
        
        if (hse)
          htm=&hse->getHmmTransitionMatrix();
        if (htm)
        {
          vector<size_t> vInd=htm->sample(numberOfSites);
          for (size_t i=0;i<numberOfSites;i++)
          {
            vMap_[i]=nProc[vInd[i]];
          }
          
        }
        else
          throw Exception("EvolutionSequenceSimulator::resetSiteSimulators : unknow Sequence Evolution.");
      }
    }
  }
}



