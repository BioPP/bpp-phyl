//
// File: EvolutionSequenceSimulator.h
// Created by: Laurent Guéguen
// Created on: vendredi 22 mai 2015, à 23h 03
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

#ifndef _EVOLUTION_SEQUENCESIMULATOR_H_
#define _EVOLUTION_SEQUENCESIMULATOR_H_


#include "SubstitutionProcessSequenceSimulator.h"


#include "../NewLikelihood/SequenceEvolution.h"

namespace bpp
{
  
  class EvolutionSequenceSimulator:
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
  

} //end of namespace bpp.

#endif //_EVOLUTION_SEQUENCESIMULATOR_H_

