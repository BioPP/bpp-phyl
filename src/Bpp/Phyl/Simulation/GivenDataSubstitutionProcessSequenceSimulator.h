//
// File: GivenDataSubstitutionProcessSequenceSimulator.h
// Created by: Laurent Guéguen
// Created on: mercredi 21 octobre 2020, à 09h 45
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

#ifndef _GIVEN_DATA_SUBSTITUTION_PROCESS_SEQUENCE_SIMULATOR_H_
#define _GIVEN_DATA_SUBSTITUTION_PROCESS_SEQUENCE_SIMULATOR_H_

#include "GivenDataSubstitutionProcessSiteSimulator.h"

#include "SequenceSimulator.h"

namespace bpp
{

/**
 * @brief Sequences simulation under a unique substitution process, but with site specific
 * posterior probabilities.
 *
 */

  class GivenDataSubstitutionProcessSequenceSimulator:
    public virtual SequenceSimulator
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
    
  public:    
    GivenDataSubstitutionProcessSequenceSimulator(std::shared_ptr<LikelihoodCalculationSingleProcess> calcul) :
      calcul_(calcul), vSiteSim_()
    {
      for (size_t i = 0; i<calcul_->getNumberOfDistinctSites(); i++)
        vSiteSim_.push_back(std::make_shared<GivenDataSubstitutionProcessSiteSimulator>(calcul_, i, true));
    }

    virtual ~GivenDataSubstitutionProcessSequenceSimulator()
    {
    }

    GivenDataSubstitutionProcessSequenceSimulator(const GivenDataSubstitutionProcessSequenceSimulator& nhss) :
      calcul_(nhss.calcul_), vSiteSim_(nhss.vSiteSim_)
    {};
    
    GivenDataSubstitutionProcessSequenceSimulator* clone() const { return new GivenDataSubstitutionProcessSequenceSimulator(*this); }

  public:
  
    /**
     * @name The SequenceSimulator interface
     * Here the numberOfSites is unused (awkward inheritance...)
     *
     */

    std::shared_ptr<SiteContainer> simulate() const
    {
      return simulate(0);
    }

    std::shared_ptr<SiteContainer> simulate(size_t numberOfSites) const;

    const SiteSimulator& getSiteSimulator(size_t pos) const
    {
      return *vSiteSim_[calcul_->getRootArrayPosition(pos)];
    }    
    
    /**
     * @name SiteSimulator and SequenceSimulator interface
     *
     * @{
     */
    const Alphabet* getAlphabet() const { return vSiteSim_[0]->getAlphabet(); }

    std::vector<std::string> getSequencesNames() const {
      return vSiteSim_[0]->getSequencesNames();
    }


    /** @} */

    /**
     * @brief Sets whether we will output the internal sequences or not.
     *
     *
     * @param yn Tell if we should output internal sequences.
     */
    void outputInternalSequences(bool yn)
    {
      for (auto& siteSim : vSiteSim_)
        siteSim->outputInternalSites(yn);
    }

  };


} //end of namespace bpp.

#endif //_GIVEN_DATA_SUBSTITUTION_PROCESS_SEQUENCE_SIMULATOR_H_

