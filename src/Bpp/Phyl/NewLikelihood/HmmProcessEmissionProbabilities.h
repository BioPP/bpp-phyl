//
// File: HmmProcessEmissionProbabilities.h
// Created by: Laurent Guéguen
// Created on: mardi 24 septembre 2013, à 10h 00
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

#ifndef _HMM_PROCESS_EMISSION_PROBABILITIES_H_
#define _HMM_PROCESS_EMISSION_PROBABILITIES_H_

#include "HmmProcessAlphabet.h"
#include "PhyloLikelihoods/MultiProcessSequencePhyloLikelihood.h"

#include <Bpp/Numeric/Hmm/HmmEmissionProbabilities.h>
#include <Bpp/Numeric/AbstractParametrizable.h>

namespace bpp
{
  
/**
 * @brief Emission probabilities in the context of substitution
 * process.
 *
 */

  class HmmProcessEmissionProbabilities:
    public virtual HmmEmissionProbabilities,
    public AbstractParametrizable
  {
  private:
    const HmmProcessAlphabet* procAlph_;
    
    const MultiProcessSequencePhyloLikelihood* multiPL_;

    mutable std::vector<std::vector<double> > emProb_;

    mutable std::vector<std::vector<double> > dEmProb_;

    mutable std::vector<std::vector<double> > d2EmProb_;

    mutable bool upToDate_;

    void updateEmissionProbabilities_() const;

  public:
    HmmProcessEmissionProbabilities(const HmmProcessAlphabet* alphabet,
                                    const MultiProcessSequencePhyloLikelihood* multiPL);

    HmmProcessEmissionProbabilities(const HmmProcessEmissionProbabilities& hEP) :
      AbstractParametrizable(hEP),
      procAlph_(hEP.procAlph_),
      multiPL_(hEP.multiPL_),
      emProb_(hEP.emProb_),
      dEmProb_(hEP.dEmProb_),
      d2EmProb_(hEP.d2EmProb_),
      upToDate_(hEP.upToDate_)
    {}

    HmmProcessEmissionProbabilities& operator=(const HmmProcessEmissionProbabilities& hEP)
    {
      AbstractParametrizable::operator=(hEP);
      procAlph_=hEP.procAlph_;
      multiPL_=hEP.multiPL_;
      emProb_=hEP.emProb_;
      dEmProb_=hEP.dEmProb_;
      d2EmProb_=hEP.d2EmProb_;
      upToDate_=hEP.upToDate_;

      return *this;
    }

    HmmProcessEmissionProbabilities* clone() const { return new HmmProcessEmissionProbabilities(*this);}

    const HmmStateAlphabet* getHmmStateAlphabet() const
    {
      return procAlph_;
    }

    void fireParameterChanged(const ParameterList& parameters) {
      upToDate_=false;
    };

    /**
     * @brief Set the new hidden state alphabet.
     * @param stateAlphabet The new state alphabet.
     * @throw UnvalidStateAlphabetException if the new alphabet is uncorrect (for instance is NULL pointer).
     */
    
    void setHmmStateAlphabet(const HmmStateAlphabet* stateAlphabet);
    
    /**@} */

    /**
     * @brief Operator access to the emission probabilities.
     *
     * This is the fastest way to get the values, but no checking is performed on the indices.
     * For debugging purpose, the getPhyloEmissionProbability would be a safer use.
     *
     * @param pos The position of the sequential data to consider.
     * @param state The index of the hidden state to consider, as defined by the HmmStateAlphabet object associated to this class.
     */

    double operator()(size_t pos, size_t state) const
    {
      if (!upToDate_)
        updateEmissionProbabilities_();
      
      return emProb_[pos][state];
    }

    // void computeDEmissionProbabilities(std::string& variable) const;
  
    // void computeD2EmissionProbabilities(std::string& variable) const;
  
    // const std::vector<double>& getDEmissionProbabilities(size_t pos) const
    // {
    //   return dEmProb_[pos];
    // }
  
    // const std::vector<double>& getD2EmissionProbabilities(size_t pos) const
    // {
    //   return d2EmProb_[pos];
    // }


    /**
     * @brief Operator access to the emission probabilities.
     *
     * This is the fastest way to get the values, but no checking is performed on the indices.
     * For debugging purpose, the getPhyloEmissionProbability would be a safer use.
     *
     * @param pos The position of the sequential data to consider.
     * @return A vector of probabilities, whose size is the number of hidden states.
     */
    
    const std::vector<double>& operator()(size_t pos) const
    {
      if (!upToDate_)
        updateEmissionProbabilities_();

      return emProb_[pos];
    }
    
    /**
     * @return The number of positions in the data.
     */
    size_t getNumberOfPositions() const
    {
      return multiPL_->getNumberOfSites();
    }
    
  };

} //end of namespace bpp.

#endif //_HMM_PROCESS_EMISSION_PROBABILITIES_H_

