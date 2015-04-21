//
// File: HmmProcessAlphabet.h
// Created by: Laurent Guéguen
// Created on: vendredi 20 septembre 2013, à 23h 46
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _HMMPROCESSALPHABET_H_
#define _HMMPROCESSALPHABET_H_

//From Numeric
#include <Bpp/Numeric/Hmm/HmmStateAlphabet.h>

#include "SubstitutionProcessCollection.h"


namespace bpp {

  /**
   * @brief Hidden states alphabet.
   *
   * Implementation of HmmStateAlphabet where Alphabet States are
   * Substitution Process belonging to a collection.
   *
   */
  
  class HmmProcessAlphabet:
    public virtual HmmStateAlphabet,
    public AbstractParametrizable
  {
  private:
    const SubstitutionProcessCollection* processColl_;

    /*
     * @brief the vector of the substitution process numbers.
     *
     */
    
    std::vector<size_t> nProc_;
    
  public:
    HmmProcessAlphabet(const SubstitutionProcessCollection* pSub, std::vector<size_t> nProc) :
      AbstractParametrizable(""),
      processColl_(pSub),
      nProc_(nProc)
    {
    }

    HmmProcessAlphabet(const HmmProcessAlphabet& hpa) :
      AbstractParametrizable(hpa),
      processColl_(hpa.processColl_),
      nProc_(hpa.nProc_)
    {};

    HmmProcessAlphabet& operator=(const HmmProcessAlphabet& hpa){
      AbstractParametrizable::operator=(*this);
      processColl_=hpa.processColl_;
      nProc_=hpa.nProc_;
        
      return *this;
    }

    HmmProcessAlphabet* clone() const {return new HmmProcessAlphabet(*this);}

    ~HmmProcessAlphabet() {};

    void fireParameterChanged(const ParameterList& pl) {};
    
    /**
     * @param stateIndex The index of a hidden state.
     * @return The corresponding hidden state.
     * @see getNumberOfStates
     */

    const Clonable& getState(size_t stateIndex) const throw (HmmBadStateException)
    {
      return processColl_->getSubstitutionProcess(nProc_[stateIndex]);
    }
      
    size_t getNumberOfStates() const
    {
      return nProc_.size();
    }

    /**
     * @brief Tell if this instance can work with the instance of alphabet given as input.
     *
     * In many case, this will return true if the pointer provided as argument refers to this object.
     *
     * @param stateAlphabet The alphabet to check.
     * @return true if the matrix is compatible with the given alphabet.
     */
    
    bool worksWith(const HmmStateAlphabet* stateAlphabet) const
    {
      return stateAlphabet==this;
    }
 
  };

}

#endif //_HMMPROCESSALPHABET_H_

