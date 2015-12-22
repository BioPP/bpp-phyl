//
// File: HmmPhyloAlphabet.h
// Created by: Laurent Guéguen
// Created on: mardi 27 octobre 2015, à 19h 09
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

#ifndef _HMM_PHYLO_ALPHABET_H_
#define _HMM_PHYLO_ALPHABET_H_

//From Numeric
#include <Bpp/Numeric/Hmm/HmmStateAlphabet.h>

#include "SetOfAlignedPhyloLikelihood.h"

namespace bpp {

  /**
   * @brief Hidden states alphabet.
   *
   * Implementation of HmmStateAlphabet where Alphabet States are all
   * the PhyloLikelihoods belonging to a SetOfAlignedPhyloLikelihood.
   *
   */
  
  class HmmPhyloAlphabet:
    public virtual HmmStateAlphabet,
    public AbstractParametrizable
  {
  private:
    /*
     * @brief vector of aligned phylolikelihoods. These
     * phylolikelihoods are not owned by the alphabet.
     *
     */
    
    std::vector<const AlignedPhyloLikelihood*> vAP_;

    size_t nbSites_;
    
  public:
    HmmPhyloAlphabet(const SetOfAlignedPhyloLikelihood& soap):
      AbstractParametrizable(""),
      vAP_(),
      nbSites_(0)
    {
      const std::vector<size_t>& nphyl=soap.getNumbersOfPhyloLikelihoods();
      
      for (size_t i=0; i<nphyl.size(); i++)
      {
        const AlignedPhyloLikelihood* ap=soap.getAbstractPhyloLikelihood(nphyl[i]);
        vAP_.push_back(ap);
        includeParameters_(ap->getParameters());
      }
      
      nbSites_=soap.getNumberOfSites();
    }
    
    
    HmmPhyloAlphabet(const HmmPhyloAlphabet& hpa) :
      AbstractParametrizable(hpa),
      vAP_(hpa.vAP_),
      nbSites_(hpa.nbSites_)
    {}
    
    HmmPhyloAlphabet& operator=(const HmmPhyloAlphabet& hpa)
    {
      AbstractParametrizable::operator=(hpa);

      vAP_=hpa.vAP_;
      nbSites_=hpa.nbSites_;

      return *this;
      
    }

    HmmPhyloAlphabet* clone() const {return new HmmPhyloAlphabet(*this);}

    ~HmmPhyloAlphabet() {};

    /**
     * @param stateIndex The index of a hidden state.
     * @return The corresponding hidden state.
     * @see getNumberOfStates
     */

    const Clonable& getState(size_t stateIndex) const throw (HmmBadStateException)
    {
      return *vAP_[stateIndex];
    }

    const AlignedPhyloLikelihood& getPhyloLikelihood(size_t stateIndex) const
    {
      return *vAP_[stateIndex];
    }

    size_t getNumberOfStates() const
    {
      return vAP_.size();
    }

    size_t getNumberOfSites() const
    {
      return nbSites_;
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

    void fireParameterChanged(const ParameterList& parameters)
    {
    }
    

    /**
     * @{
     *
     * Functions to perform computations
     *
     */

    void updateLikelihood() const
    {
      for (size_t i=0; i<vAP_.size(); i++)
        vAP_[i]->updateLikelihood();
    }

    void computeLikelihood() const
    {
      for (size_t i=0; i<vAP_.size(); i++)
        vAP_[i]->computeLikelihood();
    }

    void computeDLogLikelihood(const std::string& variable) const
    {
      for (size_t i=0; i<vAP_.size(); i++)
        vAP_[i]->computeDLogLikelihood_(variable);
    }

    void computeD2LogLikelihood(const std::string& variable) const
    {
      for (size_t i=0; i<vAP_.size(); i++)
        vAP_[i]->computeD2LogLikelihood_(variable);
    }

    /**
     * @}
     *
     */

  };

}

#endif //_HMMPROCESSALPHABET_H_

