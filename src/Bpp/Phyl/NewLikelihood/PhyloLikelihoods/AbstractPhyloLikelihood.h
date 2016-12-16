//
// File: AbstractPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: lundi 25 avril 2016, à 23h 32
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

#ifndef _ABSTRACTPHYLOLIKELIHOOD_H_
#define _ABSTRACTPHYLOLIKELIHOOD_H_

#include "PhyloLikelihood.h"

namespace bpp
{
  class AbstractPhyloLikelihood :
    public virtual PhyloLikelihood
  {
  protected:
      
    /**
     * @brief the value
     *
     **/
      
    mutable double minusLogLik_;

    /**
     * @brief sey if derivatives should be computed
     *
     */
      
    bool computeFirstOrderDerivatives_;
    bool computeSecondOrderDerivatives_;

    // say if the Likelihoods should be recomputed
      
    mutable bool computeLikelihoods_;

    // maps of the updated values for derivatives of the logLik

    mutable std::map<std::string, double> dValues_;
    mutable std::map<std::string, double> d2Values_;
    
    // say if initialized
      
    mutable bool initialized_;

  public:
    AbstractPhyloLikelihood() :
      minusLogLik_(0),
      computeFirstOrderDerivatives_(true),
      computeSecondOrderDerivatives_(true),
      computeLikelihoods_(true),
      dValues_(),
      d2Values_(),
      initialized_(false)
    {
    }

    AbstractPhyloLikelihood(const AbstractPhyloLikelihood& asd) :
      minusLogLik_(asd.minusLogLik_),
      computeFirstOrderDerivatives_(asd.computeFirstOrderDerivatives_),
      computeSecondOrderDerivatives_(asd.computeSecondOrderDerivatives_),
      computeLikelihoods_(asd.computeLikelihoods_),
      dValues_(asd.dValues_),
      d2Values_(asd.d2Values_),
      initialized_(asd.initialized_)
    {
    }
      
    AbstractPhyloLikelihood& operator=(const AbstractPhyloLikelihood& asd)
    {
      minusLogLik_                   = asd.minusLogLik_;
      computeFirstOrderDerivatives_  = asd.computeFirstOrderDerivatives_;
      computeSecondOrderDerivatives_ = asd.computeSecondOrderDerivatives_;

      computeLikelihoods_ = asd.computeLikelihoods_;
      dValues_ = asd.dValues_;
      d2Values_ = asd.d2Values_;

      initialized_ = asd.initialized_;
        
      return *this;
    }

    virtual ~AbstractPhyloLikelihood() {}

    AbstractPhyloLikelihood* clone() const = 0;

  public:
    /**
     * @brief Sets the computeLikelihoods_ to true.
     *
     */
      
    void update()
    {
      computeLikelihoods_ = true;
      dValues_.clear();
      d2Values_.clear();
    }

    virtual void initialize()
    {
      initialized_=true;
    }

  public:

    void setParameters(const ParameterList& parameters) throw (ParameterNotFoundException, ConstraintException)
    {
      setParametersValues(parameters);
    }

    virtual void enableDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
    virtual void enableFirstOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = yn; }
    virtual void enableSecondOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
    bool enableFirstOrderDerivatives() const { return computeFirstOrderDerivatives_; }
    bool enableSecondOrderDerivatives() const { return computeSecondOrderDerivatives_; }

    bool isInitialized() const { return initialized_; }

    /*
     * @brief return the value, ie -loglikelihood
     *
     * !!! check on computeLikelihoods_ is not done here.
     *
     */
      
    double getValue() const throw (Exception)
    {
      if (!isInitialized())
        throw Exception("AbstractPhyloLikelihood::getValue(). Instance is not initialized.");

      minusLogLik_=-getLogLikelihood();

      return minusLogLik_;
    }


  };

} //end of namespace bpp.

#endif  //_ABSTRACTPHYLOLIKELIHOOD_H_

