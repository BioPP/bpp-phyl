//
// File: LikelihoodCalculation.h
// Authors: François Gindraud, Laurent Guéguen (2018)
// Created: jeudi 28 février 2019, à 07h 22
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef LIKELIHOOD_CALCULATION_H
#define LIKELIHOOD_CALCULATION_H

#include "Bpp/Phyl/NewLikelihood/DataFlow/DataFlow.h"
#include <Bpp/Numeric/AbstractParametrizable.h>

namespace bpp {

  /** Base class for Likelihood calcucations
   *  
   */


  class LikelihoodCalculation :
    public AbstractParametrizable
  {
  private:
    
    Context& context_;

  protected:
    /******************************************/
    /** Likelihoods  **/
          
    ValueRef<double> likelihood_;

  public:
    LikelihoodCalculation(Context & context) :
      AbstractParametrizable(""),
      context_(context)
    {}

    LikelihoodCalculation(Context & context,
                          ParameterList& paramList) : 
      AbstractParametrizable(""),
      context_(context)
    {  
      shareParameters_(paramList);
    }


    LikelihoodCalculation(const LikelihoodCalculation& lik) :
      AbstractParametrizable(lik),
      context_(*std::shared_ptr<Context>().get())
    {};
    
    virtual ValueRef<double> getLikelihoodNode()
    {
      makeLikelihoods();
      return likelihood_;
    }

    double getLogLikelihood() 
    {
      return getLikelihoodNode()->getTargetValue();
    }

    virtual bool isInitialized() const = 0;

    const Context& getContext() const {
      return context_;
    }

    virtual void makeLikelihoods() = 0;
    
    /*
     * @brief Return likelihood_ without any computation, used to build likelihood_.
     *
     */
    
    void setLikelihoodNode(ValueRef<double> ll)
    {
      likelihood_=ll;
    }

  protected:

    Context& getContext_() {
      return context_;
    }

    // ValueRef<double> getLikelihoodNode()
    // {
    //   return likelihood_;
    // }
  };

} // namespace bpp

#endif // LIKELIHOOD_CALCULATION_H

