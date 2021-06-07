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
#include "../DataFlow/DataFlowNumeric.h"
#include "../DataFlow/Parametrizable.h"
#include "../DataFlow/LikelihoodCalculation.h"

namespace bpp
{
  class LikelihoodCalculation;
  
  class AbstractPhyloLikelihood :
    public virtual PhyloLikelihood
  {
  public:
    // Cache generated nodes representing derivatives, to avoid recreating them every time.
    // Using the mutable keyword because the table must be changed even in const methods.
    struct StringPairHash {
      std::size_t operator() (const std::pair<std::string, std::string> & p) const {
        std::hash<std::string> strHash{};
        return strHash (p.first) ^ (strHash (p.second) << 1);
      }
    };
      
  protected:
    Context & context_;
      
    /**
     * @brief the value
     *
     **/
      
    mutable double minusLogLik_;

    /**
     * @brief For Dataflow computing
     *
     */
      
    mutable std::unordered_map<std::string, ValueRef<double>> firstOrderDerivativeNodes_;

    mutable std::unordered_map<std::pair<std::string, std::string>, ValueRef<double>,
                               StringPairHash>
    secondOrderDerivativeNodes_;


  public:
    AbstractPhyloLikelihood(Context& context) :
      context_(context),
      minusLogLik_(0)
    {
    }

    AbstractPhyloLikelihood(const AbstractPhyloLikelihood& asd) :
      context_(asd.context_),
      minusLogLik_(asd.minusLogLik_)
    {
      shareParameters(asd.getParameters());
    }
      
    virtual ~AbstractPhyloLikelihood() {}

    AbstractPhyloLikelihood* clone() const override = 0;

    Context& getContext()
    {
      return context_;
    }

    /**
     * @brief Sets the computeLikelihoods_ to true.
     *
     */
      
    virtual bool isInitialized() const  override
    {
      return false;
    }

  public:

    /*
     *
     *@ brief Share Parameters, that are DF_parameters
     *
     */
     
    void shareParameters(const ParameterList& variableNodes)
    {
      this->getParameters_().shareParameters(variableNodes);
    }

    void setParameters(const ParameterList& parameters) override
    {
      setParametersValues(parameters);
    }

    /*
     *@ Return the LikelihoodCalculation.
     *
     */
    
    virtual std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation () const = 0;

    /*
     *@ Return the LikDF node where the Likelihood is computed.
     *
     */
    
    ValueRef<double> getLikelihoodNode() const {
      return getLikelihoodCalculation()->getLikelihoodNode();
    }


    // bpp::Function

    /**
     * @brief Tell if derivatives must be computed: for Function
     * inheritance.
     *
     */
    
    virtual void enableFirstOrderDerivatives(bool yn)  override {};
    virtual void enableSecondOrderDerivatives(bool yn)  override {};
    bool enableFirstOrderDerivatives() const  override{ return true; }
    bool enableSecondOrderDerivatives() const  override{ return true; }

    /*
     * @brief return the value, ie -loglikelihood
     *
     * !!! check on computeLikelihoods_ is not done here.
     *
     */
      
    double getValue() const  override
    {
      if (!isInitialized())
        throw Exception("AbstractPhyloLikelihood::getValue(). Instance is not initialized.");

      if (!getLikelihoodNode())
        throw Exception("AbstractPhyloLikelihood::getValue(). LikelihoodNode is not built.");

      
      minusLogLik_=-getLikelihoodNode()->getTargetValue();
      return minusLogLik_;
    }

    // bpp::DerivableFirstOrder
    double getFirstOrderDerivative (const std::string & variable) const override {
      return -firstOrderDerivativeNode (variable)->getTargetValue ();
    }

    // Get nodes of derivatives directly
    ValueRef<double> firstOrderDerivativeNode (const std::string & variable) const {
      const auto it = firstOrderDerivativeNodes_.find (variable);
      if (it != firstOrderDerivativeNodes_.end ()) {
        return it->second;
      } else {
        auto node = getLikelihoodNode()->deriveAsValue (context_, accessVariableNode (variable));
        firstOrderDerivativeNodes_.emplace (variable, node);
        return node;
      }
    }
      
    // bpp::DerivableSecondOrder
    double getSecondOrderDerivative (const std::string & variable) const override {
      return getSecondOrderDerivative (variable, variable);
    }

    double getSecondOrderDerivative (const std::string & variable1,
                                     const std::string & variable2) const override {
      return -secondOrderDerivativeNode (variable1, variable2)->getTargetValue ();
    }

    ValueRef<double> secondOrderDerivativeNode (const std::string & variable1,
                                                const std::string & variable2) const {
      const auto key = std::make_pair (variable1, variable2);
      const auto it = secondOrderDerivativeNodes_.find (key);
      if (it != secondOrderDerivativeNodes_.end ()) {
        return it->second;
      } else {
        // Reuse firstOrderDerivative() to generate the first derivative with caching
        auto node =
          firstOrderDerivativeNode (variable1)->deriveAsValue (context_, accessVariableNode (variable2));
        secondOrderDerivativeNodes_.emplace (key, node);
        return node;
      }
    }

  protected:
    static Node_DF & accessVariableNode (const Parameter & param) {
      return *dynamic_cast<const ConfiguredParameter&>(param).dependency(0);
    }
      
    Node_DF & accessVariableNode (const std::string & name) const {
      return accessVariableNode (getParameter (name));
    }

  };

} //end of namespace bpp.

#endif  //_ABSTRACTPHYLOLIKELIHOOD_H_

