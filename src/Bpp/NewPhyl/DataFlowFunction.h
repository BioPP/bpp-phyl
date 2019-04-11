//
// File: DataFlowFunction.h
// Authors: François Gindraud (2017)
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

#ifndef DATA_FLOW_FUNCTION_H
#define DATA_FLOW_FUNCTION_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>

#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/PhyloLikelihood.h>

#include "DataFlowNumeric.h"
#include "Parameter.h"
#include "LikelihoodCalculation.h"

#include <unordered_map>

/* This file contains wrappers.
 * They are used to bridge the gap between bpp::dataflow stuff and the rest of bpp.
 *
 */

namespace bpp {

  namespace dataflow
  {
    
    /* Wraps a dataflow graph as a function: resultNode = f(variableNodes).
     *
     */
  
    class DataFlowFunction : public PhyloLikelihood {
    private:
      dataflow::Context & context_;

      // Store nodes
      std::shared_ptr<dataflow::LikelihoodCalculation> likCal_;
      // ValueRef<double> resultNode_;
      
      ParameterList variableNodes_;

      // Cache generated nodes representing derivatives, to avoid recreating them every time.
      // Using the mutable keyword because the table must be changed even in const methods.
      struct StringPairHash {
        std::size_t operator() (const std::pair<std::string, std::string> & p) const {
          std::hash<std::string> strHash{};
          return strHash (p.first) ^ (strHash (p.second) << 1);
        }
      };

      
      mutable std::unordered_map<std::string, dataflow::ValueRef<double>> firstOrderDerivativeNodes_;
      mutable std::unordered_map<std::pair<std::string, std::string>, dataflow::ValueRef<double>,
                                 StringPairHash>
      secondOrderDerivativeNodes_;

    public:
      DataFlowFunction (dataflow::Context & context, std::shared_ptr<LikelihoodCalculation> likCal,
                        const ParameterList & variableNodes)
        : context_ (context), likCal_(likCal), variableNodes_ () {
        variableNodes_.shareParameters(variableNodes);
      }

      /*
       * @brief: the parameters are those of the LikelihoodCalculation
       */
      
      DataFlowFunction (dataflow::Context & context, std::shared_ptr<LikelihoodCalculation> likCal)
        : context_ (context), likCal_(likCal), variableNodes_ () {
        variableNodes_.shareParameters(likCal->getParameters());
      }

      // Legacy boilerplate
      DataFlowFunction * clone () const override { return new DataFlowFunction (*this); }

      /**
       * @return initialize the likelihood function.
       */
      
      void initialize() {};
      
      /**
       * @return 'true' is the likelihood function has been initialized.
       */
      bool isInitialized() const {
        return likCal_->getData();
      };
      
      /**
       * @}
       */

      /**
       * @brief set it arrays should be computed in log.
       *
       */

      void setUseLog(bool useLog) {};

      /**
       * @name The likelihood functions.
       *
       * @{
       */
      
      /**
       * @brief update the likelihood to get ready for computation
       *
       */

      void updateLikelihood() const {};

      /**
       * @brief compute the likelihood
       *
       */

      void computeLikelihood() const {};

      /**
       * @brief Get the logarithm of the likelihood for the whole dataset.
       *
       * @return The logarithm of the likelihood of the dataset.
       */
    
      double getLogLikelihood() const {
        return getValue();
      }
      
      /**
       * @brief Compute the derivates of the LogLikelihood.
       *
       */

      void computeDLogLikelihood_(const std::string& variable) const {};

    
      void computeD2LogLikelihood_(const std::string& variable) const {};
    
      /**
       * @brief Get the derivates of the LogLikelihood.
       *
       */

      double getDLogLikelihood(const std::string& variable) const
      {
        return getFirstOrderDerivative(variable);
      }

      double getD2LogLikelihood(const std::string& variable) const
      {
        return getSecondOrderDerivative(variable);
      }

      /** @} */

      /**
       * @name Retrieve some particular independent parameters subsets.
       *
       * @{
       */
    
      /**
       * @brief Get the independent branch lengths parameters.
       *
       * @return A ParameterList with all branch lengths.
       */

      ParameterList getBranchLengthParameters() const
      {
        return likCal_->getSubstitutionProcess().getBranchLengthParameters(true);
      }
      
      /**
       * @brief Get the independent parameters associated to substitution model(s).
       *
       * @return A ParameterList.
       */

      ParameterList getSubstitutionModelParameters() const
      {
        return likCal_->getSubstitutionProcess().getSubstitutionModelParameters(true);
      }

      /**
       * @brief Get the independent parameters associated to the rate distribution(s).
       *
       * @return A ParameterList.
       */

      ParameterList getRateDistributionParameters() const
      {
        return likCal_->getSubstitutionProcess().getRateDistributionParameters(true);
      }
      
      /**
       * @brief Get the independent parameters associated to the root
       * frequencies(s).
       *
       * @return A ParameterList.
       */
      
      ParameterList getRootFrequenciesParameters() const
      {
        return likCal_->getSubstitutionProcess().getRootFrequenciesParameters(true);
      }
      
      /**
       * @brief All independent non derivable parameters.
       *
       * Usually, this contains all substitution model parameters and rate distribution.
       *
       * @return A ParameterList.
       */

      ParameterList getNonDerivableParameters() const
      {
        return likCal_->getSubstitutionProcess().getNonDerivableParameters();
      }
     

      /** @} */

      /**
       * @brief Tell if derivatives must be computed.
       *
       * This methods calls the enableFirstOrderDerivatives and enableSecondOrderDerivatives.
       *
       * @param yn Yes or no.
       */
      void enableDerivatives(bool yn) {};
      
      
      // bpp::Parametrizable (prefix unused FIXME?)
      bool hasParameter (const std::string & name) const override { return variableNodes_.hasParameter (name); }
      const ParameterList & getParameters () const override { return variableNodes_; }
      const Parameter & getParameter (const std::string & name) const override {
        return variableNodes_.getParameter (name);
      }
      double getParameterValue (const std::string & name) const override {
        return variableNodes_.getParameterValue (name);
      }
      void setAllParametersValues (const ParameterList & params) override {
        return variableNodes_.setAllParametersValues (params);
      }
      void setParameterValue (const std::string & name, double value) override {
        variableNodes_.setParameterValue (name, value);
      }
      void setParametersValues (const ParameterList & params) override {
        variableNodes_.setParametersValues (params);
      }
      bool matchParametersValues (const ParameterList & params) override {
        return variableNodes_.matchParametersValues (params);
      }
      std::size_t getNumberOfParameters () const override { return variableNodes_.size (); }

      void setNamespace (const std::string &) override {}
      std::string getNamespace () const override { return {}; }
      std::string getParameterNameWithoutNamespace (const std::string & name) const override { return name; }

      // bpp::Function
      void setParameters (const ParameterList & params) override
      {
        variableNodes_.setParametersValues (params);
      }
      
      double getValue () const override
      {
        return likCal_->getLikelihood()->getTargetValue ();
      }

      std::shared_ptr<dataflow::LikelihoodCalculation> getLikelihoodCalculation()
      {
        return likCal_;
      }
      
      // bpp::DerivableFirstOrder
      void enableFirstOrderDerivatives (bool) override {}
      bool enableFirstOrderDerivatives () const override { return true; }
      double getFirstOrderDerivative (const std::string & variable) const override {
        return firstOrderDerivativeNode (variable)->getTargetValue ();
      }

      // bpp::DerivableSecondOrder
      void enableSecondOrderDerivatives (bool) override {}
      bool enableSecondOrderDerivatives () const override { return true; }
      double getSecondOrderDerivative (const std::string & variable) const override {
        return getSecondOrderDerivative (variable, variable);
      }

      double getSecondOrderDerivative (const std::string & variable1,
                                       const std::string & variable2) const override {
        return secondOrderDerivativeNode (variable1, variable2)->getTargetValue ();
      }

      // Get nodes of derivatives directly
      dataflow::ValueRef<double> firstOrderDerivativeNode (const std::string & variable) const {
        const auto it = firstOrderDerivativeNodes_.find (variable);
        if (it != firstOrderDerivativeNodes_.end ()) {
          return it->second;
        } else {
          auto node = likCal_->getLikelihood()->deriveAsValue (context_, accessVariableNode (variable));
          firstOrderDerivativeNodes_.emplace (variable, node);
          return node;
        }
      }
      
      dataflow::ValueRef<double> secondOrderDerivativeNode (const std::string & variable1,
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

    private:
      static Node & accessVariableNode (const Parameter & param) {
        return *dynamic_cast<const ConfiguredParameter&>(param).dependency(0);
      }
      
      Node & accessVariableNode (const std::string & name) const {
        return accessVariableNode (getParameter (name));
      }
    };

  }
  
} // namespace bpp

#endif // DATA_FLOW_FUNCTION_H
