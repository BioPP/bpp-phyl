//
// File: DataFlowWrappers.h
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

#ifndef DATA_FLOW_WRAPPERS_H
#define DATA_FLOW_WRAPPERS_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>

#include "DataFlowNumeric.h"
#include "Parameter.h"

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
     * FIXME This temporary interface to bpp::DerivableSecondOrder stuff
     * should be improved.
     *
     * Any bpp::Parameter can be given in the bpp::ParameterList, but
     * only DataFlowParameter are supported because we need to have
     * dataflow nodes.
     * No specific check is done, in case of error it will be a
     * std::bad_cast from dynamic_cast.
     *
     * In addition, as we need a context for derivation but the bpp
     * legacy API does not support it, a reference is stored in the
     * class which is dangerous with respect to lifetime.
     */
  
    class DataFlowFunction : public DerivableSecondOrder {
    private:
      dataflow::Context & context_;

      // Store nodes
      dataflow::ValueRef<double> resultNode_;
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
      DataFlowFunction (dataflow::Context & context, dataflow::ValueRef<double> resultNode,
                        const ParameterList & variableNodes)
        : context_ (context), resultNode_ (std::move (resultNode)), variableNodes_ () {
        variableNodes_.shareParameters(variableNodes);
      }

      // Legacy boilerplate
      DataFlowFunction * clone () const override { return new DataFlowFunction (*this); }

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
      void setParameters (const ParameterList & params) override {
        variableNodes_.setParametersValues (params);
      }
      
      double getValue () const override { return resultNode_->getTargetValue (); }

      ValueRef<double> getValueRef () const { return resultNode_; }

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
          auto node = resultNode_->deriveAsValue (context_, accessVariableNode (variable));
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

#endif // DATA_FLOW_WRAPPERS_H
