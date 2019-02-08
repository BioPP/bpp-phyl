//
// File: Parameter.h
// Authors: Laurent Gueguen 
// Created: mercredi 23 janvier 2019, à 06h 02
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

#pragma once
#ifndef BPP_NEWPHYL_PARAMETER_H
#define BPP_NEWPHYL_PARAMETER_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Exceptions.h>
#include <functional>
#include <unordered_map>

namespace bpp {

  namespace dataflow {

    /** @brief Data flow node representing a Frequencies Set
     * configured with parameter values.
     *
     * This class wraps a bpp::Parameter as a data flow node.
     *
     * It depends on Value<double> nodes (one for each parameter
     * declared in the freq set).
     * It provides a dummy value representing the "frequencies set
     * configured by its parameters".
     *
     * The dummy value is implemented as a pointer to the internal
     * frequencies set for simplicity.
     */
    
    class ConfiguredParameter : public Parameter, public Value<Parameter*>
    {
    private:

      const Context& context_;
      
    public:
      using Self = ConfiguredParameter;
      using Target = Parameter;
      
      /// Build a new ConfiguredParameter node.
      static std::shared_ptr<Self> create (Context & c, NodeRefVec&& deps, const Parameter& param)
      {
        checkDependenciesNotNull (typeid (Self), deps);
        checkDependencyVectorSize (typeid (Self), deps, 1);
        checkNthDependencyIsValue<double> (typeid (Self), deps, 0);

        return cachedAs<Self> (c, std::make_shared<Self> (c, std::move(deps), param));
      }

      ConfiguredParameter (const Context& context, NodeRefVec&& deps, const Parameter& parameter);

      ConfiguredParameter* clone() const {
        throw Exception("ConfiguredParameter::clone should not be called.");
      }
      
      ~ConfiguredParameter ();
      
      std::string description () const;
      std::string debugInfo () const;

      void setValue(double v)
      {
        Parameter::setValue(v);
        static_cast<NumericMutable<double>&>(*dependency(0)).setValue(Parameter::getValue());
      }
    
      bool compareAdditionalArguments (const Node & other) const;
      
      std::size_t hashAdditionalArguments () const;
      
      NodeRef derive (Context & c, const Node & node);

      NodeRef recreate (Context & c, NodeRefVec && deps);
      
      NodeRef recreate (Context & c);

    private:
      void compute ();
    };

  } // namespace dataflow
} // namespace bpp

#endif // BPP_NEWPHYL_PARAMETER_H
