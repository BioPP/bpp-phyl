//
// File: Parameter.cpp
// Authors: Laurent Guéguen
// Created: mercredi 23 janvier 2019, à 06h 04
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

#include <Bpp/Exceptions.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/Parameter.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowCWise.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/Parametrizable.h>
#include <Bpp/Numeric/Parameter.h>

using namespace std;

namespace bpp {
  namespace dataflow {
    // Parameter node
    
    ConfiguredParameter::ConfiguredParameter (const Context& context, NodeRefVec&& deps, const Parameter& parameter)
      : Parameter(parameter), Value<Parameter*> (deps, this), context_(context)
    {
    };

    ConfiguredParameter::ConfiguredParameter (const Context& context, NodeRefVec&& deps, Parameter&& parameter)
      : Parameter(std::move(parameter)), Value<Parameter*> (deps, this), context_(context)
    {
    };

    ConfiguredParameter::~ConfiguredParameter () = default;

    std::string ConfiguredParameter::description () const {
      return "Parameter(" + getName () + ")\nvalue="+std::to_string(getValue());
    }
    
    std::string ConfiguredParameter::debugInfo () const {
      return "";
    }

    std::string ConfiguredParameter::color () const
    {
      auto& name=getName();
      if (name.substr(0,5)=="BrLen")
        return "#00ff00";
      else
        return "#ff8800";
    }

    // Parameter node additional arguments = (type of bpp::Parameter).
    // Everything else is determined by the node dependencies.

    bool ConfiguredParameter::compareAdditionalArguments (const Node & other) const {
      const auto * derived = dynamic_cast<const Self *> (&other);
      if (derived == nullptr) {
        return false;
      } else {
        return this->getName() == derived->getName();
      }
    }
    
    std::size_t ConfiguredParameter::hashAdditionalArguments () const {
      const auto & bppFS = *this;
      return typeid (bppFS).hash_code ();
    }

    NodeRef ConfiguredParameter::derive (Context & c, const Node & node) {
      if (&node == this) {
        return ConstantOne<double>::create (c, Dimension<double>());
      }
      else
        return this->dependency(0)->derive(c, node);
    }
    
    NodeRef ConfiguredParameter::recreate (Context & c) {
      return ConfiguredParameter::create (c, NodeRefVec{NumericMutable<double>::create(c, accessValueConstCast<double>(*dependency(0)))}, *this);
    }

    NodeRef ConfiguredParameter::recreate (Context & c, NodeRefVec && deps) {
      return ConfiguredParameter::create (c, std::move (deps), *this);
    }

    void ConfiguredParameter::compute () {
      // Update with the numerical dependency
      auto & v = accessValueConstCast<double> (*this->dependency (0));
      if (Parameter::getValue () != v) {
        Parameter::setValue(v);
      }
    }

    ////////////////////////////////////////////////////
    // ValueFromConfiguredParameter

    ValueFromConfiguredParameter::ValueFromConfiguredParameter (
      NodeRefVec && deps)
      : Value<double> (std::move (deps)) {}

    std::string ValueFromConfiguredParameter::debugInfo () const {
      return dependency(0)->debugInfo();
    }

    // ValueFromConfiguredParameter additional arguments = ().
    bool ValueFromConfiguredParameter::compareAdditionalArguments (const Node & other) const {
      const auto * derived = dynamic_cast<const Self *> (&other);
      return derived != nullptr;
    }

    std::shared_ptr<ValueFromConfiguredParameter> ValueFromConfiguredParameter::create (Context & c, NodeRefVec && deps) {
      checkDependenciesNotNull (typeid (Self), deps);
      checkDependencyVectorSize (typeid (Self), deps, 1);
      checkNthDependencyIs<ConfiguredParameter> (typeid (Self), deps, 0);
      return cachedAs<ValueFromConfiguredParameter> (c, std::make_shared<ValueFromConfiguredParameter> (std::move (deps)));
    }

    NodeRef ValueFromConfiguredParameter::derive (Context & c, const Node & node) {
      if (&node == this) {
        return ConstantOne<double>::create (c, Dimension<double>());
      }
      else
        return this->dependency(0)->derive(c, node);
    }
      
    std::string ValueFromConfiguredParameter::color () const
    {
      auto& name = accessValueConstCast<const ConfiguredParameter *> (*this->dependency (0))->getName();
      if (name.substr(0,5)=="BrLen")
        return "#66ff66";
      else
        return "#ffcc00";
    }

    NodeRef ValueFromConfiguredParameter::recreate (Context & c, NodeRefVec && deps)
    {
      return ValueFromConfiguredParameter::create (c, {std::move (deps)});
    }

    void ValueFromConfiguredParameter::compute () {
      const auto * param = accessValueConstCast<const ConfiguredParameter *> (*this->dependency (0));
      this->accessValueMutable () = param->getValue();
    }

  } // namespace dataflow
} // namespace bpp
