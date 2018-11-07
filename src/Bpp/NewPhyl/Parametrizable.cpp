//
// File: Parametrizable.cpp
// Authors:
// Created: 2017-06-06
// Last modified: 2017-06-06
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
#include <Bpp/NewPhyl/Parametrizable.h>
#include <Bpp/Numeric/Parametrizable.h>

using namespace std;

namespace bpp {
  namespace dataflow {

    std::unordered_map<std::string, std::shared_ptr<NumericMutable<double>>>
      createParameterMap (Context & c, const Parametrizable& parametrizable) {
      const auto & parameters = parametrizable.getParameters ();
      const auto nbParameters = parameters.size ();
      std::unordered_map<std::string, std::shared_ptr<NumericMutable<double>>> map;
      for (std::size_t i = 0; i < nbParameters; ++i) {
        const auto & param = parameters[i];
        map.emplace (param.getName (),
                     NumericMutable<double>::create (c, param.getValue ()));
      }
      return map;

    }

    NodeRefVec createDependencyVector (const Parametrizable & parametrizable,
                                       const std::function<NodeRef (const std::string &)> & getParameter) {
      const auto & parameters = parametrizable.getParameters ();
      const auto nbParameters = parameters.size ();
      NodeRefVec deps (nbParameters);
      for (std::size_t i = 0; i < nbParameters; ++i) {
        auto dep = getParameter (parameters[i].getName ());
        if (!dep) {
          throw Exception ("createDependencyVector (Parametrizable): parameter not found: " + parameters[i].getName ());
        }
        deps[i] = std::move (dep);
      }
      return deps;
    }

    
    // std::shared_ptr<ConfiguredParametrizable> ConfiguredParametrizable::create (Context & c, NodeRefVec && deps,
    //                                                                             std::unique_ptr<Parametrizable> && parametrizable) {
    //   if (!parametrizable) {
    //     throw Exception ("ConfiguredParametrizable(): nullptr parametrizable");
    //   }
    //   // Check dependencies
    //   const auto nbParameters = parametrizable->getParameters ().size ();
    //   checkDependenciesNotNull (typeid (Self), deps);
    //   checkDependencyVectorSize (typeid (Self), deps, nbParameters);
    //   checkDependencyRangeIsValue<double> (typeid (Self), deps, 0, nbParameters);
    //   return cachedAs<Self> (c, std::make_shared<Self> (std::move (deps), std::move (parametrizable)));
    // }

    // ConfiguredParametrizable::ConfiguredParametrizable (NodeRefVec && deps, std::unique_ptr<Parametrizable>&& parametrizable)
    //   : Value<const Parametrizable*> (std::move (deps), parametrizable.get ()), parametrizable_ (std::move(parametrizable)) {}

    // ConfiguredParametrizable::~ConfiguredParametrizable () = default;

    // // Model node additional arguments = (type of bpp::TransitionModel).
    // // Everything else is determined by the node dependencies.
    // bool ConfiguredParametrizable::compareAdditionalArguments (const Node & other) const {
    //   const auto * derived = dynamic_cast<const Self *> (&other);
    //   if (derived == nullptr) {
    //     return false;
    //   } else {
    //     const auto & thisParam = *parametrizable_;
    //     const auto & otherParam = *derived->parametrizable_;
    //     return typeid (thisParam) == typeid (otherParam);
    //   }
    // }
    
    // std::size_t ConfiguredParametrizable::hashAdditionalArguments () const {
    //   const auto & bppParam = *parametrizable_;
    //   return typeid (bppParam).hash_code ();
    // }

    // NodeRef ConfiguredParametrizable::recreate (Context & c, NodeRefVec && deps) {
    //   auto m = Self::create (c, std::move (deps), std::unique_ptr<Parametrizable>{parametrizable_->clone ()});
    //   m->config = this->config; // Duplicate derivation config
    //   return m;
    // }

    // void ConfiguredParametrizable::compute () {
    //   // Update each internal model bpp::Parameter with the dependency
    //   auto & parameters = parametrizable_->getParameters ();
    //   const auto nbParameters = this->nbDependencies ();
    //   for (std::size_t i = 0; i < nbParameters; ++i) {
    //     auto & v = accessValueConstCast<double> (*this->dependency (i));
    //     auto & p = parameters[i];
    //     if (p.getValue () != v) {
    //       // TODO improve bpp::Parametrizable interface to change values by index.
    //       parametrizable_->setParameterValue (parametrizable_->getParameterNameWithoutNamespace (p.getName ()), v);
    //     }
    //   }
    // }

  } // namespace dataflow
} // namespace bpp
