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
#include <Bpp/Phyl/NewLikelihood/DataFlow/Parametrizable.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/Parameter.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/ParameterAliasable.h>

using namespace std;

namespace bpp {
    std::unordered_map<std::string, std::shared_ptr<ConfiguredParameter>>
      createParameterMap (Context & c, const ParameterAliasable& parametrizable) {
      const auto & parameters = parametrizable.getIndependentParameters ();
      const auto nbParameters = parameters.size ();
      std::unordered_map<std::string, std::shared_ptr<ConfiguredParameter>> map;
      for (std::size_t i = 0; i < nbParameters; ++i) {
        const auto & param = parameters[i];
        auto value = NumericMutable<double>::create (c, param.getValue ());
        map.emplace (param.getName (),
                     ConfiguredParameter::create (c, {std::move(value)}, param));
      }
      return map;
    }

    std::unordered_map<std::string, std::shared_ptr<ConfiguredParameter>>
      createParameterMap (Context & c, const Parametrizable& parametrizable) {
      const auto & parameters = parametrizable.getParameters ();
      const auto nbParameters = parameters.size ();
      std::unordered_map<std::string, std::shared_ptr<ConfiguredParameter>> map;
      for (std::size_t i = 0; i < nbParameters; ++i) {
        const auto & param = parameters[i];
        auto value = NumericMutable<double>::create (c, param.getValue ());
        map.emplace (param.getName (),
                     ConfiguredParameter::create (c, {std::move(value)}, param));
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

    NodeRefVec createDependencyVector (const ParameterAliasable & parametrizable,
                                       const std::function<NodeRef (const std::string &)> & getParameter) {
      const auto & parameters = parametrizable.getIndependentParameters ();
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
} // namespace bpp
