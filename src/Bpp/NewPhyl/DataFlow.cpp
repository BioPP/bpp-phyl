//
// File: DataFlow.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-30
// Last modified: 2017-05-30
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
#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowTemplates.h>
#include <Bpp/NewPhyl/Debug.h>
#include <typeinfo>

namespace bpp {
namespace DF {
	// Error functions TODO clean names
	void parameterFailComputeWasCalled (const std::type_info & ti) {
		throw Exception (prettyTypeName (ti) + ": compute() was called on a Parameter node");
	}

	void nodeHandleConversionFailed (const std::type_info & handle, const Node::Impl & impl) {
		throw Exception (prettyTypeName (handle) +
		                 " handle type cannot store: " + prettyTypeName (typeid (impl)));
	}

	void genericFunctionComputationCheckDependencyNum (const std::type_info & ti,
	                                                   std::size_t expected, std::size_t given) {
		if (expected != given)
			throw Exception (prettyTypeName (ti) + ": expected " + std::to_string (expected) +
			                 " dependencies, got " + std::to_string (given));
	}

	void dataFlowTemplatesDependencyTypeMismatch (const std::type_info & computationNode,
	                                              std::size_t index, const std::type_info & expected,
	                                              const Node::Impl & given) {
		throw Exception (prettyTypeName (computationNode) + ": expected " + prettyTypeName (expected) +
		                 " as " + std::to_string (index) + "-th argument, got " +
		                 prettyTypeName (typeid (given)));
	}
}
}
