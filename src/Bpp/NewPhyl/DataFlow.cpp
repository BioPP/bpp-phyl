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
#include <algorithm>
#include <stack>
#include <typeinfo>

namespace bpp {
namespace DF {
	// Error functions
	void failureParameterComputeWasCalled (const std::type_info & paramType) {
		throw Exception (prettyTypeName (paramType) + ": compute() was called on a Parameter node");
	}

	void failureNodeHandleConversion (const std::type_info & handleType, const Node::Impl & node) {
		throw Exception (prettyTypeName (handleType) +
		                 " cannot store: " + prettyTypeName (typeid (node)));
	}

	static void failureDependencyNumberMismatch (const std::type_info & computeNodeType,
	                                             SizeType expectedSize, SizeType givenSize) {
		throw Exception (prettyTypeName (computeNodeType) + ": expected " +
		                 std::to_string (expectedSize) + " dependencies, got " +
		                 std::to_string (givenSize));
	}
	void checkDependencyNumber (const std::type_info & computeNodeType, SizeType expectedSize,
	                            SizeType givenSize) {
		if (expectedSize != givenSize)
			failureDependencyNumberMismatch (computeNodeType, expectedSize, givenSize);
	}

	void failureDependencyTypeMismatch (const std::type_info & computeNodeType, IndexType depIndex,
	                                    const std::type_info & expectedType,
	                                    const Node::Impl & givenNode) {
		throw Exception (prettyTypeName (computeNodeType) + ": expected class derived from " +
		                 prettyTypeName (expectedType) + " as " + std::to_string (depIndex) +
		                 "-th dependency, got " + prettyTypeName (typeid (givenNode)));
	}

	//

	void Node::Impl::invalidate () noexcept {
		if (isValid ()) {
			isValid_ = false;
			for (auto * impl : dependentNodes_)
				impl->invalidate ();
		}
	}

	void Node::Impl::computeRecursively () {
		// Compute the current node (and dependencies recursively) if needed
		if (isValid ())
			return;

		// Discover then recompute needed nodes
		std::stack<Impl *> nodesToVisit;
		std::stack<Impl *> nodesToRecompute;
		nodesToVisit.push (this);
		while (!nodesToVisit.empty ()) {
			auto * n = nodesToVisit.top ();
			nodesToVisit.pop ();
			if (!n->isValid ())
				nodesToRecompute.push (n);
			for (auto & dep : n->dependencies ())
				nodesToVisit.push (&dep.getImpl ());
		}
		while (!nodesToRecompute.empty ()) {
			auto * n = nodesToRecompute.top ();
			nodesToRecompute.pop ();
			n->compute ();
			n->makeValid ();
		}
	}

	void Node::Impl::registerNode (Impl * n) { dependentNodes_.emplace_back (n); }
	void Node::Impl::unregisterNode (const Impl * n) {
		dependentNodes_.erase (std::remove (dependentNodes_.begin (), dependentNodes_.end (), n),
		                       dependentNodes_.end ());
	}
}
}
