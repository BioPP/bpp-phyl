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
#include <Bpp/NewPhyl/Utils.h> // TODO for now
#include <algorithm>
#include <stack>

#include <typeindex>
#include <typeinfo>

#include <Bpp/NewPhyl/Config.h>
#ifdef BPP_HAVE_DEMANGLING // TODO replace with sfinae test to remove dep on cmake ?
#include <cstdlib>
#include <cxxabi.h>
#endif
namespace bpp {
std::string demangle (const char * name) {
#ifdef BPP_HAVE_DEMANGLING
	int status{};
	std::unique_ptr<char, void (*) (void *)> res{
	    abi::__cxa_demangle (name, nullptr, nullptr, &status), std::free};
	return status == 0 ? res.get () : name;
#else
	return name;
#endif
}

std::string prettyTypeName (const std::type_info & ti) {
	return demangle (ti.name ());
}
std::string prettyTypeName (std::type_index ti) {
	return demangle (ti.name ());
}
} // namespace bpp

namespace bpp {
namespace dataflow {
	/*****************************************************************************
	 * Error & dependency check functions.
	 */

	// Error functions DataFlow.h
	void failureNodeConversion (const std::type_info & handleType, const Node & node) {
		throw Exception (prettyTypeName (handleType) +
		                 " cannot store: " + prettyTypeName (typeid (node)));
	}

	// Error functions DataFlowTemplates.h
	void failureComputeWasCalled (const std::type_info & nodeType) {
		throw Exception (prettyTypeName (nodeType) + ": compute() was called");
	}

	// Error functions DataFlowInternal.h
	void failureDependencyNumberMismatch (const std::type_info & contextNodeType,
	                                      std::size_t expectedSize, std::size_t givenSize) {
		throw Exception (prettyTypeName (contextNodeType) + ": expected " +
		                 std::to_string (expectedSize) + " dependencies, got " +
		                 std::to_string (givenSize));
	}

	void failureEmptyDependency (const std::type_info & contextNodeType, std::size_t depIndex) {
		throw Exception (prettyTypeName (contextNodeType) + ": " + std::to_string (depIndex) +
		                 "-th dependency is empty (nullptr)");
	}

	void failureDependencyTypeMismatch (const std::type_info & contextNodeType, std::size_t depIndex,
	                                    const std::type_info & expectedType,
	                                    const std::type_info & givenNodeType) {
		throw Exception (prettyTypeName (contextNodeType) + ": expected class derived from " +
		                 prettyTypeName (expectedType) + " as " + std::to_string (depIndex) +
		                 "-th dependency, got " + prettyTypeName (givenNodeType));
	}

	void checkDependencyVectorSize (const std::type_info & contextNodeType, const NodeRefVec & deps,
	                                std::size_t expectedSize) {
		auto size = deps.size ();
		if (size != expectedSize) {
			failureDependencyNumberMismatch (contextNodeType, expectedSize, size);
		}
	}

	void checkDependenciesNotNull (const std::type_info & contextNodeType, const NodeRefVec & deps) {
		for (std::size_t i = 0; i < deps.size (); ++i) {
			if (!deps[i]) {
				failureEmptyDependency (contextNodeType, i);
			}
		}
	}

	/************************************ Base DF::Node impl *******************************/

	Node::Node (const NodeRefVec & dependenciesArg) : dependencyNodes_ (dependenciesArg) {
		for (auto & n : dependencyNodes_)
			n->registerNode (this);
	}
	Node::Node (NodeRefVec && dependenciesArg) : dependencyNodes_ (std::move (dependenciesArg)) {
		for (auto & n : dependencyNodes_)
			n->registerNode (this);
	}

	Node::~Node () {
		for (auto & n : dependencyNodes_)
			n->unregisterNode (this);
	}

	std::string Node::description () const { return prettyTypeName (typeid (*this)); }

	std::string Node::debugInfo () const { return {}; }

	bool Node::hasNumericalProperty (NumericalProperty) const { return false; }

	NodeRef Node::derive (const Node &) {
		throw Exception ("Node does not support derivation: " + description ());
	}
	bool Node::isDerivable (const Node &) const { return false; }

	bool Node::isTransitivelyDependentOn (const Node & node) const {
		std::stack<const Node *> nodesToVisit;
		nodesToVisit.push (this);
		while (!nodesToVisit.empty ()) {
			auto * current = nodesToVisit.top ();
			nodesToVisit.pop ();
			if (&node == current)
				return true;
			for (auto & dep : current->dependencies ())
				nodesToVisit.push (dep.get ());
		}
		return false;
	}

	NodeRef Node::rebuild (NodeRefVec &&) const {
		throw Exception ("Node does not support rebuild(deps): " + description ());
	}

	void Node::invalidateRecursively () noexcept {
		if (!isValid ())
			return;
		std::stack<Node *> nodesToInvalidate;
		nodesToInvalidate.push (this);
		while (!nodesToInvalidate.empty ()) {
			auto * n = nodesToInvalidate.top ();
			nodesToInvalidate.pop ();
			if (n->isValid ()) {
				n->makeInvalid ();
				for (auto * dependent : n->dependentNodes_)
					nodesToInvalidate.push (dependent);
			}
		}
	}

	void Node::computeRecursively () {
		// Compute the current node (and dependencies recursively) if needed
		if (isValid ())
			return;

		// Discover then recompute needed nodes
		std::stack<Node *> nodesToVisit;
		std::stack<Node *> nodesToRecompute;
		nodesToVisit.push (this);
		while (!nodesToVisit.empty ()) {
			auto * n = nodesToVisit.top ();
			nodesToVisit.pop ();
			if (!n->isValid ())
				nodesToRecompute.push (n);
			for (auto & dep : n->dependencies ())
				nodesToVisit.push (dep.get ());
		}
		while (!nodesToRecompute.empty ()) {
			auto * n = nodesToRecompute.top ();
			nodesToRecompute.pop ();
			n->compute ();
			n->makeValid ();
		}
	}

	void Node::registerNode (Node * n) { dependentNodes_.emplace_back (n); }
	void Node::unregisterNode (const Node * n) {
		dependentNodes_.erase (std::remove (dependentNodes_.begin (), dependentNodes_.end (), n),
		                       dependentNodes_.end ());
	}

	NodeRef rebuildWithSubstitution (const NodeRef & node,
	                                 const std::map<const Node *, NodeRef> & substitutions) {
		auto it = substitutions.find (node.get ());
		if (it != substitutions.end ()) {
			// Substitute sub tree.
			return it->second;
		} else if (node->dependencies ().empty ()) {
			// Leaves do no support rebuild: just copy them
			return node;
		} else {
			// Recursion : only rebuild if dependencies have changed
			auto rebuiltDeps = mapToVector (node->dependencies (), [&substitutions](const NodeRef & dep) {
				return rebuildWithSubstitution (dep, substitutions);
			});
			if (rebuiltDeps == node->dependencies ()) {
				return node;
			} else {
				return node->rebuild (std::move (rebuiltDeps));
			}
		}
	}
} // namespace dataflow
} // namespace bpp
