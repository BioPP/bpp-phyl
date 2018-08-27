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

#include <algorithm>
#include <fstream>     // debug
#include <functional>  // std::hash
#include <ostream>     // debug
#include <stack>       // invalidate/compute recursively + debug
#include <type_traits> // DotOptions flags
#include <typeinfo>
#include <unordered_set> // debug

#include "DataFlow.h"

/* std::type_info::name() returns a "mangled" type name, not very readable.
 * Compilers can optionally provide an ABI header cxxabi.h.
 * This header contain the demangle function, which makes the type name readable.
 * By testing on godbolt, all versions of gcc and clang have this header.
 */
#if defined(__GNUC__) || defined(__clang__)
#include <cstdlib>
#include <cxxabi.h>
static std::string demangle (const char * name) {
	int status{};
	std::unique_ptr<char, void (*) (void *)> res{
	    abi::__cxa_demangle (name, nullptr, nullptr, &status), std::free};
	return status == 0 ? res.get () : name;
}
#else
static std::string demangle (const char * name) {
	return name;
}
#endif

namespace bpp {
std::string prettyTypeName (const std::type_info & ti) {
	return demangle (ti.name ());
}
} // namespace bpp

namespace bpp {
namespace dataflow {
	/*****************************************************************************
	 * Error & dependency check functions.
	 */
	void failureComputeWasCalled (const std::type_info & nodeType) {
		throw Exception (prettyTypeName (nodeType) + ": compute() was called");
	}

	void failureNodeConversion (const std::type_info & handleType, const Node & node) {
		throw Exception (prettyTypeName (handleType) +
		                 " cannot store: " + prettyTypeName (typeid (node)));
	}

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

	/*****************************************************************************
	 * Node impl.
	 */
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

	std::string Node::description () const {
		auto replaceAll = [](std::string & str, const std::string & pattern,
		                     const std::string & replacement) {
			std::string::size_type i = str.find (pattern);
			while (i != std::string::npos) {
				str.replace (i, pattern.size (), replacement);
				i = str.find (pattern, i + replacement.size ());
			}
		};
		std::string nodeType = prettyTypeName (typeid (*this));
    // Shorten displayed name by removing namespaces.
		replaceAll (nodeType, "bpp::dataflow::", "");
		replaceAll (nodeType, "Eigen::", "");
		replaceAll (nodeType, "std::", "");
		return nodeType;
	}

	std::string Node::debugInfo () const { return {}; }

	bool Node::hasNumericalProperty (NumericalProperty) const { return false; }

	NodeRef Node::derive (Context &, const Node &) {
		throw Exception ("Node does not support derivation: " + description ());
	}
	bool Node::isDerivable (const Node &) const { return false; }

	NodeRef Node::rebuild (NodeRefVec &&) const {
		throw Exception ("Node does not support rebuild(deps): " + description ());
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
			if (!n->isValid ()) {
				nodesToRecompute.push (n);
				for (auto & dep : n->dependencies ())
					nodesToVisit.push (dep.get ());
			}
		}
		while (!nodesToRecompute.empty ()) {
			auto * n = nodesToRecompute.top ();
			nodesToRecompute.pop ();
			n->compute ();
			n->makeValid ();
		}
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

	void Node::registerNode (Node * n) { dependentNodes_.emplace_back (n); }
	void Node::unregisterNode (const Node * n) {
		dependentNodes_.erase (std::remove (dependentNodes_.begin (), dependentNodes_.end (), n),
		                       dependentNodes_.end ());
	}

	/*****************************************************************************
	 * Free functions.
	 */
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
			NodeRefVec rebuiltDeps (node->nbDependencies ());
			for (std::size_t i = 0; i < rebuiltDeps.size (); ++i) {
				rebuiltDeps[i] = rebuildWithSubstitution (node->dependency (i), substitutions);
			}
			if (rebuiltDeps == node->dependencies ()) {
				return node;
			} else {
				return node->rebuild (std::move (rebuiltDeps));
			}
		}
	}

	bool isTransitivelyDependentOn (const Node & searchedDependency, const Node & node) {
		std::stack<const Node *> nodesToVisit;
		nodesToVisit.push (&node);
		while (!nodesToVisit.empty ()) {
			auto * current = nodesToVisit.top ();
			nodesToVisit.pop ();
			if (current == &searchedDependency)
				return true;
			for (auto & dep : current->dependencies ())
				nodesToVisit.push (dep.get ());
		}
		return false;
	}

	/*****************************************************************************
	 * output dataflow graph to dot format.
	 */

	// Make enum behave as flags.
	DotOptions operator| (DotOptions a, DotOptions b) {
		using IntType = typename std::underlying_type<DotOptions>::type;
		return static_cast<DotOptions> (static_cast<IntType> (a) | static_cast<IntType> (b));
	}
	bool operator& (DotOptions a, DotOptions b) {
		using IntType = typename std::underlying_type<DotOptions>::type;
		return static_cast<IntType> (a) & static_cast<IntType> (b);
	}

	// Escape text in a dot box label
	static std::string dotLabelEscape (const std::string & s) {
		std::string result;
		result.reserve (s.size ());
		static const char toEscape[] = "<>|{} ";
		for (const char c : s) {
			if (std::any_of (std::begin (toEscape), std::end (toEscape),
			                 [c](const char c2) { return c == c2; })) {
				result.push_back ('\\');
			}
			result.push_back (c);
		}
		return result;
	}

	// Dot node id for dataflow Node: 'N' + hash as a string
	static std::string dotIdentifier (const Node & node) {
		return 'N' + std::to_string (std::hash<const Node *>{}(&node));
	}

	// Write line with node representation
	static void writeDotNode (std::ostream & os, const Node & node, DotOptions opt) {
		os << '\t' << dotIdentifier (node);
		if (opt & DotOptions::DetailedNodeInfo) {
			os << " [shape=Mrecord,label=\"{" << dotLabelEscape (node.description ())
			   << "| valid=" << node.isValid () << ' ' << dotLabelEscape (node.debugInfo ()) << "}\"]";
		} else {
			os << " [shape=box,label=\"" << dotLabelEscape (node.description ()) << "\"]";
		}
		os << ";\n";
	}

	// Write line with edge representation for n-th dependency of from
	static void writeDotEdge (std::ostream & os, const Node & from, std::size_t depIndex,
	                          DotOptions opt) {
		os << '\t' << dotIdentifier (from) << " -> " << dotIdentifier (*from.dependency (depIndex));
		if (opt & DotOptions::ShowDependencyIndex) {
			os << " [label=\"" << depIndex << "\"]";
		}
		os << ";\n";
	}

	// Write dot lines for graph structure, starting from the given entry points.
	static void writeGraphStructure (std::ostream & os, const std::vector<const Node *> & entryPoints,
	                                 DotOptions opt) {
		std::stack<const Node *> nodesToVisit;
		std::unordered_set<const Node *> discoveredNodes;

		const auto discover = [&nodesToVisit, &discoveredNodes](const Node * n) {
			const bool discovered = discoveredNodes.find (n) != discoveredNodes.end ();
			if (!discovered) {
				nodesToVisit.emplace (n);
				discoveredNodes.emplace (n);
			}
		};

		for (const auto * n : entryPoints) {
			discover (n);
		}

		while (!nodesToVisit.empty ()) {
			const auto * node = nodesToVisit.top ();
			nodesToVisit.pop ();
			writeDotNode (os, *node, opt);
			if (opt & DotOptions::FollowUpwardLinks) {
				for (const auto * dependent : node->dependentNodes ()) {
					discover (dependent);
				}
			}
			for (std::size_t index = 0; index < node->nbDependencies (); ++index) {
				writeDotEdge (os, *node, index, opt);
				discover (node->dependency (index).get ());
			}
		}
	}

	void writeGraphToDot (std::ostream & os, const std::vector<const Node *> & nodes,
	                      DotOptions opt) {
		os << "digraph {\n";
		writeGraphStructure (os, nodes, opt);
		os << "}\n";
	}
	void writeGraphToDot (const std::string & filename, const std::vector<const Node *> & nodes,
	                      DotOptions opt) {
		std::ofstream file{filename};
		writeGraphToDot (file, nodes, opt);
	}
} // namespace dataflow
} // namespace bpp
