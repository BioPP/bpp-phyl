//
// File: Debug.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-28 00:00:00
// Last modified: 2017-04-28
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

#include <Bpp/NewPhyl/Config.h>
#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/FrozenPtr.h>
#include <Bpp/NewPhyl/NodeSpecification.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/NewPhyl/Topology.h>
#include <algorithm>
#include <memory>
#include <ostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <unordered_set>

#ifdef BPP_HAVE_DEMANGLING
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

namespace Topology {
	// Print tree structure
	void debugTree (std::ostream & os, FrozenSharedPtr<Tree> tree) {
		os << "digraph {\n";
		std::queue<Node> nodesToVisit;
		if (tree->rootNodeId () != invalid)
			nodesToVisit.emplace (Node (tree, tree->rootNodeId ()));

		while (!nodesToVisit.empty ()) {
			auto node = nodesToVisit.front ();
			nodesToVisit.pop ();
			os << '\t' << node.nodeId () << " [shape=box,label=\"" << node.nodeId () << "\"];\n";
			node.foreachChildBranch ([&os, &node, &nodesToVisit](Branch && branch) {
				auto childNode = std::move (branch).childNode ();
				os << '\t' << node.nodeId () << " -> " << childNode.nodeId () << ";\n";
				nodesToVisit.emplace (std::move (childNode));
			});
		}
		os << "}\n";
	}
}
namespace DF {
	namespace {
		// Dot utils
		std::string dotLabelEscape (std::string s) {
			// Escapes characters in a record type dot node label.
			const char toEscape[] = "<>|{} ";
			std::string result;
			for (auto c : s) {
				if (std::any_of (std::begin (toEscape), std::end (toEscape),
				                 [c](char c2) { return c == c2; }))
					result.push_back ('\\');
				result.push_back (c);
			}
			return result;
		}
		std::string typeToDotLabel (const std::type_index & type) {
			return dotLabelEscape (demangle (type.name ()));
		}

		// Dot node key: hash code reduced to a short
		std::string dotNodeKey (char type, std::size_t hash) {
			return type + std::to_string (std::uint16_t (hash));
		}
		std::string dotNodeKey (const Node::Impl * p) {
			return dotNodeKey ('N', std::hash<const Node::Impl *>{}(p));
		}
		std::string dotNodeKey (const Node & n) { return dotNodeKey (&n.getImpl ()); }
		std::string dotNodeKey (const Registry::Key & key) { return dotNodeKey ('K', key.hashCode ()); }
		std::string dotNodeKey (const NodeSpecification & nodeSpec) {
			return dotNodeKey ('S', std::hash<std::string>{}(nodeSpec.description ()));
		}

		// Dot pretty print of elements
		void dotNodePretty (std::ostream & os, const Node::Impl * node) {
			os << '\t' << dotNodeKey (node) << " [color=blue,shape=record,label=\"" << dotNodeKey (node)
			   << '|' << typeToDotLabel (typeid (*node)) << "\"];\n";
		}
		void dotNodePretty (std::ostream & os, const Registry::Key & key) {
			os << '\t' << dotNodeKey (key) << " [shape=Mrecord,label=\"{" << dotNodeKey (key) << "|{"
			   << typeToDotLabel (key.operation ()) << '|';
			for (auto & dep : key.dependencies ())
				os << dotNodeKey (dep) << ' ';
			os << "}}\"];\n";
		}
		void dotNodePretty (std::ostream & os, const NodeSpecification & spec) {
			os << '\t' << dotNodeKey (spec) << " [color=red,shape=record,label=\"{" << dotNodeKey (spec)
			   << "|" << dotLabelEscape (spec.description ()) << "}\"];\n";
		}

		template <typename T, typename U>
		void dotEdgePretty (std::ostream & os, const T & from, const U & to,
		                    const std::string & style) {
			os << '\t' << dotNodeKey (from) << " -> " << dotNodeKey (to) << ' ' << style << ";\n";
		}
		void dotEdgePretty (std::ostream & os, const Node::Impl * from, const Node::Impl * to) {
			dotEdgePretty (os, from, to, "[color=blue]");
		}
		void dotEdgePretty (std::ostream & os, const Registry::Key & from, const Node & to) {
			dotEdgePretty (os, from, to, "");
		}
		void dotEdgePretty (std::ostream & os, const NodeSpecification & from,
		                    const NodeSpecification & to) {
			dotEdgePretty (os, from, to, "[color=red]");
		}
		void dotEdgePretty (std::ostream & os, const NodeSpecification & from, const Node & to) {
			dotEdgePretty (os, from, to, "[color=green]");
		}

		// Print the DF dag structure (in blue).
		// Takes list of entry points.
		void debugDagStructure (std::ostream & os, std::vector<Node> entryPoints) {
			std::queue<const Node::Impl *> nodesToVisit;
			std::unordered_set<const Node::Impl *> nodesAlreadyVisited;

			for (auto & n : entryPoints)
				nodesToVisit.emplace (&n.getImpl ());

			while (!nodesToVisit.empty ()) {
				auto node = nodesToVisit.front ();
				nodesToVisit.pop ();
				if (nodesAlreadyVisited.count (node))
					continue;

				dotNodePretty (os, node);
				nodesAlreadyVisited.emplace (node);

				node->foreachDependentNode ([&](const Node::Impl * p) {
					if (nodesAlreadyVisited.count (p))
						return;
					dotEdgePretty (os, p, node);
					nodesToVisit.emplace (p);
				});
				node->foreachDependencyNode ([&](const Node::Impl * p) {
					if (nodesAlreadyVisited.count (p))
						return;
					dotEdgePretty (os, node, p);
					nodesToVisit.emplace (p);
				});
			}
		}

		// Print registry keys, and links to stored nodes (key only).
		// Returns list of pointed-to nodes.
		std::vector<Node> debugRegistryLinks (std::ostream & os, const Registry & registry) {
			std::vector<Node> entryPoints;
			registry.foreachKeyValue ([&entryPoints, &os](const Registry::Key & key, const Node & node) {
				dotNodePretty (os, key);
				dotEdgePretty (os, key, node);
				entryPoints.emplace_back (node);
			});
			return entryPoints;
		}

		// Instantiate a NodeSpec (without registry), duplicate of NodeSpec.instantiate
		// Print NodeSpec details, and links between NodeSpecs.
		// Print links to node (key only).
		Node debugPlayNodeSpecInstantiation (std::ostream & os, const NodeSpecification & nodeSpec) {
			NodeVec deps;
			dotNodePretty (os, nodeSpec);
			for (auto & depSpec : nodeSpec.computeDependencies ()) {
				deps.emplace_back (debugPlayNodeSpecInstantiation (os, depSpec));
				dotEdgePretty (os, nodeSpec, depSpec);
			}
			auto n = nodeSpec.buildNode (std::move (deps));
			dotEdgePretty (os, nodeSpec, n);
			return n;
		}

		// Replay instantiation of a NodeSpec (with registry), must have already be build in registry.
		// Print NodeSpec details, and links between NodeSpecs.
		// Print links to node (key only).
		Node debugReplayNodeSpecInstantiationInRegistry (std::ostream & os,
		                                                 const NodeSpecification & nodeSpec,
		                                                 const Registry & registry) {
			dotNodePretty (os, nodeSpec);
			auto depSpecs = nodeSpec.computeDependencies ();
			if (depSpecs.empty ()) {
				Node n = nodeSpec.buildNode ({});
				dotEdgePretty (os, nodeSpec, n);
				return n;
			} else {
				// Instantiate dependencies
				NodeVec deps;
				for (auto & depSpec : depSpecs) {
					deps.emplace_back (debugReplayNodeSpecInstantiationInRegistry (os, depSpec, registry));
					dotEdgePretty (os, nodeSpec, depSpec);
				}
				// Check the registry
				auto result = registry.get (Registry::Key{nodeSpec.nodeType (), deps});
				if (!result)
					throw std::runtime_error (
					    "debugReplayNodeSpecInstantiationWithReuse: nodeSpec was not found in registry");
				dotEdgePretty (os, nodeSpec, *result);
				return std::move (result).value ();
			}
		}
	}

	void debugDag (std::ostream & os, const Node & entryPoint) {
		os << "digraph {\n";
		debugDagStructure (os, {entryPoint});
		os << "}\n";
	}
	void debugRegistry (std::ostream & os, const Registry & registry) {
		os << "digraph {\n";
		debugDagStructure (os, debugRegistryLinks (os, registry));
		os << "}\n";
	}
	void debugNodeSpecInstantiation (std::ostream & os, const NodeSpecification & nodeSpec) {
		os << "digraph {\n";
		debugDagStructure (os, {debugPlayNodeSpecInstantiation (os, nodeSpec)});
		os << "}\n";
	}
	void debugNodeSpecInstantiationInRegistry (std::ostream & os, const NodeSpecification & nodeSpec,
	                                           const Registry & registry, bool showRegistryLinks) {
		os << "digraph {\n";
		debugDagStructure (os, {debugReplayNodeSpecInstantiationInRegistry (os, nodeSpec, registry)});
		if (showRegistryLinks)
			debugRegistryLinks (os, registry);
		os << "}\n";
	}
}
}
