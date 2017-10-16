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
	void debugTree (std::ostream & os, FrozenPtr<Tree> tree) {
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
} // namespace Topology
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
		std::string dotNodeKey (const Node * p) {
			return dotNodeKey ('N', std::hash<const Node *>{}(p));
		}
		std::string dotNodeKey (const NodeRef & p) { return dotNodeKey (p.get ()); }
		// std::string dotNodeKey (const Registry::Key & key) { return dotNodeKey ('K', key.hashCode
		// ()); }
		std::string dotNodeKey (const NamedNodeRef & namedNode) {
			return dotNodeKey ('T', std::hash<std::string>{}(namedNode.name));
		}

		// Dot pretty print of elements
		void dotNodePretty (std::ostream & os, const Node * node) {
			os << '\t' << dotNodeKey (node) << " [color=blue,shape=record,label=\"" << dotNodeKey (node)
			   << '|' << dotLabelEscape (node->description ()) << "\"];\n";
		}
		void dotNodePrettyDetailed (std::ostream & os, const Node * node) {
			auto debugInfo = (node->isValid () ? "valid " : "invalid ") + node->debugInfo ();
			os << '\t' << dotNodeKey (node) << " [color=blue,shape=Mrecord,label=\"{{"
			   << dotNodeKey (node) << '|' << dotLabelEscape (node->description ()) << "}|"
			   << dotLabelEscape (debugInfo) << "}\"];\n";
		}
		/*void dotNodePretty (std::ostream & os, const Registry::Key & key) {
		  os << '\t' << dotNodeKey (key) << " [shape=Mrecord,label=\"{" << dotNodeKey (key) << "|{"
		     << typeToDotLabel (key.operation ()) << '|';
		  for (auto & dep : key.dependencies ())
		    os << dotNodeKey (dep) << ' ';
		  os << "}}\"];\n";
		}*/
		void dotNodePretty (std::ostream & os, const NamedNodeRef & namedNode) {
			os << '\t' << dotNodeKey (namedNode) << " [color=orange,shape=record,label=\""
			   << dotLabelEscape (namedNode.name) << "\"];\n";
		}

		template <typename T, typename U>
		void dotEdgePretty (std::ostream & os, const T & from, const U & to,
		                    const std::string & style) {
			os << '\t' << dotNodeKey (from) << " -> " << dotNodeKey (to) << ' ' << style << ";\n";
		}
		void dotEdgePretty (std::ostream & os, const Node * from, const Node * to) {
			dotEdgePretty (os, from, to, "[color=blue]");
		}
		void dotEdgePretty (std::ostream & os, const Node * from, const Node * to, SizeType num) {
			dotEdgePretty (os, from, to, "[color=blue,label=\"" + std::to_string (num) + "\"]");
		}
		/*void dotEdgePretty (std::ostream & os, const Registry::Key & from, const Node * to) {
		  dotEdgePretty (os, from, to, "");
		}*/
		void dotEdgePretty (std::ostream & os, const NamedNodeRef & from, const Node * to) {
			dotEdgePretty (os, from, to, "[color=orange]");
		}

		// Print the DF dag structure (in blue).
		// Takes list of entry points.
		void debugDagStructure (std::ostream & os, NodeRefVec entryPoints, DebugOptions opt) {
			std::queue<const Node *> nodesToVisit;
			std::unordered_set<const Node *> nodesAlreadyVisited;

			for (auto & n : entryPoints)
				nodesToVisit.emplace (n.get ());

			while (!nodesToVisit.empty ()) {
				auto * node = nodesToVisit.front ();
				nodesToVisit.pop ();
				if (nodesAlreadyVisited.count (node))
					continue;

				if (opt & DebugOptions::DetailedNodeInfo) {
					dotNodePrettyDetailed (os, node);
				} else {
					dotNodePretty (os, node);
				}
				nodesAlreadyVisited.emplace (node);

				if (opt & DebugOptions::FollowUpwardLinks)
					for (auto * p : node->dependentNodes ())
						if (!nodesAlreadyVisited.count (p))
							nodesToVisit.emplace (p);

				for (auto index : index_range (node->dependencies ())) {
					auto * dep = node->dependencies ()[index].get ();
					if (opt & DebugOptions::ShowDependencyIndex) {
						dotEdgePretty (os, node, dep, index);
					} else {
						dotEdgePretty (os, node, dep);
					}
					if (!nodesAlreadyVisited.count (dep))
						nodesToVisit.emplace (dep);
				}
			}
		}

		// Print name tags pointing to nodes.
		// Returns list of nodes.
		NodeRefVec debugNamedNodeRefs (std::ostream & os, const Vector<NamedNodeRef> & namedNodes) {
			NodeRefVec nodes;
			for (auto & namedNodeRef : namedNodes) {
				dotNodePretty (os, namedNodeRef);
				dotEdgePretty (os, namedNodeRef, namedNodeRef.nodeRef.get ());
				nodes.emplace_back (namedNodeRef.nodeRef);
			}
			return nodes;
		}

		// Print registry keys, and links to stored nodes (key only).
		// Returns list of pointed-to nodes.
		/*NodeRefVec debugRegistryLinks (std::ostream & os, const Registry & registry) {
		  NodeRefVec entryPoints;
		  registry.foreachKeyValue (
		      [&entryPoints, &os](const Registry::Key & key, const NodeRef & node) {
		        dotNodePretty (os, key);
		        dotEdgePretty (os, key, node.get ());
		        entryPoints.emplace_back (node);
		      });
		  return entryPoints;
		}*/
	} // namespace

	void debugDag (std::ostream & os, const std::shared_ptr<Node> & entryPoint, DebugOptions opt) {
		os << "digraph {\n";
		debugDagStructure (os, {entryPoint}, opt);
		os << "}\n";
	}
	void debugDag (std::ostream & os, const Vector<NamedNodeRef> & namedNodes, DebugOptions opt) {
		os << "digraph {\n";
		debugDagStructure (os, debugNamedNodeRefs (os, namedNodes), opt);
		os << "}\n";
	}
	/*void debugRegistry (std::ostream & os, const Registry & registry, DebugOptions opt) {
	  os << "digraph {\n";
	  debugDagStructure (os, debugRegistryLinks (os, registry), opt);
	  os << "}\n";
	}*/
} // namespace DF
} // namespace bpp
