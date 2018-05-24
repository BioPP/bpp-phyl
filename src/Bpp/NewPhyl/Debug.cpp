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
#include <Bpp/NewPhyl/IntegerRange.h>
#include <Bpp/NewPhyl/Phylogeny.h>
#include <algorithm>
#include <fstream>
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

namespace {
	// Escape text in a dot box label
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

	// Generic dot node id generation: a char + a hash code as a uint16_t to keep it short
	std::string dotNodeKey (char type, std::size_t hash) {
		return type + std::to_string (std::uint16_t (hash));
	}

	// Craft node keys
	std::string dotNodeKey (TopologyNodeIndex id) { return 'n' + std::to_string (id.value); }
	std::string dotNodeKey (const DF::Node * p) {
		return dotNodeKey ('N', std::hash<const DF::Node *>{}(p));
	}
	std::string dotNodeKey (const DF::NamedNodeRef & namedNode) {
		return dotNodeKey ('T', std::hash<std::string>{}(namedNode.name));
	}

	// Generic pretty print of edge (with a style argument)
	template <typename T, typename U>
	void dotEdgePretty (std::ostream & os, const T & from, const U & to, const std::string & style) {
		os << '\t' << dotNodeKey (from) << " -> " << dotNodeKey (to) << ' ' << style << ";\n";
	}

	// Pretty print nodes
	void dotNodePretty (std::ostream & os, TopologyNodeIndex id) {
		os << '\t' << dotNodeKey (id) << " [shape=box,label=\"" << id.value << "\"];\n";
	}
	void dotNodePretty (std::ostream & os, const DF::Node * node) {
		os << '\t' << dotNodeKey (node) << " [color=blue,shape=record,label=\"" << dotNodeKey (node)
		   << '|' << dotLabelEscape (node->description ()) << "\"];\n";
	}
	void dotNodePrettyDetailed (std::ostream & os, const DF::Node * node) {
		auto debugInfo = (node->isValid () ? "valid " : "invalid ") + node->debugInfo ();
		os << '\t' << dotNodeKey (node) << " [color=blue,shape=Mrecord,label=\"{{" << dotNodeKey (node)
		   << '|' << dotLabelEscape (node->description ()) << "}|" << dotLabelEscape (debugInfo)
		   << "}\"];\n";
	}
	void dotNodePretty (std::ostream & os, const DF::NamedNodeRef & namedNode) {
		os << '\t' << dotNodeKey (namedNode) << " [color=orange,shape=record,label=\""
		   << dotLabelEscape (namedNode.name) << "\"];\n";
	}

	// Pretty print edges
	void dotEdgePretty (std::ostream & os, TopologyNodeIndex parent, TopologyBranchIndex branch,
	                    TopologyNodeIndex child) {
		dotEdgePretty (os, parent, child, "[label=\"" + std::to_string (branch.value) + "\"]");
	}
	void dotEdgePretty (std::ostream & os, const DF::Node * from, const DF::Node * to) {
		dotEdgePretty (os, from, to, "[color=blue]");
	}
	void dotEdgePretty (std::ostream & os, const DF::Node * from, const DF::Node * to,
	                    std::size_t num) {
		dotEdgePretty (os, from, to, "[color=blue,label=\"" + std::to_string (num) + "\"]");
	}
	void dotEdgePretty (std::ostream & os, const DF::NamedNodeRef & from, const DF::Node * to) {
		dotEdgePretty (os, from, to, "[color=orange]");
	}
} // namespace

// Print tree structure
void debugTree (std::ostream & os, const TreeTopologyInterface & tree) {
	os << "digraph {\n";
	std::queue<TopologyNodeIndex> nodesToVisit;
	nodesToVisit.emplace (tree.rootNode ());
	while (!nodesToVisit.empty ()) {
		auto nodeId = nodesToVisit.front ();
		nodesToVisit.pop ();
		dotNodePretty (os, nodeId);
		for (auto branch : tree.childBranches (nodeId)) {
			auto child = tree.childNode (branch);
			dotEdgePretty (os, nodeId, branch, child);
			nodesToVisit.emplace (child);
		}
	}
	os << "}\n";
}

// Print the DF dag structure (in blue).
// Takes list of entry points.
void debugDagStructure (std::ostream & os, std::vector<const DF::Node *> entryPoints,
                        DF::DebugOptions opt) {
	std::queue<const DF::Node *> nodesToVisit;
	std::unordered_set<const DF::Node *> nodesAlreadyVisited;

	for (const auto * n : entryPoints)
		nodesToVisit.emplace (n);

	while (!nodesToVisit.empty ()) {
		auto * node = nodesToVisit.front ();
		nodesToVisit.pop ();
		if (nodesAlreadyVisited.count (node))
			continue;

		if (opt & DF::DebugOptions::DetailedNodeInfo) {
			dotNodePrettyDetailed (os, node);
		} else {
			dotNodePretty (os, node);
		}
		nodesAlreadyVisited.emplace (node);

		if (opt & DF::DebugOptions::FollowUpwardLinks)
			for (auto * p : node->dependentNodes ())
				if (!nodesAlreadyVisited.count (p))
					nodesToVisit.emplace (p);

		for (auto index : range (node->dependencies ().size ())) {
			auto * dep = node->dependencies ()[index].get ();
			if (opt & DF::DebugOptions::ShowDependencyIndex) {
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
std::vector<const DF::Node *>
debugNamedNodeRefs (std::ostream & os, const std::vector<DF::NamedNodeRef> & namedNodes) {
	std::vector<const DF::Node *> nodes;
	for (auto & namedNodeRef : namedNodes) {
		dotNodePretty (os, namedNodeRef);
		dotEdgePretty (os, namedNodeRef, namedNodeRef.nodeRef.get ());
		nodes.emplace_back (namedNodeRef.nodeRef.get ());
	}
	return nodes;
}

void debugDag (std::ostream & os, const DF::Node & entryPoint, DF::DebugOptions opt) {
	os << "digraph {\n";
	debugDagStructure (os, {&entryPoint}, opt);
	os << "}\n";
}
void debugDag (const std::string & filename, const DF::Node & entryPoint, DF::DebugOptions opt) {
	std::ofstream file{filename};
	debugDag (file, entryPoint, opt);
}

void debugDag (std::ostream & os, const std::vector<DF::NamedNodeRef> & namedNodes,
               DF::DebugOptions opt) {
	os << "digraph {\n";
	debugDagStructure (os, debugNamedNodeRefs (os, namedNodes), opt);
	os << "}\n";
}
void debugDag (const std::string & filename, const std::vector<DF::NamedNodeRef> & namedNodes,
               DF::DebugOptions opt) {
	std::ofstream file{filename};
	debugDag (file, namedNodes, opt);
}
} // namespace bpp

#ifndef NDEBUG
/// Simpler version of debugDag, for use in gdb interpreter (only compiled in Debug mode)
void bpp_df_print_dag (const bpp::DF::Node & node) {
	bpp::debugDag ("df_debug", node, bpp::DF::DebugOptions::DetailedNodeInfo);
}
#endif
