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

namespace bpp {

namespace {
	// Generic dot node id generation: a char + a hash code as a uint16_t to keep it short
	std::string dotNodeKey (char type, std::size_t hash) {
		return type + std::to_string (std::uint16_t (hash));
	}

	// Craft node keys
	std::string dotNodeKey (TopologyNodeIndex id) { return 'n' + std::to_string (id.value); }
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
	void dotNodePretty (std::ostream & os, const DF::NamedNodeRef & namedNode) {
		os << '\t' << dotNodeKey (namedNode) << " [color=orange,shape=record,label=\""
		   << dotLabelEscape (namedNode.name) << "\"];\n";
	}

	// Pretty print edges
	void dotEdgePretty (std::ostream & os, TopologyNodeIndex parent, TopologyBranchIndex branch,
	                    TopologyNodeIndex child) {
		dotEdgePretty (os, parent, child, "[label=\"" + std::to_string (branch.value) + "\"]");
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

