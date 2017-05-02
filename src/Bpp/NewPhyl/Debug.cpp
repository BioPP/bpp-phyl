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
#include <Bpp/NewPhyl/Registry.h>
#include <Bpp/NewPhyl/Topology.h>
#include <algorithm>
#include <ostream>
#include <queue>
#include <string>
#include <typeinfo>
#include <unordered_set>

#ifdef BPP_HAVE_DEMANGLING
#include <cstdlib>
#include <cxxabi.h>
#include <memory>
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

namespace Topology {
	void debugTree (std::ostream & os, const Tree & tree) {
		os << "digraph {\n";

		std::queue<IndexType> nodesToVisit;
		if (tree.rootId () != invalid)
			nodesToVisit.emplace (tree.rootId ());

		while (!nodesToVisit.empty ()) {
			auto nodeId = nodesToVisit.front ();
			nodesToVisit.pop ();
			auto & node = tree.node (nodeId);
			os << '\t' << nodeId << " [shape=box,label=\"" << nodeId << '-' << node.nodeName_ << "\"];\n";
			for (auto childId : node.childrenIds_) {
				nodesToVisit.emplace (childId);
				os << '\t' << nodeId << " -> " << childId << ";\n";
			}
		}

		os << "}\n";
	}
}
namespace DF {
	// Dot node key: hash code reduced to a short
	static std::string dotNodeKey (char type, std::size_t hash) {
		return type + std::to_string (std::uint16_t (hash));
	}
	static std::string dotNodeKey (const Node::Impl * p) {
		return dotNodeKey ('N', std::hash<const Node::Impl *>{}(p));
	}
	static std::string dotNodeKey (const RegistryKey & key) {
		return dotNodeKey ('R', key.hashCode ());
	}

	static std::string dotLabelEscape (std::string s) {
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
	static std::string typeToDotLabel (const std::type_index & type) {
		return dotLabelEscape (demangle (type.name ()));
	}

	static void debugDagStructure (std::ostream & os, const Node & entryPoint) {
		std::queue<const Node::Impl *> nodesToVisit;
		std::unordered_set<const Node::Impl *> nodesAlreadyVisited;

		nodesToVisit.emplace (&entryPoint.getImpl ());
		while (!nodesToVisit.empty ()) {
			auto node = nodesToVisit.front ();
			os << '\t' << dotNodeKey (node) << " [shape=record,label=\"" << dotNodeKey (node) << '|'
			   << typeToDotLabel (typeid (*node)) << "\"];\n";
			nodesAlreadyVisited.emplace (node);
			nodesToVisit.pop ();

			node->foreachDependentNode ([&](const Node::Impl * p) {
				if (nodesAlreadyVisited.count (p))
					return;
				os << '\t' << dotNodeKey (p) << " -> " << dotNodeKey (node) << ";\n";
				nodesToVisit.emplace (p);
			});
			node->foreachDependencyNode ([&](const Node::Impl * p) {
				if (nodesAlreadyVisited.count (p))
					return;
				os << '\t' << dotNodeKey (node) << " -> " << dotNodeKey (p) << ";\n";
				nodesToVisit.emplace (p);
			});
		}
	}

	static void debugRegistryLinks (std::ostream & os, const Registry & registry) {
		for (auto & it : registry.rawAccess ()) {
			auto & key = it.first;
			os << '\t' << dotNodeKey (key) << " [shape=Mrecord,label=\"{" << dotNodeKey (key) << "|{";
			os << typeToDotLabel (key.operation ()) << '|';
			if (key.element ().type () == bpp::Topology::Element::Node)
				os << 'N' << key.element ().asNodeRef ().nodeId ();
			else
				os << 'B' << key.element ().asBranchRef ().childNodeId ();
			os << "}}\"];\n";
			os << '\t' << dotNodeKey (key) << " -> " << dotNodeKey (&it.second.getImpl ()) << ";\n";
		}
	}

	void debugDag (std::ostream & os, const Node & entryPoint) {
		os << "digraph {\n";
		debugDagStructure (os, entryPoint);
		os << "}\n";
	}
	void debugRegistry (std::ostream & os, const Registry & registry) {
		os << "digraph {\n";
		if (!registry.rawAccess ().empty ()) {
			debugDagStructure (os, registry.rawAccess ().begin ()->second);
			debugRegistryLinks (os, registry);
		}
		os << "}\n";
	}
}
}
