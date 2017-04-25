//
// File: DataFlow.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-19
// Last modified: 2017-04-19
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

#include <Bpp/NewPhyl/DataFlow.h>
#include <ostream>
#include <queue>
#include <typeinfo>
#include <unordered_set>

namespace bpp {
namespace DF {

	static std::uintptr_t debugFormat (const Node::Impl * p) { return std::uintptr_t (p); }

	void debugDagStructure (std::ostream & os, const Node & entryPoint) {
		std::queue<const Node::Impl *> nodesToVisit;
		std::unordered_set<const Node::Impl *> nodesAlreadyVisited;

		nodesToVisit.push (&entryPoint.get ());
		while (!nodesToVisit.empty ()) {
			auto node = nodesToVisit.front ();
			os << '\t' << debugFormat (node) << " [shape=box,label=\"" << debugFormat (node) << '-'
			   << typeid (*node).name () << "\"];\n";
			nodesAlreadyVisited.emplace (node);
			nodesToVisit.pop ();

			node->foreachDependentNode ([&](const Node::Impl * p) {
				if (nodesAlreadyVisited.count (p))
					return;
				os << '\t' << debugFormat (p) << " -> " << debugFormat (node) << ";\n";
				nodesToVisit.emplace (p);
			});
			node->foreachDependencyNode ([&](const Node::Impl * p) {
				if (nodesAlreadyVisited.count (p))
					return;
				os << '\t' << debugFormat (node) << " -> " << debugFormat (p) << ";\n";
				nodesToVisit.emplace (p);
			});
		}
	}
	void debugDag (std::ostream & os, const Node & entryPoint) {
		os << "digraph {\n";
		debugDagStructure (os, entryPoint);
		os << "}\n";
	}
}
}
