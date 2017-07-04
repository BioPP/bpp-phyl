//
// File: ImportNewlik.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-06-29
// Last modified: 2017-06-29
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
#include <Bpp/NewPhyl/ImportNewlik.h>
#include <Bpp/Phyl/Tree/PhyloTree.h>

namespace bpp {
namespace Phyl {
	ConvertedPhyloTreeData convertPhyloTree (const bpp::PhyloTree & fromTree) {
		if (!fromTree.isRooted ())
			throw Exception ("PhyloTree is not rooted");

		// Build topology : create nodes and an index conversion map
		// Assumes our indexes are densely allocated (in sequence from 0)
		Topology::IndexMapBase<IndexType> fromNodeIdMap{
		    static_cast<SizeType> (fromTree.getNumberOfNodes ())};
		auto convertId = [&fromNodeIdMap](const bpp::PhyloTree::NodeIndex & index) {
			return fromNodeIdMap.index (static_cast<IndexType> (index)).value ();
		};
		auto tmpTree = make_freezable<Topology::Tree> ();
		auto allNodeIds = fromTree.getAllNodesIndexes ();
		for (auto i : allNodeIds) {
			auto ourId = tmpTree->createNode ();
			fromNodeIdMap.set (ourId, static_cast<IndexType> (i));
		}
		for (auto i : allNodeIds) {
			if (fromTree.hasFather (i)) {
				auto fatherId = fromTree.getNodeIndex (fromTree.getFather (fromTree.getNode (i)));
				tmpTree->createEdge (convertId (fatherId), convertId (i));
			}
		}
		tmpTree->setRootNodeId (convertId (fromTree.getRootIndex ()));
		auto tree = std::move (tmpTree).freeze ();

		// Data
		auto brLens = make_freezable<Topology::BranchValueMap<double>> (tree);
		auto nodeNames = make_freezable<Topology::NodeIndexMap<std::string>> (tree);
		for (auto i : allNodeIds) {
			auto node = tree->node (convertId (i));
			// Branch length
			if (fromTree.hasFather (i)) {
				auto branch = fromTree.getEdgeToFather (i);
				brLens->access (node.fatherBranch ()) = branch->getLength ();
			}
			// Leaf name
			auto phyloNode = fromTree.getNode (i);
			if (phyloNode->hasName ())
				nodeNames->set (node, phyloNode->getName ());
		}

		return {std::move (tree),
		        make_frozen<Topology::NodeIndexMap<IndexType>> (tree, std::move (fromNodeIdMap)),
		        std::move (brLens).freeze (), std::move (nodeNames).freeze ()};
	}
}
}
