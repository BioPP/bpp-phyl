//
// File: Topology.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-19
// Last modified: 2017-05-19
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

#include <Bpp/NewPhyl/Topology.h>
#include <Bpp/NewPhyl/TopologyAnnotation.h>
#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <memory>
#include <stdexcept>
#include <unordered_map>

namespace bpp {
namespace Topology {
	// Convert from PhyloTree
	ConvertedPhyloTreeData convertPhyloTree (const bpp::PhyloTree & phyloTree) {
		ConvertedPhyloTreeData returnData{};
		{
			auto tree = Tree::create ();

			// Build topology
			std::unordered_map<bpp::PhyloTree::NodeIndex, IndexType> phyloNodeIdToOurIds;
			for (auto phyloNodeId : phyloTree.getAllNodesIndexes ()) {
				// Create all nodes
				auto ourId = tree->createNode ();
				phyloNodeIdToOurIds[phyloNodeId] = ourId;
			}
			for (auto phyloNodeId : phyloTree.getAllNodesIndexes ()) {
				if (phyloTree.hasFather (phyloNodeId)) {
					// Link them by using the father link
					auto phyloFatherId =
					    phyloTree.getNodeIndex (phyloTree.getFather (phyloTree.getNode (phyloNodeId)));
					tree->createEdge (phyloNodeIdToOurIds.at (phyloFatherId),
					                  phyloNodeIdToOurIds.at (phyloNodeId));
				}
			}
			if (!phyloTree.isRooted ())
				throw std::runtime_error ("PhyloTree is not rooted");
			tree->setRootNodeId (phyloNodeIdToOurIds.at (phyloTree.getRootIndex ()));

			returnData.topology = Tree::finalize (std::move (tree));
		}
		return returnData;
	}
}
}
