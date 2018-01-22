//
// File: ImportNewlik.h
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

#ifndef BPP_NEWPHYL_IMPORTNEWLIK_H
#define BPP_NEWPHYL_IMPORTNEWLIK_H

#include <Bpp/NewPhyl/Phylogeny.h>
#include <Bpp/Phyl/Tree/PhyloTree.h> // FIXME forward declare (use node/edge types ?)

namespace bpp {

// Forward declaration
class PhyloTree;

class PhyloTreeView : public TreeTopologyInterface,
                      public BranchLengthValueAccess,
                      public SequenceNameValueAccess {
public:
	PhyloTreeView (const PhyloTree & tree) : tree_ (tree) {}

	// Convert between TreeTopologyInterface indexes and PhyloTree indexes
	static inline PhyloTree::NodeIndex convertNode (TopologyNodeIndex id) {
		return static_cast<PhyloTree::NodeIndex> (id.value);
	}
	static inline TopologyNodeIndex convertNode (PhyloTree::NodeIndex id) {
		return TopologyNodeIndex (static_cast<std::intptr_t> (id));
	}
	static inline PhyloTree::EdgeIndex convertBranch (TopologyBranchIndex id) {
		return static_cast<PhyloTree::EdgeIndex> (id.value);
	}
	static inline TopologyBranchIndex convertBranch (PhyloTree::EdgeIndex id) {
		return TopologyBranchIndex (static_cast<std::intptr_t> (id));
	}

	// TreeTopologyInterface
	TopologyNodeIndex rootNode () const final;
	TopologyNodeIndex fatherNode (TopologyBranchIndex id) const final;
	TopologyNodeIndex childNode (TopologyBranchIndex id) const final;
	TopologyBranchIndex fatherBranch (TopologyNodeIndex id) const final;
	Vector<TopologyBranchIndex> childBranches (TopologyNodeIndex id) const final;

	// BranchLengthValueAccess
	double getBranchLengthValue (TopologyBranchIndex id) const final;

	// SequenceNameValueAccess
	std::string getSequenceName (TopologyNodeIndex id) const final;

private:
	const PhyloTree & tree_;
};
} // namespace bpp

#endif // BPP_NEWPHYL_IMPORTNEWLIK_H
