//
// File: ImportMaster.h
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

#ifndef BPP_NEWPHYL_IMPORTMASTER_H
#define BPP_NEWPHYL_IMPORTMASTER_H

#include <Bpp/NewPhyl/FrozenPtr.h>
#include <Bpp/NewPhyl/Phylogeny.h>
#include <Bpp/NewPhyl/Signed.h>
#include <Bpp/NewPhyl/Topology.h>
#include <Bpp/NewPhyl/TopologyMap.h>

namespace bpp {

// Forward declarations
class Node;
template <typename N> class TreeTemplate;
class VectorSiteContainer;

class TreeTemplateView : public TreeTopologyView,
                         public BranchLengthValueAccess,
                         public SequenceNameValueAccess {
public:
	TreeTemplateView (const TreeTemplate<Node> & tree) : tree_ (tree) {}

	// Convert between TreeTopologyView node Index and TreeTemplate node id
	static inline int convert (NodeIndex id) { return static_cast<int> (id.value); }
	static inline NodeIndex convert (int nodeId) {
		return NodeIndex (static_cast<std::intptr_t> (nodeId));
	}

	// TreeTopologyView
	bool valid (NodeIndex id) const override final;
	bool valid (BranchIndex id) const override final;
	NodeIndex rootNode () const override final;
	NodeIndex fatherNode (BranchIndex id) const override final;
	NodeIndex childNode (BranchIndex id) const override final;
	BranchIndex fatherBranch (NodeIndex id) const override final;
	Vector<BranchIndex> childBranches (NodeIndex id) const override final;

	// BranchLengthValueAccess
	double getBranchLengthValue (BranchIndex id) const override final;

  // SequenceNameValueAccess
  std::string getSequenceName (NodeIndex id) const override final;

private:
	const TreeTemplate<Node> & tree_;
};

namespace Phyl {
	// Extract information from a TreeTemplate<Node> structure (master)
	struct ConvertedTreeTemplateData {
		const FrozenPtr<Topology::Tree> topology;
		const FrozenPtr<Topology::NodeIndexMap<IndexType>> treeTemplateNodeIndexes;
		const FrozenPtr<Topology::BranchValueMap<double>> branchLengths;
		const FrozenPtr<Topology::NodeIndexMap<std::string>> nodeNames;
	};
	ConvertedTreeTemplateData convertTreeTemplate (const TreeTemplate<Node> & tree);

	// Create sequence map from bpp VectorSiteContainer and a name mapping.
	SequenceMap makeSequenceMap (const Topology::NodeIndexMap<std::string> & nodeNames,
	                             const VectorSiteContainer & sequences);
} // namespace Phyl
} // namespace bpp

#endif // BPP_NEWPHYL_IMPORTMASTER_H
