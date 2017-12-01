//
// File: ImportMaster.cpp
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
#include <Bpp/NewPhyl/ImportMaster.h>
#include <Bpp/Phyl/Tree/Node.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>

namespace bpp {

/* In TreeTemplate there is no notion of edges, only nodes.
 * In this view class, the index of a branch is defined as the index of its child node.
 */
TreeTopologyView::NodeIndex TreeTemplateView::rootNode () const {
	// isRooted method checks if 2 sons, fails on newphyl example that has 3 sons...
	if (tree_.getRootNode () == nullptr)
		throw Exception ("TreeTemplateView: tree has no root node");
	return convert (tree_.getRootId ());
}
TreeTopologyView::NodeIndex TreeTemplateView::fatherNode (BranchIndex id) const {
	return convert (tree_.getFatherId (convert (childNode (id))));
}
TreeTopologyView::NodeIndex TreeTemplateView::childNode (BranchIndex id) const {
	return NodeIndex (id.value);
}
TreeTopologyView::BranchIndex TreeTemplateView::fatherBranch (NodeIndex id) const {
	return BranchIndex (id.value);
}
Vector<TreeTopologyView::BranchIndex> TreeTemplateView::childBranches (NodeIndex id) const {
	return mapToVector (tree_.getSonsId (convert (id)),
	                    [](int i) { return BranchIndex (static_cast<std::intptr_t> (i)); });
}

double TreeTemplateView::getBranchLengthValue (BranchIndex id) const {
	return tree_.getDistanceToFather (convert (childNode (id)));
}

std::string TreeTemplateView::getSequenceName (NodeIndex id) const {
	return tree_.getNodeName (convert (id));
}
} // namespace bpp
