// File: PhylogenyTree.h
// Authors:
//   Francois Gindraud (2017)
// Created: 15/03/2017

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _PHYLOGENY_TREE_H_
#define _PHYLOGENY_TREE_H_

//#include <Bpp/Phyl/NewLikelihood/SubstitutionProcess.h>

#include <Bpp/Phyl/DF/DataFlowBaseClasses.h>
#include <Bpp/Phyl/DF/ForRange.h>

#include <cstdint>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

namespace bpp
{
  // Forward declare manipulator
  template <typename NodeType, typename BranchType, typename StorageType>
  class TreeManipulator;

  // Base class of extended trees.
  // Contains a tree topology only.
  class ExtendableTree
  {
  public:
    using IndexType = std::uint32_t;
    static constexpr IndexType invalidIndex = std::numeric_limits<IndexType>::max();

    // Base node type
    class Node
    {
    private:
      IndexType id_;
      IndexType fatherBranch_{invalidIndex};
      std::vector<IndexType> childBranches_;

    public:
      virtual ~Node() = default;
    };
    // Base branch type
    class Branch
    {
    private:
      IndexType id_;
      IndexType fatherNode_;
      IndexType childNode_;

    public:
      virtual ~Branch() = default;
    };
    // Store tree.
    // Ptr to node and branches to allow virtual dispatching.
    // Also because tree elements will ultimately store DF nodes that are non movable.
    struct Storage
    {
      std::vector<std::unique_ptr<Node>> nodes_;
      std::vector<std::unique_ptr<Branch>> branches_;
    };

  protected:
    Storage tree_;

  public:
    // Method that allow to access the tree.
    // Should be overriden in each extended tree subclass to return a manipulator with extended Node and Branch.
    TreeManipulator<Node, Branch, Storage> tree(void);
    TreeManipulator<const Node, const Branch, const Storage> tree(void) const;
  };

  // Class that manipulates the tree (access, modification, etc).
  // Access the extensions of basic Node and Branch given in the templates.
  template <typename NodeType, typename BranchType, typename StorageType>
  class TreeManipulator
  {
  public:
    using IndexType = ExtendableTree::IndexType;
    static_assert(std::is_base_of<ExtendableTree::Node, NodeType>::value,
                  "NodeType must inherit from ExtendableTree::Node");
    static_assert(std::is_base_of<ExtendableTree::Branch, BranchType>::value,
                  "BranchType must inherit from ExtendableTree::Branch");
    static_assert(std::is_same<ExtendableTree::Storage, typename std::decay<StorageType>::type>::value,
                  "StorageType must be a const or non-const ExtendableTree::Storage");

  private:
    StorageType& tree_;

  public:
    TreeManipulator(StorageType& storage)
      : tree_(storage)
    {
    }

    NodeType& node(IndexType index) const { return *tree_.nodes_[index]; }
    BranchType& branch(IndexType index) const { return *tree_.branches_[index]; }

    IndexType addNode(void) const { return 0; }
  };

  // Need to put it after TreeManipulator definition
  auto ExtendableTree::tree(void) -> TreeManipulator<Node, Branch, Storage> { return {tree_}; }
  auto ExtendableTree::tree(void) const -> TreeManipulator<const Node, const Branch, const Storage> { return {tree_}; }

  class PhylogenyTree : public virtual ExtendableTree
  {
  public:
    class Node : public virtual ExtendableTree::Node
    {
    };
    class Branch : public virtual ExtendableTree::Branch
    {
    };

    TreeManipulator<Node, Branch, Storage> tree(void) { return {tree_}; }
  };

  class PhylogenyProcess : public virtual PhylogenyTree
  {
    // Add a model per branch
  public:
    class Node : public virtual PhylogenyTree::Node
    {
    };
    class Branch : public virtual PhylogenyTree::Branch
    {
    };

    TreeManipulator<Node, Branch, Storage> tree(void) { return {tree_}; }
  };

} // end of namespace bpp.

#endif //_PHYLOGENY_TREE_H_
