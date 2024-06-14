// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../Likelihood/ParametrizablePhyloTree.h"
#include "PhyloTree.h"

using namespace bpp;
using namespace std;

PhyloTree::PhyloTree(bool rooted) :
  AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranch>(rooted),
  name_("")
{}

PhyloTree::PhyloTree(const PhyloTree* tree) :
  AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranch>(tree ? *tree : false),
  name_(tree ? tree->name_ : "")
{}

PhyloTree::PhyloTree(const ParametrizablePhyloTree& tree) :
  AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranch>(tree),
  name_("")
{}

void PhyloTree::resetNodesId()
{
  std::vector<shared_ptr<PhyloNode>> nodes = getAllNodes();
  for (unsigned int i = 0; i < nodes.size(); i++)
  {
    setNodeIndex(nodes[i], i);

    if (hasFather(nodes[i]))
      setEdgeIndex(getEdgeToFather(nodes[i]), i);
  }
}

std::shared_ptr<PhyloNode> PhyloTree::getPhyloNode(const std::string& name) const
{
  vector<shared_ptr<PhyloNode>> vpn = getAllNodes();

  for (auto it:vpn)
  {
    if (it->hasName() && it->getName() == name)
      return it;
  }

  return std::make_shared<PhyloNode>();
}

std::vector<std::string> PhyloTree::getAllLeavesNames() const
{
  vector<string> vn;

  vector<shared_ptr<PhyloNode>> vpn = getAllLeaves();

  for (vector<shared_ptr<PhyloNode>>::const_iterator it = vpn.begin(); it != vpn.end(); it++)
  {
    vn.push_back((*it)->getName());
  }

  return vn;
}

void PhyloTree::scaleTree(shared_ptr<PhyloNode> node, double factor)
{
  vector<shared_ptr<PhyloBranch>> branches = getSubtreeEdges(node);
  for (vector<shared_ptr<PhyloBranch>>::iterator currBranch = branches.begin(); currBranch != branches.end(); currBranch++)
  {
    if ((*currBranch)->hasLength())
      (*currBranch)->setLength((*currBranch)->getLength() * factor);
    else
      throw PhyloBranchPException("PhyloTree::scaleTree : undefined length", (*currBranch).get());
  }
}

void PhyloTree::scaleTree(double factor)
{
  scaleTree(getRoot(), factor);
}

void PhyloTree::pruneTree(std::vector<string> leaves)
{
  vector<shared_ptr<PhyloNode>> vpn = getAllLeaves();

  for (auto& leaf:vpn)
  {
    if (std::find(leaves.begin(), leaves.end(), leaf->getName()) == leaves.end())
    {
      while (leaf && isLeaf(leaf)) // to get rid of internal nodes as leaves
      {
        auto fat = getFatherOfNode(leaf);
        deleteNode(leaf);
        leaf = fat;
      }
    }
  }
}

void PhyloTree::setBranchLengths(double l)
{
  vector<shared_ptr<PhyloBranch>> vpn = getAllEdges();

  for (auto& it: vpn)
  {
    it->setLength(l);
  }
}

Vdouble PhyloTree::getBranchLengths() const
{
  Vdouble vl;

  vector<shared_ptr<PhyloBranch>> vpn = getAllEdges();

  for (auto& it: vpn)
  {
    vl.push_back(it->getLength());
  }
  return vl;
}


PhyloTree& PhyloTree::operator+=(const PhyloTree& phylotree)
{
  vector<shared_ptr<PhyloBranch>> vpn = getAllEdges();

  for (auto& it: vpn)
  {
    uint ei = getEdgeIndex(it);

    if (!phylotree.hasEdge(ei))
      throw Exception("Phylotree::operator+= : argument tree does not have edge " + TextTools::toString(ei));
    if (!it->hasLength() || !phylotree.getEdge(ei)->hasLength())
      throw Exception("Phylotree::operator+= : no summing of branches without length.");

    it->setLength(it->getLength() + phylotree.getEdge(ei)->getLength());
  }

  return *this;
}

PhyloTree& PhyloTree::operator-=(const PhyloTree& phylotree)
{
  vector<shared_ptr<PhyloBranch>> vpn = getAllEdges();

  for (auto& it: vpn)
  {
    uint ei = getEdgeIndex(it);

    if (!phylotree.hasEdge(ei))
      throw Exception("Phylotree::operator+= : argument tree does not have edge " + TextTools::toString(ei));
    if (!it->hasLength() || !phylotree.getEdge(ei)->hasLength())
      throw Exception("Phylotree::operator+= : no summing of branches without length.");

    it->setLength(it->getLength() - phylotree.getEdge(ei)->getLength());
  }

  return *this;
}

PhyloTree& PhyloTree::operator/=(const PhyloTree& phylotree)
{
  vector<shared_ptr<PhyloBranch>> vpn = getAllEdges();

  for (auto& it: vpn)
  {
    uint ei = getEdgeIndex(it);

    if (!phylotree.hasEdge(ei))
      throw Exception("Phylotree::operator/= : argument tree does not have edge " + TextTools::toString(ei));
    if (!it->hasLength() || !phylotree.getEdge(ei)->hasLength())
      throw Exception("Phylotree::operator/= : no summing of branches without length.");

    it->setLength(it->getLength() / phylotree.getEdge(ei)->getLength());
  }

  return *this;
}

PhyloTree& PhyloTree::operator*=(const PhyloTree& phylotree)
{
  vector<shared_ptr<PhyloBranch>> vpn = getAllEdges();

  for (auto& it: vpn)
  {
    uint ei = getEdgeIndex(it);

    if (!phylotree.hasEdge(ei))
      throw Exception("Phylotree::operator/= : argument tree does not have edge " + TextTools::toString(ei));
    if (!it->hasLength() || !phylotree.getEdge(ei)->hasLength())
      throw Exception("Phylotree::operator/= : no summing of branches without length.");

    it->setLength(it->getLength() * phylotree.getEdge(ei)->getLength());
  }

  return *this;
}

void PhyloTree::addSubTree(std::shared_ptr<PhyloNode> phyloNode, const Node& node)
{
  for (int i = 0; i < static_cast<int>(node.getNumberOfSons()); i++)
  {
    const Node& fils = *node[i];

    // the son
    auto soni = std::make_shared<PhyloNode>(fils.hasName() ? fils.getName() : "");
    setNodeIndex(soni, (uint)fils.getId());

    auto propi = fils.getNodePropertyNames();
    for (const auto& prop:propi)
    {
      soni->setProperty(prop, *fils.getNodeProperty(prop));
    }

    auto branchi = std::make_shared<PhyloBranch> ();
    if (fils.hasDistanceToFather())
      branchi->setLength(fils.getDistanceToFather());
    setEdgeIndex(branchi, (uint)fils.getId());

    // the branch to the son
    propi = fils.getBranchPropertyNames();
    for (const auto& prop:propi)
    {
      branchi->setProperty(prop, *fils.getBranchProperty(prop));
    }

    // the link
    createNode(phyloNode, soni, branchi);

    // recursion
    addSubTree(soni, fils);
  }
}
