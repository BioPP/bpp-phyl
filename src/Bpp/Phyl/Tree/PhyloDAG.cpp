// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

// #include "../Likelihood/ParametrizablePhyloDAG.h"
#include "PhyloDAG.h"

using namespace bpp;
using namespace std;

PhyloDAG::PhyloDAG() :
  AssociationDAGlobalGraphObserver<PhyloNode, PhyloBranch>(),
  name_("")
{}

PhyloDAG::PhyloDAG(const PhyloDAG* dag) :
  AssociationDAGlobalGraphObserver<PhyloNode, PhyloBranch>(*dag),
  name_(dag ? dag->name_ : "")
{}

// PhyloDAG::PhyloDAG(const ParametrizablePhyloDAG& dag) :
//   AssociationDAGlobalGraphObserver<PhyloNode, PhyloBranch>(dag),
//   name_("")
// {}

void PhyloDAG::resetNodesId()
{
  std::vector<shared_ptr<PhyloNode>> nodes = getAllNodes();
  std::vector<shared_ptr<PhyloBranch>> branches = getAllEdges();

  for (unsigned int i = 0; i < nodes.size(); i++)
  {
    setNodeIndex(nodes[i], i);

    // set arbitrary edge node nmber to an incoming edge
    if (hasFather(nodes[i]))
    {
      const auto vInEdges = getIncomingEdges(nodes[i]);
      setEdgeIndex(vInEdges[0], i);
    }
  }

  size_t eid = 0;
  for (auto j = nodes.size(); j < branches.size(); j++)
  {
    while (hasEdgeIndex(branches[eid]))
      eid++;
    setEdgeIndex(branches[eid], (uint)j);
  }
}

std::shared_ptr<PhyloNode> PhyloDAG::getPhyloNode(const std::string& name) const
{
  vector<shared_ptr<PhyloNode>> vpn = getAllNodes();

  for (auto it:vpn)
  {
    if (it->hasName() && it->getName() == name)
      return it;
  }

  return std::make_shared<PhyloNode>();
}

std::vector<std::string> PhyloDAG::getAllLeavesNames() const
{
  vector<string> vn;

  vector<shared_ptr<PhyloNode>> vpn = getAllLeaves();

  for (vector<shared_ptr<PhyloNode>>::const_iterator it = vpn.begin(); it != vpn.end(); it++)
  {
    vn.push_back((*it)->getName());
  }

  return vn;
}

void PhyloDAG::scaleDAG(shared_ptr<PhyloNode> node, double factor)
{
  vector<shared_ptr<PhyloBranch>> branches = getBelowEdges(node);
  for (vector<shared_ptr<PhyloBranch>>::iterator currBranch = branches.begin(); currBranch != branches.end(); currBranch++)
  {
    if ((*currBranch)->hasLength())
      (*currBranch)->setLength((*currBranch)->getLength() * factor);
    else
      throw PhyloBranchPException("PhyloDAG::scaleDAG : undefined length", (*currBranch).get());
  }
}

void PhyloDAG::scaleDAG(double factor)
{
  scaleDAG(getRoot(), factor);
}

void PhyloDAG::pruneDAG(std::vector<string> leaves)
{
  vector<shared_ptr<PhyloNode>> vpn = getAllLeaves();

  for (auto& leaf:vpn)
  {
    if (std::find(leaves.begin(), leaves.end(), leaf->getName()) == leaves.end())
    {
      std::vector<shared_ptr<PhyloNode>> vfat({leaf});  // one vector  per leaf removed to avoid too large vector

      std::vector<shared_ptr<PhyloNode>>::iterator vfatit;  // one vector  per leaf removed to avoid too large vector

      for (vfatit = vfat.begin(); vfatit != vfat.end();)
      {
        if (hasNode(*vfatit) && isLeaf(*vfatit))
        {
          auto vfat2 = getFathers(*vfatit);
          vfat.insert(vfat.end(), vfat2.begin(), vfat2.end() );
          deleteNode(*vfatit);
        }
      }
    }
  }
}

void PhyloDAG::setBranchLengths(double l)
{
  vector<shared_ptr<PhyloBranch>> vpn = getAllEdges();

  for (auto& it: vpn)
  {
    it->setLength(l);
  }
}

Vdouble PhyloDAG::getBranchLengths() const
{
  Vdouble vl;

  vector<shared_ptr<PhyloBranch>> vpn = getAllEdges();

  for (auto& it: vpn)
  {
    vl.push_back(it->getLength());
  }
  return vl;
}


PhyloDAG& PhyloDAG::operator+=(const PhyloDAG& phylodag)
{
  vector<shared_ptr<PhyloBranch>> vpn = getAllEdges();

  for (auto& it: vpn)
  {
    uint ei = getEdgeIndex(it);

    if (!phylodag.hasEdge(ei))
      throw Exception("Phylodag::operator+= : argument dag does not have edge " + TextTools::toString(ei));
    if (!it->hasLength() || !phylodag.getEdge(ei)->hasLength())
      throw Exception("Phylodag::operator+= : no summing of branches without length.");

    it->setLength(it->getLength() + phylodag.getEdge(ei)->getLength());
  }

  return *this;
}

PhyloDAG& PhyloDAG::operator-=(const PhyloDAG& phylodag)
{
  vector<shared_ptr<PhyloBranch>> vpn = getAllEdges();

  for (auto& it: vpn)
  {
    uint ei = getEdgeIndex(it);

    if (!phylodag.hasEdge(ei))
      throw Exception("Phylodag::operator+= : argument dag does not have edge " + TextTools::toString(ei));
    if (!it->hasLength() || !phylodag.getEdge(ei)->hasLength())
      throw Exception("Phylodag::operator+= : no summing of branches without length.");

    it->setLength(it->getLength() - phylodag.getEdge(ei)->getLength());
  }

  return *this;
}

PhyloDAG& PhyloDAG::operator/=(const PhyloDAG& phylodag)
{
  vector<shared_ptr<PhyloBranch>> vpn = getAllEdges();

  for (auto& it: vpn)
  {
    uint ei = getEdgeIndex(it);

    if (!phylodag.hasEdge(ei))
      throw Exception("Phylodag::operator/= : argument dag does not have edge " + TextTools::toString(ei));
    if (!it->hasLength() || !phylodag.getEdge(ei)->hasLength())
      throw Exception("Phylodag::operator/= : no summing of branches without length.");

    it->setLength(it->getLength() / phylodag.getEdge(ei)->getLength());
  }

  return *this;
}

PhyloDAG& PhyloDAG::operator*=(const PhyloDAG& phylodag)
{
  vector<shared_ptr<PhyloBranch>> vpn = getAllEdges();

  for (auto& it: vpn)
  {
    uint ei = getEdgeIndex(it);

    if (!phylodag.hasEdge(ei))
      throw Exception("Phylodag::operator/= : argument dag does not have edge " + TextTools::toString(ei));
    if (!it->hasLength() || !phylodag.getEdge(ei)->hasLength())
      throw Exception("Phylodag::operator/= : no summing of branches without length.");

    it->setLength(it->getLength() * phylodag.getEdge(ei)->getLength());
  }

  return *this;
}
