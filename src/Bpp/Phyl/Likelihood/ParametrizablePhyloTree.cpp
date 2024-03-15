// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "ParametrizablePhyloTree.h"

using namespace bpp;
using namespace std;

ParametrizablePhyloTree::ParametrizablePhyloTree(const PhyloTree& tree, const std::string& prefix) :
  AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchParam>(tree),
  AbstractParametrizable(prefix),
  minimumBrLen_(0.000001),
  maximumBrLen_(10000),
  brLenConstraint_(new IntervalConstraint(minimumBrLen_, maximumBrLen_, true, true))
{
  vector<shared_ptr<PhyloBranchParam> > vB = getAllEdges();

  for (auto& it:vB)
  {
    if (hasEdgeIndex(it))
      it->getParameter_(0).setName("BrLen" + TextTools::toString(getEdgeIndex(it)));
    shareParameter_(it->getParameter(0));
  }
}

ParametrizablePhyloTree::ParametrizablePhyloTree(const ParametrizablePhyloTree& pTree) :
  AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchParam>(pTree),
  AbstractParametrizable(pTree),
  minimumBrLen_(pTree.minimumBrLen_),
  maximumBrLen_(pTree.maximumBrLen_),
  brLenConstraint_(pTree.brLenConstraint_ ? pTree.brLenConstraint_->clone() : 0)
{}


ParametrizablePhyloTree& ParametrizablePhyloTree::operator=(const ParametrizablePhyloTree& pTree)
{
  AssociationTreeGlobalGraphObserver<PhyloNode, PhyloBranchParam>::operator=(pTree);
  AbstractParametrizable::operator=(pTree);
  minimumBrLen_ = pTree.minimumBrLen_;
  maximumBrLen_ = pTree.maximumBrLen_;
  if (pTree.brLenConstraint_)
    brLenConstraint_.reset(pTree.brLenConstraint_->clone());
  else
    brLenConstraint_ = 0;


  return *this;
}

std::vector<std::string> ParametrizablePhyloTree::getAllLeavesNames() const
{
  vector<string> vn;

  vector<shared_ptr<PhyloNode> > vpn = getAllLeaves();

  for (vector<shared_ptr<PhyloNode> >::const_iterator it = vpn.begin(); it != vpn.end(); it++)
  {
    vn.push_back((*it)->getName());
  }

  return vn;
}

Vdouble ParametrizablePhyloTree::getBranchLengths() const
{
  Vdouble vl;

  vector<shared_ptr<PhyloBranchParam> > vpn = getAllEdges();

  for (auto& it: vpn)
  {
    vl.push_back(it->getLength());
  }
  return vl;
}

void ParametrizablePhyloTree::fireParameterChanged(const ParameterList& parameters)
{}
