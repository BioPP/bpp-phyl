//
// File: ParametrizablePhyloTree.cpp
// Authors:
//   Laurent GuÃÂ©guen
// Created: 2012-07-11 20:36:00
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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
    shareParameter_(it->getSharedParameter(0));
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
