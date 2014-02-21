//
// File: ParametrizableTree.cpp
// Created by: Julien Dutheil
// Created on: Wed Jul 11 20:36 2012
//

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

#include "ParametrizableTree.h"

using namespace bpp;
using namespace std;

ParametrizableTree::ParametrizableTree(const Tree& tree, bool reparametrizeRoot, bool liveIndex, const std::string& prefix): 
  AbstractParametrizable(prefix),
  tree_(),
  liveIndex_(liveIndex),
  index_(),
  reverseIndex_(),
  isSynchronized_(true),
  minimumBrLen_(0.000001),
  maximumBrLen_(10000),
  brLenConstraint_(new IntervalConstraint(minimumBrLen_, maximumBrLen_, true, true))
{
  if (!dynamic_cast<const TreeTemplate<Node>*>(&tree))
      throw Exception("ParametrizableTree::constructor. Input tree should be TreeTemplate<Node>.");

  tree_=*(dynamic_cast<const TreeTemplate<Node>*>(&tree));
  
  //TODO allow root reparametrization
  if (reparametrizeRoot)
    throw Exception("ParametrizableTree::constructor. Reparametrization of root is not implemented yet.");

  buildIndex_(*tree_.getRootNode()); 
  if (liveIndex_)
    buildReverseIndex_(tree_.getRootNode());
}

ParametrizableTree::ParametrizableTree(const ParametrizableTree& pTree): 
  AbstractParametrizable(pTree),
  tree_(pTree.tree_),
  liveIndex_(),
  index_(pTree.index_),
  reverseIndex_(),
  isSynchronized_(true),
  minimumBrLen_(pTree.minimumBrLen_),
  maximumBrLen_(pTree.maximumBrLen_),
  brLenConstraint_(new IntervalConstraint(pTree.minimumBrLen_, pTree.maximumBrLen_, true, true))
{
  if (liveIndex_)
    buildReverseIndex_(tree_.getRootNode()); 
}

ParametrizableTree& ParametrizableTree::operator=(const ParametrizableTree& pTree)
{
  AbstractParametrizable::operator=(pTree);
  tree_ = pTree.tree_;
  liveIndex_ = pTree.liveIndex_;
  index_ = pTree.index_;
  reverseIndex_.clear();
  if (liveIndex_) {
    buildReverseIndex_(tree_.getRootNode());
  }
  isSynchronized_ = true;
  minimumBrLen_ = pTree.minimumBrLen_;
  maximumBrLen_ = pTree.maximumBrLen_;
  brLenConstraint_.reset(new IntervalConstraint(pTree.minimumBrLen_, pTree.maximumBrLen_, true, true));
  return *this;
}

size_t ParametrizableTree::buildIndex_(Node& node, size_t nPar)
{
  size_t npar = nPar;
  
  if (node.hasFather()) {
    index_[node.getId()] = npar;
  
    double d = minimumBrLen_;
    if (!node.hasDistanceToFather())
    {
      ApplicationTools::displayWarning("Missing branch length " + TextTools::toString(node.getId()) + ". Value is set to " + TextTools::toString(minimumBrLen_));
      node.setDistanceToFather(minimumBrLen_);
    }
    else
    {
      d = node.getDistanceToFather();
      if (d < minimumBrLen_)
      {
        ApplicationTools::displayWarning("Branch length " + TextTools::toString(node.getId()) + " is too small: " + TextTools::toString(d) + ". Value is set to " + TextTools::toString(minimumBrLen_));
        node.setDistanceToFather(minimumBrLen_);
        d = minimumBrLen_;
      }
      if (d > maximumBrLen_)
      {
        ApplicationTools::displayWarning("Branch length " + TextTools::toString(node.getId()) + " is too big: " + TextTools::toString(d) + ". Value is set to " + TextTools::toString(maximumBrLen_));
        node.setDistanceToFather(maximumBrLen_);
        d = maximumBrLen_;
      }
    }
    //if (!hasParameter("BrLen" + TextTools::toString(node.getId()))){
    addParameter_(new Parameter("BrLen" + TextTools::toString(node.getId()), d, brLenConstraint_->clone(), true)); // Attach constraint to avoid clonage problems!
    //}
    npar++;
  }
    
  //Now apply recursively:
  for (unsigned int i = 0; i < node.getNumberOfSons(); ++i)
    npar = buildIndex_(*node[i], npar);

  return npar;
}

void ParametrizableTree::buildReverseIndex_(Node* node)
{
  if (node->hasFather()) {
    reverseIndex_["BrLen" + TextTools::toString(node->getId())] = node; 
  }

  for (unsigned int i = 0; i < node->getNumberOfSons(); ++i)
    buildReverseIndex_((*node)[i]);
}

const TreeTemplate<Node>& ParametrizableTree::getTree() const
{
  if (liveIndex_ && !isSynchronized_)
    updateTreeFromParameters_();
  return tree_;
}

void ParametrizableTree::updateTreeFromParameters_() const
{
  for (unsigned int i = 0; i < getNumberOfParameters(); ++i) {
    const Parameter& param = getParameter_(i);
    reverseIndex_[param.getName()]->setDistanceToFather(param.getValue());
  }
  isSynchronized_ = true;
}

void ParametrizableTree::fireParameterChanged(const ParameterList& parameters)
{
  if (liveIndex_) {
    for (size_t i = 0; i < parameters.size(); ++i) {
      const Parameter& param = parameters[i];
      if (reverseIndex_.find(param.getName()) != reverseIndex_.end())
        reverseIndex_[param.getName()]->setDistanceToFather(param.getValue());
    }
    isSynchronized_ = true;
  } else {
    isSynchronized_ = false;
  }
}

