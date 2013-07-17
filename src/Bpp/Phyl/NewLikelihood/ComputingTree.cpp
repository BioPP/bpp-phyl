//
// File: ComputingTree.cpp
// Created by: Laurent Guéguen
// Created on: mardi 9 juillet 2013, à 15h 37
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

#include "ComputingTree.h"

#include <Bpp/Exceptions.h>

using namespace bpp;
using namespace std;

ComputingTree::ComputingTree(const ComputingTree& tree) :
  AbstractParametrizable(tree),
  parTree_(tree.parTree_),
  pDist_(tree.pDist_),
  vTree_(),
  vvBrMod_(tree.vvBrMod_)
{
  for (size_t i=0; i<tree.vTree_.size(); i++)
    vTree_.push_back(tree.vTree_[i]->clone());
}

ComputingTree& ComputingTree::operator=(const ComputingTree& tree)
{
  AbstractParametrizable::operator=(tree);
  parTree_=tree.parTree_;
  pDist_=tree.pDist_;
  
  for (size_t i=0; i<tree.vTree_.size(); i++)
    vTree_.push_back(tree.vTree_[i]->clone());

  vvBrMod_=tree.vvBrMod_;
  
  return *this;
}

ComputingTree::~ComputingTree()
{
  for (size_t i=0;i<vTree_.size();i++)
    TreeTemplateTools::deleteSubtree(dynamic_cast<ComputingNode*>(vTree_[i]->getRootNode()));

  vTree_.clear();
}

ComputingTree::ComputingTree(const ParametrizableTree& ptree, const DiscreteDistribution& dist) :
  AbstractParametrizable(""),
  parTree_(&ptree),
  pDist_(&dist),
  vTree_(),
  vvBrMod_()
{
  TreeTemplate<Node> tree=ptree.getTree();
  
  ComputingNode* rCN= TreeTemplateTools::cloneSubtree<ComputingNode>(*tree.getRootNode());
  TreeTemplate<ComputingNode>* pTC=new TreeTemplate<ComputingNode>(rCN);

  size_t nbCl=dist.getNumberOfCategories();
  
  for (size_t i=0; i<nbCl; i++){
    TreeTemplate<ComputingNode>* pTC2=pTC->clone();
    pTC2->scaleTree(dist.getCategory(i));
    vTree_.push_back(pTC2);
  }
  
  delete pTC;

  addParameters_(ptree.getParameters());
  addParameters_(dist.getParameters());
}

void ComputingTree::addModel(const SubstitutionModel* pSubMod, std::vector<int>  vBr)
{
  for (size_t i=0; i< getNumberOfClasses(); i++)
    for (size_t j=0; j<vBr.size(); j++){
      vTree_[i]->getNode(vBr[j])->setSubstitutionModel(pSubMod);
    }

  vvBrMod_.push_back(vBr);
  
  addParameters_(pSubMod->getParameters());
}


void ComputingTree::fireParameterChanged(const ParameterList& pl)
{
  if (vTree_.size()==0)
    return;

  // In this case, all trees have the same models on the nodes
  
  for (size_t i=0; i<vvBrMod_.size(); i++){
    if (vTree_[0]->getNode(vvBrMod_[i][0])->getParameters().testParametersValues(pl)){
      for (size_t j=0; j<vvBrMod_[i].size(); j++)
        for (size_t k=0; k<vTree_.size(); k++)
          vTree_[k]->getNode(vvBrMod_[i][j])->fireParameterChanged(pl);
    }
  }

  // test the 
}

