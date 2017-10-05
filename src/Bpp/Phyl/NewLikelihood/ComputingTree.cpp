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

#include "SubstitutionProcessCollection.h"


#include <Bpp/Exceptions.h>

using namespace bpp;
using namespace std;

ComputingTree::ComputingTree(const ParametrizablePhyloTree& ptree, const DiscreteDistribution& dist) :
  AbstractParametrizable(""),  
  parTree_(&ptree),
  pDist_(&dist),
  vTree_(),
  isReadyToCompute_(false)
{
  size_t nbCl=dist.getNumberOfCategories();

  for (size_t i=0; i<nbCl; i++){

    std::shared_ptr<SpecCompTree> pTC2(new SpecCompTree(ptree));

    std::vector<unsigned int> vCN=pTC2->getAllNodesIndexes();
    
    for (size_t j=0; j<vCN.size(); j++){
      pTC2->getNode(vCN[j])->setParameterValue("scale",dist.getCategory(i));
      if (pTC2->hasFather(vCN[j]))
        pTC2->getNode(vCN[j])->setDistanceToFather(ptree.getEdgeToFather(vCN[j])->getLength());
    }

    std::shared_ptr<CompTree> pTC(new CompTree(*pTC2));

    for (size_t j=0; j<vCN.size(); j++){
      pTC->getNode(vCN[j])->updateTree(pTC.get(), vCN[j]);
    }
    
    vTree_.push_back(pTC);
  }
  
  addParameters_(ptree.getParameters());
  addParameters_(dist.getIndependentParameters());
}

ComputingTree::ComputingTree(const ParametrizablePhyloTree& ptree) :
  AbstractParametrizable(""),  
  parTree_(&ptree),
  pDist_(0),
  vTree_(),
  isReadyToCompute_(false)
{
  std::shared_ptr<SpecCompTree> pTC2(new SpecCompTree(ptree));

  std::vector<unsigned int> vCN=pTC2->getAllNodesIndexes();
    
  for (size_t j=0; j<vCN.size(); j++){
    pTC2->getNode(vCN[j])->setParameterValue("scale",1.);
    if (pTC2->hasFather(vCN[j]))
      pTC2->getNode(vCN[j])->setDistanceToFather(ptree.getEdgeToFather(vCN[j])->getLength());
  }

  std::shared_ptr<CompTree> pTC(new CompTree(*pTC2));

  for (size_t j=0; j<vCN.size(); j++){
    pTC->getNode(vCN[j])->updateTree(pTC.get(), vCN[j]);
  }

  vTree_.push_back(pTC);

  addParameters_(ptree.getParameters());
}

ComputingTree::ComputingTree(const SubstitutionProcessCollection* pSubProColl, size_t nTree, size_t nDist) :
  AbstractParametrizable(""),  
  parTree_(&pSubProColl->getTree(nTree)),
  pDist_(&pSubProColl->getRateDistribution(nDist)),
  vTree_(),
  isReadyToCompute_(false)
{
  size_t nbCl=pDist_->getNumberOfCategories();

  for (size_t i=0; i<nbCl; i++){

    std::shared_ptr<SpecCompTree> pTC2(new SpecCompTree(*parTree_));
    std::vector<unsigned int> vCN=pTC2->getAllNodesIndexes();
    
    for (size_t j=0; j<vCN.size(); j++){
      pTC2->getNode(vCN[j])->setParameterValue("scale",pDist_->getCategory(i));
      if (pTC2->hasFather(vCN[j]))
        pTC2->getNode(vCN[j])->setDistanceToFather(parTree_->getEdgeToFather(vCN[j])->getLength());
    }

    std::shared_ptr<CompTree> pTC(new CompTree(*pTC2));

    for (size_t j=0; j<vCN.size(); j++){
      pTC->getNode(vCN[j])->updateTree(pTC.get(), vCN[j]);
    }
    
    vTree_.push_back(pTC);
  }
  
  ParameterList pl=pSubProColl->getTree(nTree).getParameters();

  for (size_t i=0; i<pl.size(); i++)
    pl[i].setName(pl[i].getName()+"_"+TextTools::toString(nTree));

  addParameters_(pl);
  
  pl=pSubProColl->getRateDistribution(nDist).getParameters();

  for (size_t i=0; i<pl.size(); i++)
    pl[i].setName(pl[i].getName()+"_"+TextTools::toString(nDist));

  addParameters_(pl);
}



ComputingTree::ComputingTree(const ComputingTree& tree) :
  AbstractParametrizable(tree),
  parTree_(tree.parTree_),
  pDist_(tree.pDist_),
  vTree_(),
  isReadyToCompute_(false)
{
  for (size_t i=0; i<tree.vTree_.size(); i++){
    
    std::shared_ptr<CompTree> pTC(tree.vTree_[i]->clone());

    std::vector<std::shared_ptr<ComputingNode> > vCN2=pTC->getAllNodes();
    for (size_t j=0; j<vCN2.size(); j++){
      vCN2[j]->updateTree(pTC.get(), pTC->getNodeIndex(vCN2[j]));
    }
    
    vTree_.push_back(pTC);
  }

}

ComputingTree& ComputingTree::operator=(const ComputingTree& tree)
{
  AbstractParametrizable::operator=(tree);
  parTree_=tree.parTree_;
  pDist_=tree.pDist_;
  isReadyToCompute_=false;
  
  for (size_t i=0; i<tree.vTree_.size(); i++){
    
    std::shared_ptr<CompTree> pTC(tree.vTree_[i]->clone());

    std::vector<std::shared_ptr<ComputingNode> > vCN2=pTC->getAllNodes();
    for (size_t j=0; j<vCN2.size(); j++){
      vCN2[j]->updateTree(pTC.get(), pTC->getNodeIndex(vCN2[j]));
    }
    
    vTree_.push_back(pTC);
  }
  
  return *this;
}

ComputingTree::~ComputingTree()
{
  vTree_.clear();
}

void ComputingTree::checkModelOnEachNode()
{
  for (size_t i=0; i<vTree_.size(); i++)
  {
    vector<shared_ptr<ComputingNode> > vn=vTree_[i]->getAllNodes();
    
    for (size_t j=0; j< vn.size(); j++)
    {
      if (vn[j]==vTree_[i]->getRoot())
        continue;

      if (dynamic_cast<SpeciationComputingNode*>(vn[j].get()))
      {
        if (dynamic_cast<SpeciationComputingNode&>(*vn[j]).getModel()==NULL)
        {
          isReadyToCompute_=false;
          return;
        }
      }
    }
  }
  
  isReadyToCompute_=true;
}
  
void ComputingTree::addModel(const TransitionModel* pSubMod, std::vector<unsigned int>  vBr)
{
  for (size_t i=0; i< getNumberOfClasses(); i++)
    for (size_t j=0; j<vBr.size(); j++)
    {
      shared_ptr<ComputingNode> cn=vTree_[i]->getNode(vBr[j]);
      if (dynamic_cast<SpeciationComputingNode*>(cn.get()))
        dynamic_cast<SpeciationComputingNode&>(*cn).setTransitionModel(pSubMod);
    }
  
  checkModelOnEachNode();
}

void ComputingTree::addModel(const TransitionModel* pSubMod)
{
  if (pSubMod==0)
    return;
  
  vector<unsigned int> vId=vTree_[0]->getAllNodesIndexes();
  
  for (size_t i=0; i< getNumberOfClasses(); i++)
    for (size_t j=0; j< vId.size(); j++)
    {
      shared_ptr<ComputingNode> cn=vTree_[i]->getNode(vId[j]);
      if (dynamic_cast<SpeciationComputingNode*>(cn.get()))
        dynamic_cast<SpeciationComputingNode&>(*cn).setTransitionModel(pSubMod);
    }
  
  isReadyToCompute_=true;
}

void ComputingTree::fireParameterChanged(const ParameterList& pl)
{
  if (!isReadyToCompute_)
    throw Exception("ComputingTree::fireParameterChanged : some nodes do not have a model.");

  bool chDist=false;

  for (size_t i=0; i<pl.size(); i++)
  {
    if (hasParameter(pl[i].getName()))
    {
      string n=pl[i].getName();
      if (n.substr(0,5)!="BrLen"){
        chDist=true;
        continue;
      }
      size_t pt=n.find("_");
      unsigned int nBr=(unsigned int)atoi(n.substr(5,pt-5).c_str());
      for (size_t i2=0; i2<vTree_.size(); i2++)
        vTree_[i2]->getNode(nBr)->setDistanceToFather(pl[i].getValue());
    }
  }

  if (chDist)
  {
    size_t nbCl=pDist_?pDist_->getNumberOfCategories():1;
    for (size_t i=0; i<nbCl; i++){
      std::vector<shared_ptr<ComputingNode> > vCN=vTree_[i]->getAllNodes();
      for (size_t j=0; j<vCN.size(); j++)
        dynamic_cast<SpeciationComputingNode&>(*vCN[j]).setParameterValue("scale",pDist_?pDist_->getCategory(i):1);
    }
  }
}

void ComputingTree::update(vector<unsigned int>& vId, bool flag)
{
  if (!isReadyToCompute_)
    throw Exception("ComputingTree::update : some nodes do not have a model.");

  for (size_t i=0; i<vTree_.size(); i++)
      for (size_t j=0; j<vId.size(); j++)
        vTree_[i]->getNode(vId[j])->update(flag);
}

void ComputingTree::update(unsigned int id, bool flag)
{
  if (!isReadyToCompute_)
    throw Exception("ComputingTree::update :  node " + TextTools::toString(id) + " is not ready.");

  for (size_t i=0; i<vTree_.size(); i++)
    vTree_[i]->getNode(id)->update(flag);
}


Vuint ComputingTree::toBeUpdatedNodes() const
{
  Vuint lId;
  vTree_[0]->getRoot()->toBeUpdatedBelow(lId);

  return lId;
}

void ComputingTree::updateAll()
{
  if (!isReadyToCompute_)
    throw Exception("ComputingTree::updateAll : some nodes are not ready.");

  for (size_t i=0; i<vTree_.size(); i++)
    vTree_[i]->getRoot()->updateAllBelow();
}

