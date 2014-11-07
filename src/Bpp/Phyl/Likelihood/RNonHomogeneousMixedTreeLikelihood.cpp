//
// File: RNonHomogeneousMixedTreeLikelihood.cpp
// Created by: Laurent Gueguen
// Created on: jeudi 11 novembre 2010, à 07h 56
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

#include "RNonHomogeneousMixedTreeLikelihood.h"
#include "../PatternTools.h"
#include "../Model/MixedSubstitutionModel.h"
#include "../Tree/TreeTools.h"

#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

RNonHomogeneousMixedTreeLikelihood::RNonHomogeneousMixedTreeLikelihood(const Tree& tree,
                                                                       MixedSubstitutionModelSet* modelSet,
                                                                       DiscreteDistribution* rDist,
                                                                       bool verbose,
                                                                       bool usePatterns)
throw (Exception) :
  RNonHomogeneousTreeLikelihood(tree, modelSet, rDist, verbose, usePatterns),
  mvTreeLikelihoods_(),
  hyperNode_(modelSet),
  upperNode_(tree.getRootId()),
  main_(true)
{
  if (!modelSet->isFullySetUpFor(tree))
    throw Exception("RNonHomogeneousMixedTreeLikelihood(constructor). Model set is not fully specified.");

  for (size_t i=0;i<modelSet->getNumberOfHyperNodes();i++){
    mvTreeLikelihoods_[tree.getRootId()].push_back(new RNonHomogeneousMixedTreeLikelihood(tree, modelSet, modelSet->getHyperNode(i), upperNode_, rDist, false, usePatterns));
  }

}

/******************************************************************************/

RNonHomogeneousMixedTreeLikelihood::RNonHomogeneousMixedTreeLikelihood(const Tree& tree,
                                                                       const SiteContainer& data,
                                                                       MixedSubstitutionModelSet* modelSet,
                                                                       DiscreteDistribution* rDist,
                                                                       bool verbose,
                                                                       bool usePatterns)
throw (Exception) :
  RNonHomogeneousTreeLikelihood(tree, data, modelSet, rDist, verbose, usePatterns),
  mvTreeLikelihoods_(),
  hyperNode_(modelSet),
  upperNode_(tree.getRootId()),
  main_(true)
{
  if (!modelSet->isFullySetUpFor(tree))
    throw Exception("RNonHomogeneousMixedTreeLikelihood(constructor). Model set is not fully specified.");

  for (size_t i=0;i<modelSet->getNumberOfHyperNodes();i++)
    mvTreeLikelihoods_[tree.getRootId()].push_back(new RNonHomogeneousMixedTreeLikelihood(tree, data, modelSet, modelSet->getHyperNode(i), upperNode_, rDist, false, usePatterns));
}

/******************************************************************************/

RNonHomogeneousMixedTreeLikelihood::RNonHomogeneousMixedTreeLikelihood(const Tree& tree,
                                                                       MixedSubstitutionModelSet* modelSet,
                                                                       const MixedSubstitutionModelSet::HyperNode& hyperNode,
                                                                       int upperNode,
                                                                       DiscreteDistribution* rDist,
                                                                       bool verbose,
                                                                       bool usePatterns) :
  RNonHomogeneousTreeLikelihood(tree, modelSet, rDist, verbose, usePatterns),
  mvTreeLikelihoods_(),
  hyperNode_(hyperNode),
  upperNode_(upperNode),
  main_(false)
{
  if (!modelSet->isFullySetUpFor(tree))
    throw Exception("RNonHomogeneousMixedTreeLikelihood(constructor). Model set is not fully specified.");

  init(usePatterns);
}

/******************************************************************************/

RNonHomogeneousMixedTreeLikelihood::RNonHomogeneousMixedTreeLikelihood(const Tree& tree,
                                                                       const SiteContainer& data,
                                                                       MixedSubstitutionModelSet* modelSet,
                                                                       const MixedSubstitutionModelSet::HyperNode& hyperNode,
                                                                       int upperNode,
                                                                       DiscreteDistribution* rDist,
                                                                       bool verbose,
                                                                       bool usePatterns) :
  RNonHomogeneousTreeLikelihood(tree, data, modelSet, rDist, verbose, usePatterns),
  mvTreeLikelihoods_(),
  hyperNode_(hyperNode),
  upperNode_(upperNode),
  main_(false)
{
  if (!modelSet->isFullySetUpFor(tree))
    throw Exception("RNonHomogeneousMixedTreeLikelihood(constructor). Model set is not fully specified.");

  init(usePatterns);
}

/******************************************************************************/

void RNonHomogeneousMixedTreeLikelihood::init(bool usePatterns)
{
  std::vector<int> vDesc; // vector of the explorated descendents
  int desc;
  vector<int> vn;
  size_t nbmodels = modelSet_->getNumberOfModels();

  const SiteContainer* pdata = getData();
  
  const Tree& tree = getTree();
  
  vDesc.push_back(upperNode_); // start of the exploration

  while (vDesc.size() != 0)  {
    desc = vDesc.back();
    vDesc.pop_back();
    vector<int> vExpMod; // vector of the ids of the MixedModels which
                         // nodes are not in only one subtree under desc

    vector<int> vson = tree.getSonsId(desc);
    std::map<int, vector<int> > mdesc; // map of the subtree nodes for
                                       // each son of desc
    for (size_t i = 0; i < vson.size(); i++)
    {
      std::vector<int> vi;
      TreeTools::getNodesId(tree, vson[i], vi);
      mdesc[vson[i]] = vi;
    }

    for (size_t i = 0; i < nbmodels; i++)
    {
      const MixedSubstitutionModelSet::HyperNode::Node& node = hyperNode_.getNode(i);
      
      if (node.size()>1)
      {
        vn = modelSet_->getNodesWithModel(i); // tree nodes associated to model

        /* Check if the vn members are in the same subtree */
        size_t flag = 0; // count of the subtrees that have vn members
        std::map<int, vector<int> >::iterator it;
        for (it = mdesc.begin(); it != mdesc.end(); it++)
        {
          for (size_t j = 0; j < it->second.size(); j++)
            {
            if (it->second[j] != it->first)
            {
              if (find(vn.begin(), vn.end(), it->second[j]) != vn.end())
              {
                flag += (find(vn.begin(), vn.end(), it->first) != vn.end()) ? 2 : 1; // check if the son
                // has this model too
                break;
              }
            }
            else if (find(vn.begin(), vn.end(), it->first) != vn.end())
              flag++;
          }
          if (flag >= 2)
            break;
        }
        if (flag >= 2)
          vExpMod.push_back(static_cast<int>(i));  // mixed model that must be expanded
      }
    }

    if (vExpMod.size() != 0)
    {
      std::map<int, int> mapmodels;
      size_t ttmodels = 1;
      for (vector<int>::iterator it = vExpMod.begin(); it != vExpMod.end(); it++)
      {
        mapmodels[*it] = static_cast<int>(hyperNode_.getNode(static_cast<size_t>(*it)).size());
        ttmodels *= static_cast<size_t>(mapmodels[*it]);
      }

      for (size_t i = 0; i < ttmodels; i++)
      {
        int s = static_cast<int>(i);
        MixedSubstitutionModelSet::HyperNode hn(hyperNode_);
        
        for (size_t j = 0; j < nbmodels; j++)
        {
          if ((hyperNode_.getNode(j).size() >= 1) && find(vExpMod.begin(), vExpMod.end(), static_cast<int>(j)) != vExpMod.end())
          {
            hn.setModel(j, Vint(1, hyperNode_.getNode(j)[static_cast<size_t>(s % mapmodels[static_cast<int>(j)])]));
            s /= mapmodels[static_cast<int>(j)];
          }
        }
        hn.setProbability((dynamic_cast<MixedSubstitutionModelSet*>(modelSet_))->getHyperNodeProbability(hn));
        RNonHomogeneousMixedTreeLikelihood* pr;

        if (pdata)
          pr = new RNonHomogeneousMixedTreeLikelihood(tree, *pdata, dynamic_cast<MixedSubstitutionModelSet*>(modelSet_), hn, desc, rateDistribution_, false, usePatterns);
        else
          pr = new RNonHomogeneousMixedTreeLikelihood(tree, dynamic_cast<MixedSubstitutionModelSet*>(modelSet_), hn, desc, rateDistribution_, false, usePatterns);
        pr->resetParameters_();
        mvTreeLikelihoods_[desc].push_back(pr);
      }
    }
    else
      for (size_t i = 0; i < vson.size(); i++)
      {
        vDesc.push_back(vson[i]);
      }
  }
}


/******************************************************************************/

RNonHomogeneousMixedTreeLikelihood::RNonHomogeneousMixedTreeLikelihood(
  const RNonHomogeneousMixedTreeLikelihood& lik) :
  RNonHomogeneousTreeLikelihood(lik),
  mvTreeLikelihoods_(),
  hyperNode_(lik.hyperNode_),
  upperNode_(lik.upperNode_),
  main_(lik.main_)
{
  map<int, vector<RNonHomogeneousMixedTreeLikelihood*> >::const_iterator it;
  for (it = lik.mvTreeLikelihoods_.begin(); it != lik.mvTreeLikelihoods_.end(); it++)
  {
    for (size_t i = 0; i < it->second.size(); i++)
    {
      mvTreeLikelihoods_[it->first].push_back(new RNonHomogeneousMixedTreeLikelihood(*it->second[i]));
    }
  }
}

/******************************************************************************/

RNonHomogeneousMixedTreeLikelihood& RNonHomogeneousMixedTreeLikelihood::operator=(
  const RNonHomogeneousMixedTreeLikelihood& lik)
{
  RNonHomogeneousTreeLikelihood::operator=(lik);

  mvTreeLikelihoods_.clear();

  upperNode_ = lik.upperNode_;
  main_ = lik.main_;

  map<int, vector<RNonHomogeneousMixedTreeLikelihood*> >::const_iterator it;
  for (it = lik.mvTreeLikelihoods_.begin(); it != lik.mvTreeLikelihoods_.end(); it++)
  {
    for (size_t i = 0; i < it->second.size(); i++)
    {
      mvTreeLikelihoods_[it->first].push_back(new RNonHomogeneousMixedTreeLikelihood(*it->second[i]));
    }
  }

  hyperNode_=lik.hyperNode_;

  return *this;
}

/******************************************************************************/

RNonHomogeneousMixedTreeLikelihood::~RNonHomogeneousMixedTreeLikelihood()
{
  map<int, vector<RNonHomogeneousMixedTreeLikelihood*> >::iterator it;
  for (it = mvTreeLikelihoods_.begin(); it != mvTreeLikelihoods_.end(); it++)
  {
    for (size_t i = 0; i < it->second.size(); i++)
    {
      delete it->second[i];
    }
  }
}

/******************************************************************************/
 void RNonHomogeneousMixedTreeLikelihood::initialize() throw (Exception)
{
  if (main_)
    initParameters();
  else {
    initBranchLengthsParameters(false);
    includeParameters_(brLenParameters_);
  }

  map<int, vector<RNonHomogeneousMixedTreeLikelihood*> >::iterator it;
  for (it = mvTreeLikelihoods_.begin(); it != mvTreeLikelihoods_.end(); it++)
  {
    for (size_t i = 0; i < it->second.size(); i++)
    {
      it->second[i]->initialize();
    }
  }

  RNonHomogeneousTreeLikelihood::initialize();
}

/******************************************************************************/

 void RNonHomogeneousMixedTreeLikelihood::fireParameterChanged(const ParameterList& params)
{
  if (main_)
    applyParameters();
  else {
    for (size_t i = 0; i < nbNodes_; i++)
      {
        int id = nodes_[i]->getId();
        if (reparametrizeRoot_ && id == root1_)
          {
            const Parameter* rootBrLen = &getParameter("BrLenRoot");
            const Parameter* rootPos = &getParameter("RootPosition");
            nodes_[i]->setDistanceToFather(rootBrLen->getValue() * rootPos->getValue());
          }
        else if (reparametrizeRoot_ && id == root2_)
          {
            const Parameter* rootBrLen = &getParameter("BrLenRoot");
            const Parameter* rootPos = &getParameter("RootPosition");
            nodes_[i]->setDistanceToFather(rootBrLen->getValue() * (1. - rootPos->getValue()));
          }
        else
          {
            const Parameter* brLen = &getParameter(string("BrLen") + TextTools::toString(i));
            if (brLen) nodes_[i]->setDistanceToFather(brLen->getValue());
          }
      }
  }

  map<int, vector<RNonHomogeneousMixedTreeLikelihood*> >::const_iterator it2;
  for (it2 = mvTreeLikelihoods_.begin(); it2 != mvTreeLikelihoods_.end(); it2++)
    for (size_t i = 0; i < it2->second.size(); i++){
      (it2->second)[i]->setProbability((dynamic_cast<MixedSubstitutionModelSet*>(modelSet_))->getHyperNodeProbability((it2->second)[i]->getHyperNode()));
    }

  if (main_){
    for (size_t i=0;i< mvTreeLikelihoods_[upperNode_].size(); i++)
      mvTreeLikelihoods_[upperNode_][i]->matchParametersValues(params);
    rootFreqs_ = modelSet_->getRootFrequencies();
  }
  else {
    if (params.getCommonParametersWith(rateDistribution_->getIndependentParameters()).size() > 0)
      {
        computeAllTransitionProbabilities();
      }
    else
      {
        vector<int> ids;
        vector<string> tmp = params.getCommonParametersWith(modelSet_->getNodeParameters()).getParameterNames();
        for (size_t i = 0; i < tmp.size(); i++)
          {
            vector<int> tmpv = modelSet_->getNodesWithParameter(tmp[i]);
            ids = VectorTools::vectorUnion(ids, tmpv);
          }
        tmp = params.getCommonParametersWith(brLenParameters_).getParameterNames();
        vector<const Node*> nodes;
        for (size_t i = 0; i < ids.size(); i++)
          {
            nodes.push_back(idToNode_[ids[i]]);
          }
        vector<const Node*> tmpv;
        bool test = false;
        for (size_t i = 0; i < tmp.size(); i++)
          {
            if (tmp[i] == "BrLenRoot" || tmp[i] == "RootPosition")
              {
                if (!test)
                  {
                    tmpv.push_back(tree_->getRootNode()->getSon(0));
                    tmpv.push_back(tree_->getRootNode()->getSon(1));
                    test = true; // Add only once.
                  }
              }
            else
              tmpv.push_back(nodes_[TextTools::to < size_t > (tmp[i].substr(5))]);
          }
        nodes = VectorTools::vectorUnion(nodes, tmpv);
        
        for (size_t i = 0; i < nodes.size(); i++){
          computeTransitionProbabilitiesForNode(nodes[i]);
        }
      }

    map<int, vector<RNonHomogeneousMixedTreeLikelihood*> >::iterator it;
    for (it = mvTreeLikelihoods_.begin(); it != mvTreeLikelihoods_.end(); it++)
      {
        for (size_t i = 0; i < it->second.size(); i++)
          {
            it->second[i]->matchParametersValues(params);
          }
      }
  }
  
  if (main_)
    {
      computeTreeLikelihood();
      minusLogLik_ = -getLogLikelihood();
    }
}

/******************************************************************************/
void RNonHomogeneousMixedTreeLikelihood::setData(const SiteContainer& sites) throw (Exception)
{
  RNonHomogeneousTreeLikelihood::setData(sites);
  map<int, vector<RNonHomogeneousMixedTreeLikelihood*> >::iterator it;
  for (it = mvTreeLikelihoods_.begin(); it != mvTreeLikelihoods_.end(); it++)
  {
    for (size_t i = 0; i < it->second.size(); i++)
    {
      it->second[i]->setData(sites);
    }
  }
}


/******************************************************************************/
double RNonHomogeneousMixedTreeLikelihood::getProbability() const
{
  return hyperNode_.getProbability();
}

/******************************************************************************/
void RNonHomogeneousMixedTreeLikelihood::setProbability(double x)
{
  return hyperNode_.setProbability(x);
}

/******************************************************************************/

void RNonHomogeneousMixedTreeLikelihood::computeSubtreeLikelihood(const Node* node)
{
  // if the subtree is divided in several RNonHomogeneousMixedTreeLikelihood*
  if (node->isLeaf())
    return;

  int nodeId=main_?upperNode_:node->getId();
  if (mvTreeLikelihoods_.find(nodeId) != mvTreeLikelihoods_.end()) {

    size_t nbSites  = likelihoodData_->getLikelihoodArray(nodeId).size();
    
    // Must reset the likelihood array first (i.e. set all of them to 0):
    VVVdouble* _likelihoods_node = &likelihoodData_->getLikelihoodArray(nodeId);
    for (size_t i = 0; i < nbSites; i++)
      {
        // For each site in the sequence,
        VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
        for (size_t c = 0; c < nbClasses_; c++)
          {
              // For each rate classe,
            Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
            for (size_t x = 0; x < nbStates_; x++)
              {
                // For each initial state,
                (*_likelihoods_node_i_c)[x] = 0.;
              }
          }
      }

    if (getProbability()==0)
      return;
  
    vector<RNonHomogeneousMixedTreeLikelihood* > vr = mvTreeLikelihoods_[nodeId];
    for (size_t t = 0; t < vr.size(); t++)
      vr[t]->computeSubtreeLikelihood(node);

    // for each specific subtree
    for (size_t t = 0; t < vr.size(); t++)
    {
      VVVdouble* _vt_likelihoods_node = &vr[t]->likelihoodData_->getLikelihoodArray(nodeId);
      for (size_t i = 0; i < nbSites; i++)
      {
        // For each site in the sequence,
        VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
        for (size_t c = 0; c < nbClasses_; c++)
        {
          // For each rate classe,
          Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
          Vdouble* _vt_likelihoods_node_i_c = &(*_vt_likelihoods_node)[i][c];
          for (size_t x = 0; x < nbStates_; x++)
          {
            (*_likelihoods_node_i_c)[x] +=  (*_vt_likelihoods_node_i_c)[x] * vr[t]->getProbability()/getProbability();
          }
        }
      }
    }
  }
  

  // otherwise...

  // nb: if the subtree is made of independent branches the computing is
  // as in the non mixed case, where the mean of the probas of
  // transition of a mixed model are taken.

  else
    RNonHomogeneousTreeLikelihood::computeSubtreeLikelihood(node); 
}

/******************************************************************************
*                           First Order Derivatives                          *
******************************************************************************/
void RNonHomogeneousMixedTreeLikelihood::computeTreeDLikelihood(const string& variable)
{
  const Node* father, father2;
    
    
  if (main_)
    father = tree_->getRootNode();
  else {
    if ((variable == "BrLenRoot") ||  (variable == "RootPosition"))
      father = tree_->getRootNode();
    else
      {
        size_t brI = TextTools::to<size_t>(variable.substr(5));
        const Node* branch = nodes_[brI];
        father = branch->getFather();
      }
  }
  
  bool flok = 0;
  while (father){
    if (mvTreeLikelihoods_.find(father->getId()) != mvTreeLikelihoods_.end()) {
      flok = 1;
      break;
    }
    if (father->getId() == upperNode_)
      break;
    father = father->getFather();
  }
    
  if (flok) {  // there is an expanded model above the derivated branch
    int fatherId = father->getId();
    // Compute dLikelihoods array for the father node.
    // Fist initialize to 0:
    VVVdouble* _dLikelihoods_father = &likelihoodData_->getDLikelihoodArray(fatherId);
    size_t nbSites  = _dLikelihoods_father->size();
    for (size_t i = 0; i < nbSites; i++) {
      VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
      for (size_t c = 0; c < nbClasses_; c++) {
        Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
        for (size_t s = 0; s < nbStates_; s++){
          (*_dLikelihoods_father_i_c)[s] = 0.;
        }
      }
    }

    if (getProbability()!=0){
      vector<RNonHomogeneousMixedTreeLikelihood* > vr = mvTreeLikelihoods_[fatherId];
      for (size_t t = 0; t < vr.size(); t++)
        vr[t]->computeTreeDLikelihood(variable);
      
    
      // for each specific subtree
      for (size_t t = 0; t < vr.size(); t++) {
        VVVdouble* _vt_dLikelihoods_father = &vr[t]->likelihoodData_->getDLikelihoodArray(fatherId);
        for (size_t i = 0; i < nbSites; i++){
          // For each site in the sequence,
          VVdouble* _dLikelihoods_father_i = &(*_dLikelihoods_father)[i];
          for (size_t c = 0; c < nbClasses_; c++){
            // For each rate classe,
            Vdouble* _dLikelihoods_father_i_c = &(*_dLikelihoods_father_i)[c];
            Vdouble* _vt_dLikelihoods_father_i_c = &(*_vt_dLikelihoods_father)[i][c];
            for (size_t x = 0; x < nbStates_; x++) {
              (*_dLikelihoods_father_i_c)[x] +=  (*_vt_dLikelihoods_father_i_c)[x] * vr[t]->getProbability()/getProbability();
            }
          }
        }
      }
    }
    computeDownSubtreeDLikelihood(father);
  }
  else
    RNonHomogeneousTreeLikelihood::computeTreeDLikelihood(variable); 
}

void RNonHomogeneousMixedTreeLikelihood::computeDownSubtreeDLikelihood(const Node* node)
{
  const Node* father = node->getFather();
  // // We assume that the _dLikelihoods array has been filled for the current node 'node'.
  // // We will evaluate the array for the father node.
  if (father == 0)
    return;  // We reached the up!

  if (node->getId() == upperNode_)
    return;  // We reached the top of the subtree

  RNonHomogeneousTreeLikelihood::computeDownSubtreeDLikelihood(node);
}

/******************************************************************************
*                           Second Order Derivatives                         *
******************************************************************************/
void RNonHomogeneousMixedTreeLikelihood::computeTreeD2Likelihood(const string& variable)
{
  const Node* father, father2;

  if (main_)
    father = tree_->getRootNode();
  else {
    if ((variable == "BrLenRoot") ||  (variable == "RootPosition"))
      father = tree_->getRootNode();
    else
      {
        size_t brI = TextTools::to<size_t>(variable.substr(5));
        const Node* branch = nodes_[brI];
        father = branch->getFather();
      }
  }
  
  bool flok = 0;
  while (father){
    if (mvTreeLikelihoods_.find(father->getId()) != mvTreeLikelihoods_.end())
      {
        flok = 1;
        break;
      }
    if (father->getId() == upperNode_)
      break;
    father = father->getFather();
  }

  if (flok)  // there is an expanded model above the derivated branch
    {
      int fatherId = father->getId();
      // Compute d2Likelihoods array for the father node.
      // Fist initialize to 0:
      VVVdouble* _d2Likelihoods_father = &likelihoodData_->getD2LikelihoodArray(fatherId);
      size_t nbSites  = _d2Likelihoods_father->size();
      for (size_t i = 0; i < nbSites; i++) {
        VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
        for (size_t c = 0; c < nbClasses_; c++) {
          Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
          for (size_t s = 0; s < nbStates_; s++) {
            (*_d2Likelihoods_father_i_c)[s] = 0.;
          }
        }
      }

      if (getProbability()!=0){
        
        vector<RNonHomogeneousMixedTreeLikelihood* > vr = mvTreeLikelihoods_[fatherId];
        for (size_t t = 0; t < vr.size(); t++)
          vr[t]->computeTreeD2Likelihood(variable);
      
        // for each specific subtree
        for (size_t t = 0; t < vr.size(); t++) {
          VVVdouble* _vt_d2Likelihoods_father = &vr[t]->likelihoodData_->getD2LikelihoodArray(fatherId);
          for (size_t i = 0; i < nbSites; i++) {
            // For each site in the sequence,
            VVdouble* _d2Likelihoods_father_i = &(*_d2Likelihoods_father)[i];
            for (size_t c = 0; c < nbClasses_; c++){
              // For each rate classe,
              Vdouble* _d2Likelihoods_father_i_c = &(*_d2Likelihoods_father_i)[c];
              Vdouble* _vt_d2Likelihoods_father_i_c = &(*_vt_d2Likelihoods_father)[i][c];
              for (size_t x = 0; x < nbStates_; x++) {
                (*_d2Likelihoods_father_i_c)[x] +=  (*_vt_d2Likelihoods_father_i_c)[x] * vr[t]->getProbability() / getProbability();
              }
            }
          }
        }
      }
      computeDownSubtreeD2Likelihood(father);
    }
  else
    RNonHomogeneousTreeLikelihood::computeTreeD2Likelihood(variable);
}


void RNonHomogeneousMixedTreeLikelihood::computeDownSubtreeD2Likelihood(const Node* node)
{
  const Node* father = node->getFather();
  // // We assume that the _dLikelihoods array has been filled for the current node 'node'.
  // // We will evaluate the array for the father node.
  if (father == 0)
    return;  // We reached the up!

  if (node->getId() == upperNode_)
    return;  // We reached the top of the subtree

  RNonHomogeneousTreeLikelihood::computeDownSubtreeD2Likelihood(node);
}


/*******************************************************************************/

void RNonHomogeneousMixedTreeLikelihood::computeTransitionProbabilitiesForNode(const Node* node)
{
  const SubstitutionModel* model = modelSet_->getModelForNode(node->getId());
  size_t modelnum = modelSet_->getModelIndexForNode(node->getId());

  vector<const SubstitutionModel*> vModel;
  vector<double> vProba;
  
  const MixedSubstitutionModelSet::HyperNode::Node& nd = hyperNode_.getNode(modelnum);
  if (nd.size() == 0) {
    vModel.push_back(model);
    vProba.push_back(1);
  }
  else {
    const MixedSubstitutionModel* mmodel = dynamic_cast<const MixedSubstitutionModel*>(model);
    double x = 0;
    for (size_t i = 0; i < nd.size(); ++i){
      vModel.push_back(mmodel->getNModel(static_cast<size_t>(nd[i])));
      vProba.push_back(mmodel->getNProbability(static_cast<size_t>(nd[i])));
      x += vProba[i];
    }
    if (x != 0)
      for (size_t i = 0; i < nd.size(); ++i)
        vProba[i] /= x;
  }

  double l = node->getDistanceToFather();
  // Computes all pxy and pyx once for all:
  VVVdouble* pxy__node = &pxy_[node->getId()];
  for (size_t c = 0; c < nbClasses_; c++) {
    VVdouble* pxy__node_c = &(*pxy__node)[c];
    for (size_t x = 0; x < nbStates_; x++){
      Vdouble* pxy__node_c_x = &(*pxy__node_c)[x];
      for (size_t y = 0; y < nbStates_; y++){
        (*pxy__node_c_x)[y] = 0;
      }
    }
    
    for (size_t i=0;i<vModel.size();i++){
      RowMatrix<double> Q = vModel[i]->getPij_t(l * rateDistribution_->getCategory(c));
      for (size_t x = 0; x < nbStates_; x++){
        Vdouble* pxy__node_c_x = &(*pxy__node_c)[x];
        for (size_t y = 0; y < nbStates_; y++){
          (*pxy__node_c_x)[y] += vProba[i] * Q(x, y);
        }
      }
    }
  }
  
  if (computeFirstOrderDerivatives_) {
    // Computes all dpxy/dt once for all:
    VVVdouble* dpxy__node = &dpxy_[node->getId()];

    for (size_t c = 0; c < nbClasses_; c++){
      VVdouble* dpxy__node_c = &(*dpxy__node)[c];
      double rc = rateDistribution_->getCategory(c);
      for (size_t x = 0; x < nbStates_; x++){
        Vdouble* dpxy__node_c_x = &(*dpxy__node_c)[x];
        for (size_t y = 0; y < nbStates_; y++){
          (*dpxy__node_c_x)[y] = 0;
        }
      }

      for (size_t i=0;i<vModel.size();i++){
        RowMatrix<double> dQ = vModel[i]->getdPij_dt(l * rc);

        for (size_t x = 0; x < nbStates_; x++){
          Vdouble* dpxy__node_c_x = &(*dpxy__node_c)[x];
          for (size_t y = 0; y < nbStates_; y++){
            (*dpxy__node_c_x)[y] +=  vProba[i] * rc * dQ(x, y);
          }
        }
      }
    }
  }
  
  if (computeSecondOrderDerivatives_) {
    // Computes all d2pxy/dt2 once for all:
    VVVdouble* d2pxy__node = &d2pxy_[node->getId()];
    for (size_t c = 0; c < nbClasses_; c++){
      VVdouble* d2pxy__node_c = &(*d2pxy__node)[c];
      for (size_t x = 0; x < nbStates_; x++){
        Vdouble* d2pxy__node_c_x = &(*d2pxy__node_c)[x];
        for (size_t y = 0; y < nbStates_; y++){
          (*d2pxy__node_c_x)[y] = 0;
        }
      }
      
      double rc =  rateDistribution_->getCategory(c);
      for (size_t i=0;i<vModel.size();i++){
        RowMatrix<double> d2Q = vModel[i]->getd2Pij_dt2(l * rc);
        for (size_t x = 0; x < nbStates_; x++){
          Vdouble* d2pxy__node_c_x = &(*d2pxy__node_c)[x];
          for (size_t y = 0; y < nbStates_; y++){
            (*d2pxy__node_c_x)[y] +=  vProba[i] * rc * rc * d2Q(x, y);
          }
        }
      }
    }
  }
  
}

/*******************************************************************************/

