//
// File: DRHomogeneousTreeLikelihood.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 17 18:14:51 2003
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#include "DRHomogeneousTreeLikelihood.h"
#include "PatternTools.h"

//From SeqLib:
#include <Seq/SiteTools.h>
#include <Seq/AlignedSequenceContainer.h>
#include <Seq/SequenceContainerTools.h>

// From Utils:
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

DRHomogeneousTreeLikelihood::DRHomogeneousTreeLikelihood(
  const Tree & tree,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose)
throw (Exception):
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  _likelihoodData(NULL)
{
  _init();
}

/******************************************************************************/

DRHomogeneousTreeLikelihood::DRHomogeneousTreeLikelihood(
  const Tree & tree,
  const SiteContainer & data,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose)
throw (Exception):
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  _likelihoodData(NULL)
{
  _init(); 
  setData(data);
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::_init() throw (Exception)
{
  _likelihoodData = new DRASDRTreeLikelihoodData(*_tree, _rateDistribution->getNumberOfCategories());
}

/******************************************************************************/

DRHomogeneousTreeLikelihood::DRHomogeneousTreeLikelihood(const DRHomogeneousTreeLikelihood & lik):
  AbstractHomogeneousTreeLikelihood(lik),
  _likelihoodData(NULL)
{
  _likelihoodData = dynamic_cast<DRASDRTreeLikelihoodData *>(lik._likelihoodData->clone());
  _likelihoodData->setTree(*_tree);
}

/******************************************************************************/

DRHomogeneousTreeLikelihood & DRHomogeneousTreeLikelihood::operator=(const DRHomogeneousTreeLikelihood & lik)
{
  AbstractHomogeneousTreeLikelihood::operator=(lik);
  if(_likelihoodData) delete _likelihoodData;
  _likelihoodData = dynamic_cast<DRASDRTreeLikelihoodData *>(lik._likelihoodData->clone());
  _likelihoodData->setTree(*_tree);
  return *this;
}

/******************************************************************************/

DRHomogeneousTreeLikelihood::~DRHomogeneousTreeLikelihood()
{
  delete _likelihoodData;
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::setData(const SiteContainer & sites) throw (Exception)
{
  if(_data) delete _data;
  _data = PatternTools::getSequenceSubset(sites, *_tree->getRootNode());
   if(verbose_) ApplicationTools::displayTask("Initializing data structure");
  _likelihoodData->initLikelihoods(*_data, *model_);
  if(verbose_) ApplicationTools::displayTaskDone();

  nbSites_ = _likelihoodData->getNumberOfSites();
  nbDistinctSites_ = _likelihoodData->getNumberOfDistinctSites();
  nbStates_ = _likelihoodData->getNumberOfStates();
  
  if(verbose_) ApplicationTools::displayResult("Number of distinct sites",
      TextTools::toString(nbDistinctSites_));
  _initialized = false;
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLikelihood() const
{
  double l = 1.;
  Vdouble * lik = & _likelihoodData->getRootRateSiteLikelihoodArray(); 
  const vector<unsigned int> * w = & _likelihoodData->getWeights();
  for(unsigned int i = 0; i < nbDistinctSites_; i++)
  {
    l *= std::pow((*lik)[i], (int)(* w)[i]);
  }
  return l;
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihood() const
{
  double ll = 0;
  Vdouble * lik = & _likelihoodData->getRootRateSiteLikelihoodArray(); 
  const vector<unsigned int> * w = & _likelihoodData->getWeights();
  vector<double> la(nbDistinctSites_);
  for(unsigned int i = 0; i < nbDistinctSites_; i++)
  {
    la[i] = (* w)[i] * log((* lik)[i]);
  }
  sort(la.begin(), la.end());
  for(unsigned int i = nbDistinctSites_; i > 0; i--)
    ll += la[i-1];
  return ll;
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLikelihoodForASite(unsigned int site) const
{
  return _likelihoodData->getRootRateSiteLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)];
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihoodForASite(unsigned int site) const
{
  return log(_likelihoodData->getRootRateSiteLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)]);
}

/******************************************************************************/
double DRHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
  return _likelihoodData->getRootSiteLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)][rateClass];
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
  return log(_likelihoodData->getRootSiteLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)][rateClass]);
}

/******************************************************************************/  

double DRHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
  return _likelihoodData->getRootLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)][rateClass][state];
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
  return log(_likelihoodData->getRootLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)][rateClass][state]);
}

/******************************************************************************/  

void DRHomogeneousTreeLikelihood::setParameters(const ParameterList & parameters)
  throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::fireParameterChanged(const ParameterList & params)
{
  applyParameters();

  if(_rateDistribution->getParameters().getCommonParametersWith(params).size() > 0
  || model_->getParameters().getCommonParametersWith(params).size() > 0)
  {
    //Rate parameter changed, need to recompute all probs:
    computeAllTransitionProbabilities();
  }
  else if(params.size() > 0)
  {
    //We may save some computations:
    for(unsigned int i = 0; i < params.size(); i++)
    {
      string s = params[i]->getName();
      if(s.substr(0,5) == "BrLen")
      {
        //Branch length parameter:
        computeTransitionProbabilitiesForNode(nodes_[TextTools::to<unsigned int>(s.substr(5))]);
      }
    }
  }

  computeTreeLikelihood();
  if(_computeFirstOrderDerivatives)
  {
    computeTreeDLikelihoods();  
  }
  if(_computeSecondOrderDerivatives)
  {
    computeTreeD2Likelihoods();
  }
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getValue() const
throw (Exception)
{
  if(!isInitialized()) throw Exception("DRHomogeneousTreeLikelihood::getValue(). Instance is not initialized.");
  return - getLogLikelihood();
}

/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/  

void DRHomogeneousTreeLikelihood::computeTreeDLikelihoodAtNode(const Node * node)
{
  const Node * father = node->getFather();
  VVVdouble * _likelihoods_father_node = & _likelihoodData->getLikelihoodArray(father->getId(), node->getId());
  Vdouble * _dLikelihoods_node = & _likelihoodData->getDLikelihoodArray(node->getId());
  VVVdouble *  pxy__node = &  pxy_[node->getId()];
  VVVdouble * dpxy__node = & dpxy_[node->getId()];
  VVVdouble larray;
  computeLikelihoodAtNode(father->getId(), larray);
  Vdouble * rootLikelihoodsSR = & _likelihoodData->getRootRateSiteLikelihoodArray();

  double dLi, dLic, dLicx, numerator, denominator;
  
  for(unsigned int i = 0; i < nbDistinctSites_; i++)
  {
    VVdouble * _likelihoods_father_node_i = & (* _likelihoods_father_node)[i];
    VVdouble * larray_i = & larray[i];
    dLi = 0;
    for(unsigned int c = 0; c < nbClasses_; c++)
    {
      Vdouble * _likelihoods_father_node_i_c = & (* _likelihoods_father_node_i)[c];
      Vdouble * larray_i_c = & (* larray_i)[c];
      VVdouble *  pxy__node_c = & (*  pxy__node)[c];
      VVdouble * dpxy__node_c = & (* dpxy__node)[c];
      dLic = 0;
      for(unsigned int x = 0; x < nbStates_; x++)
      {
        numerator = 0;
        denominator = 0;
        Vdouble *  pxy__node_c_x = & (*  pxy__node_c)[x];
        Vdouble * dpxy__node_c_x = & (* dpxy__node_c)[x];
        dLicx = 0;
        for(unsigned int y = 0; y < nbStates_; y++)
        {
          numerator   += (* dpxy__node_c_x)[y] * (* _likelihoods_father_node_i_c)[y];
          denominator += (*  pxy__node_c_x)[y] * (* _likelihoods_father_node_i_c)[y];
        }
        dLicx = (* larray_i_c)[x] * numerator / denominator;
        dLic += dLicx;  
      }
      dLi += _rateDistribution->getProbability(c) * dLic;
    }
    (* _dLikelihoods_node)[i] = dLi / (* rootLikelihoodsSR)[i]; 
  }
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeTreeDLikelihoods()
{
  for(unsigned int k = 0; k < nbNodes_; k++)
  {
    computeTreeDLikelihoodAtNode(nodes_[k]);
  }
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getFirstOrderDerivative(const string & variable) const
throw (Exception)
{ 
  const Parameter * p = &getParameter(variable);
  if(getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameters are not implemented.");
  }
  if(getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }
  
  //
  // Computation for branch lengths:
  //
  
  // Get the node with the branch whose length must be derivated:
  unsigned int brI = TextTools::to<unsigned int>(variable.substr(5));
  const Node * branch = nodes_[brI];
  Vdouble * _dLikelihoods_branch = & _likelihoodData->getDLikelihoodArray(branch->getId());
  double d = 0;
  const vector<unsigned int> * w = & _likelihoodData->getWeights();
  for(unsigned int i = 0; i < nbDistinctSites_; i++)
    d += (* w)[i] * (* _dLikelihoods_branch)[i];
  return -d;
}

/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/  

void DRHomogeneousTreeLikelihood::computeTreeD2LikelihoodAtNode(const Node * node)
{
  const Node * father = node->getFather();
  VVVdouble * _likelihoods_father_node = & _likelihoodData->getLikelihoodArray(father->getId(), node->getId());
  Vdouble * _d2Likelihoods_node = & _likelihoodData->getD2LikelihoodArray(node->getId());  
  VVVdouble *   pxy__node = &   pxy_[node->getId()];
  VVVdouble * d2pxy__node = & d2pxy_[node->getId()];
  VVVdouble larray;
  computeLikelihoodAtNode(father->getId(), larray);
  Vdouble * rootLikelihoodsSR = & _likelihoodData->getRootRateSiteLikelihoodArray();
  
  double d2Li, d2Lic, d2Licx, numerator, denominator;
  
  for(unsigned int i = 0; i < nbDistinctSites_; i++)
  {
    VVdouble * _likelihoods_father_node_i = & (* _likelihoods_father_node)[i];
    VVdouble * larray_i = & larray[i];
    d2Li = 0;
    for(unsigned int c = 0; c < nbClasses_; c++)
    {
      Vdouble * _likelihoods_father_node_i_c = & (* _likelihoods_father_node_i)[c];
      Vdouble * larray_i_c = & (* larray_i)[c];
      VVdouble *   pxy__node_c = & (*   pxy__node)[c];
      VVdouble * d2pxy__node_c = & (* d2pxy__node)[c];
      d2Lic = 0;
      for(unsigned int x = 0; x < nbStates_; x++)
      {
        numerator = 0;
        denominator = 0;
        Vdouble *   pxy__node_c_x = & (*   pxy__node_c)[x];
        Vdouble * d2pxy__node_c_x = & (* d2pxy__node_c)[x];
        d2Licx = 0;
        for(unsigned int y = 0; y < nbStates_; y++)
        {
          numerator   += (* d2pxy__node_c_x)[y] * (* _likelihoods_father_node_i_c)[y];
          denominator += (*   pxy__node_c_x)[y] * (* _likelihoods_father_node_i_c)[y];
        }
        d2Licx = (* larray_i_c)[x] * numerator / denominator;
        d2Lic += d2Licx;
      }
      d2Li += _rateDistribution->getProbability(c) * d2Lic;
    }
    (* _d2Likelihoods_node)[i] = d2Li / (* rootLikelihoodsSR)[i]; 
  }
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeTreeD2Likelihoods()
{
  for(unsigned int k = 0; k < nbNodes_; k++)
  {
    computeTreeD2LikelihoodAtNode(nodes_[k]);
  }
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getSecondOrderDerivative(const string & variable) const
throw (Exception)
{
  const Parameter * p = &getParameter(variable);
  if(getRateDistributionParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to rate distribution parameters are not implemented.");
  }
  if(getSubstitutionModelParameters().hasParameter(variable))
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }
  
  //
  // Computation for branch lengths:
  //
  
  // Get the node with the branch whose length must be derivated:
  unsigned int brI = TextTools::to<unsigned int>(variable.substr(5));
  const Node * branch = nodes_[brI];
  Vdouble * _dLikelihoods_branch = & _likelihoodData->getDLikelihoodArray(branch->getId());
  Vdouble * _d2Likelihoods_branch = & _likelihoodData->getD2LikelihoodArray(branch->getId());
  double d2 = 0;
  const vector<unsigned int> * w = & _likelihoodData->getWeights();
  for(unsigned int i = 0; i < nbDistinctSites_; i++)
    d2 += (* w)[i] * ((* _d2Likelihoods_branch)[i] - pow((* _dLikelihoods_branch)[i], 2));
  return -d2;
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::resetLikelihoodArrays(const Node * node)
{
  for(unsigned int n = 0; n < node->getNumberOfSons(); n++)
  {
    const Node * subNode = node->getSon(n);
    resetLikelihoodArray(_likelihoodData->getLikelihoodArray(node->getId(), subNode->getId()));
  }
  if(node->hasFather())
  {
    const Node * father = node->getFather();
    resetLikelihoodArray(_likelihoodData->getLikelihoodArray(node->getId(), father->getId()));
  }
}

/******************************************************************************/
  
void DRHomogeneousTreeLikelihood::computeTreeLikelihood()
{
  computeSubtreeLikelihoodPostfix(_tree->getRootNode());
  computeSubtreeLikelihoodPrefix(_tree->getRootNode());
  computeRootLikelihood();
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeSubtreeLikelihoodPostfix(const Node * node)
{
//  if(node->isLeaf()) return;
  //cout << node->getId() << "\t" << (node->hasName()?node->getName():"") << endl;
  if(node->getNumberOfSons() == 0) return;

  // Set all likelihood arrays to 1 for a start:
  resetLikelihoodArrays(node);
  
  map<int, VVVdouble> * _likelihoods_node = & _likelihoodData->getLikelihoodArrays(node->getId());
  unsigned int nbNodes = node->getNumberOfSons();
  for(unsigned int l = 0; l < nbNodes; l++)
  {
    //For each son node...  

    const Node * son = node->getSon(l);
    VVVdouble * _likelihoods_node_son = & (* _likelihoods_node)[son->getId()];
    
    if(son->isLeaf())
    {
      VVdouble * _likelihoods_leaf = & _likelihoodData->getLeafLikelihoods(son->getId());
      for(unsigned int i = 0; i < nbDistinctSites_; i++)
      {
        //For each site in the sequence,
        Vdouble * _likelihoods_leaf_i = & (* _likelihoods_leaf)[i];
        VVdouble * _likelihoods_node_son_i = & (* _likelihoods_node_son)[i];
        for(unsigned int c = 0; c < nbClasses_; c++)
        {
          //For each rate classe,
          Vdouble * _likelihoods_node_son_i_c = & (* _likelihoods_node_son_i)[c];
          for(unsigned int x = 0; x < nbStates_; x++)
          {
            //For each initial state,
            (* _likelihoods_node_son_i_c)[x] = (* _likelihoods_leaf_i)[x];
          }
        }
      }
    }
    else
    {
      computeSubtreeLikelihoodPostfix(son); //Recursive method:
      unsigned int nbSons = son->getNumberOfSons();
      map<int, VVVdouble> * _likelihoods_son = & _likelihoodData->getLikelihoodArrays(son->getId());
      
      vector<const VVVdouble *> iLik(nbSons);
      vector<const VVVdouble *> tProb(nbSons);
      for(unsigned int n = 0; n < nbSons; n++)
      {
        const Node * sonSon = son->getSon(n);
        tProb[n] = & pxy_[sonSon->getId()];
        iLik[n] = & (* _likelihoods_son)[sonSon->getId()];
      }
      computeLikelihoodFromArrays(iLik, tProb, *_likelihoods_node_son, nbSons, nbDistinctSites_, nbClasses_, nbStates_, false); 
    }
  }
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeSubtreeLikelihoodPrefix(const Node * node)
{
  if(! node->hasFather())
  { 
    // 'node' is the root of the tree.  
    // Just call the method on each son node:
    unsigned int nbSons = node->getNumberOfSons();
    for(unsigned int n = 0; n < nbSons; n++)
      computeSubtreeLikelihoodPrefix(node->getSon(n));
    return;
  }
  else
  {
    const Node * father = node->getFather();
    map<int, VVVdouble> * _likelihoods_node = & _likelihoodData->getLikelihoodArrays(node->getId());
    map<int, VVVdouble> * _likelihoods_father = & _likelihoodData->getLikelihoodArrays(father->getId());
    VVVdouble * _likelihoods_node_father = & (* _likelihoods_node)[father->getId()];
    if(node->isLeaf())
    {
      resetLikelihoodArray(*_likelihoods_node_father);
    }
  
    if(father->isLeaf())
    { 
      // If the tree is rooted by a leaf
      VVdouble * _likelihoods_leaf = & _likelihoodData->getLeafLikelihoods(father->getId());
      for(unsigned int i = 0; i < nbDistinctSites_; i++)
      {
        //For each site in the sequence,
        Vdouble * _likelihoods_leaf_i = & (* _likelihoods_leaf)[i];
        VVdouble * _likelihoods_node_father_i = & (* _likelihoods_node_father)[i];
        for(unsigned int c = 0; c < nbClasses_; c++)
        {
          //For each rate classe,
          Vdouble * _likelihoods_node_father_i_c = & (* _likelihoods_node_father_i)[c];
          for(unsigned int x = 0; x < nbStates_; x++)
          {
            //For each initial state,
            (* _likelihoods_node_father_i_c)[x] = (* _likelihoods_leaf_i)[x];
          }
        }
      }
    }
    else 
    {
      vector<const Node *> nodes;
      // Add brothers:
      unsigned int nbFatherSons = father->getNumberOfSons();
      for(unsigned int n = 0; n < nbFatherSons; n++)
      {
        const Node * son = father->getSon(n);
        if(son->getId() != node->getId()) nodes.push_back(son); //This is a real brother, not current node!
      }
      // Now the real stuff... We've got to compute the likelihoods for the
      // subtree defined by node 'father'.
      // This is the same as postfix method, but with different subnodes.
  
      unsigned int nbSons = nodes.size(); // In case of a bifurcating tree, this is equal to 1, excepted for the root.
      
      vector<const VVVdouble *> iLik(nbSons);
      vector<const VVVdouble *> tProb(nbSons);
      for(unsigned int n = 0; n < nbSons; n++)
      {
        const Node * fatherSon = nodes[n];
        tProb[n] = & pxy_[fatherSon->getId()];
        iLik[n] = & (* _likelihoods_father)[fatherSon->getId()];
      }
    
      if(father->hasFather())
      {
        const Node * fatherFather = father->getFather();
        computeLikelihoodFromArrays(iLik, tProb, & (* _likelihoods_father)[fatherFather->getId()], & pxy_[father->getId()], *_likelihoods_node_father, nbSons, nbDistinctSites_, nbClasses_, nbStates_, false); 
      }
      else
      {
        computeLikelihoodFromArrays(iLik, tProb, *_likelihoods_node_father, nbSons, nbDistinctSites_, nbClasses_, nbStates_, false); 
      }
    }

    if(!father->hasFather())
    {
      //We have to account for the root frequencies:
      for(unsigned int i = 0; i < nbDistinctSites_; i++)
      {
        VVdouble * _likelihoods_node_father_i = & (*_likelihoods_node_father)[i];
        for(unsigned int c = 0; c < nbClasses_; c++)
        {
          Vdouble * _likelihoods_node_father_i_c = & (* _likelihoods_node_father_i)[c];
          for(unsigned int x = 0; x < nbStates_; x++)
          {
            (* _likelihoods_node_father_i_c)[x] *= rootFreqs_[x];
          }
        }
      }
    }

    // Call the method on each son node:
    unsigned int nbNodeSons = node->getNumberOfSons();
    for(unsigned int i = 0; i < nbNodeSons; i++)
      computeSubtreeLikelihoodPrefix(node->getSon(i)); //Recursive method.
  }
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeRootLikelihood()
{
  const Node * root = _tree->getRootNode();
  VVVdouble * rootLikelihoods = & _likelihoodData->getRootLikelihoodArray();
  // Set all likelihoods to 1 for a start:
  if(root->isLeaf())
  {
    VVdouble * leavesLikelihoods_root = & _likelihoodData->getLeafLikelihoods(root->getId());
    for(unsigned int i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble * rootLikelihoods_i = & (* rootLikelihoods)[i];
      Vdouble * leavesLikelihoods_root_i = & (* leavesLikelihoods_root)[i];
      for(unsigned int c = 0; c < nbClasses_; c++)
      {
        Vdouble * rootLikelihoods_i_c = & (* rootLikelihoods_i)[c];
        for(unsigned int x = 0; x < nbStates_; x++)
        {
          (* rootLikelihoods_i_c)[x] = (* leavesLikelihoods_root_i)[x];
        }
      }
    }
  }
  else
  {
    resetLikelihoodArray(* rootLikelihoods);
  }
  
  map<int, VVVdouble> * likelihoods_root = & _likelihoodData->getLikelihoodArrays(root->getId());
  unsigned int nbNodes = root->getNumberOfSons();
  vector<const VVVdouble *> iLik(nbNodes);
  vector<const VVVdouble *> tProb(nbNodes);
  for(unsigned int n = 0; n < nbNodes; n++)
  {
    const Node * son = root->getSon(n);
    tProb[n] = & pxy_[son->getId()];
    iLik[n] = & (* likelihoods_root)[son->getId()];
  }
  computeLikelihoodFromArrays(iLik, tProb, *rootLikelihoods, nbNodes, nbDistinctSites_, nbClasses_, nbStates_, false);

  Vdouble p = _rateDistribution->getProbabilities();
  VVdouble * rootLikelihoodsS  = & _likelihoodData->getRootSiteLikelihoodArray();
  Vdouble  * rootLikelihoodsSR = & _likelihoodData->getRootRateSiteLikelihoodArray();
  for(unsigned int i = 0; i < nbDistinctSites_; i++)
  {
    //For each site in the sequence,
    VVdouble * rootLikelihoods_i = & (* rootLikelihoods)[i];
    Vdouble * rootLikelihoodsS_i = & (* rootLikelihoodsS)[i];
    (* rootLikelihoodsSR)[i] = 0;
    for(unsigned int c = 0; c < nbClasses_; c++)
    {
      //For each rate classe,
      Vdouble * rootLikelihoods_i_c = & (* rootLikelihoods_i)[c];
      double * rootLikelihoodsS_i_c = & (* rootLikelihoodsS_i)[c];
      (* rootLikelihoodsS_i_c) = 0;
      for(unsigned int x = 0; x < nbStates_; x++)
      {
        //For each initial state,
        (* rootLikelihoodsS_i_c) += rootFreqs_[x] * (* rootLikelihoods_i_c)[x];
      }
      (* rootLikelihoodsSR)[i] += p[c] * (* rootLikelihoodsS_i_c);
    }

    //Final checking (for numerical errors):
    if((* rootLikelihoodsSR)[i] < 0) (* rootLikelihoodsSR)[i] = 0.;
  }
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeLikelihoodAtNode(int nodeId, VVVdouble& likelihoodArray) const
{
  const Node * node = _tree->getNode(nodeId);

  likelihoodArray.resize(nbDistinctSites_);
  map<int, VVVdouble> * likelihoods_node = & _likelihoodData->getLikelihoodArrays(nodeId);
  
  //Initialize likelihood array:
  if(node->isLeaf())
  {
    VVdouble * leavesLikelihoods_node = & _likelihoodData->getLeafLikelihoods(nodeId);
    for(unsigned int i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble * likelihoodArray_i = & likelihoodArray[i];
      Vdouble * leavesLikelihoods_node_i = & (* leavesLikelihoods_node)[i];
      likelihoodArray_i->resize(nbClasses_);
      for(unsigned int c = 0; c < nbClasses_; c++)
      {
        Vdouble * likelihoodArray_i_c = & (* likelihoodArray_i)[c];
        likelihoodArray_i_c->resize(nbStates_);
        for(unsigned int x = 0; x < nbStates_; x++)
        {
          (* likelihoodArray_i_c)[x] = (* leavesLikelihoods_node_i)[x];
        }
      }
    }
  }
  else
  {
    // Otherwise:
    // Set all likelihoods to 1 for a start:
    for(unsigned int i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble * likelihoodArray_i = & likelihoodArray[i];
      likelihoodArray_i->resize(nbClasses_);
      for(unsigned int c = 0; c < nbClasses_; c++)
      {
        Vdouble * likelihoodArray_i_c = & (* likelihoodArray_i)[c];
        likelihoodArray_i_c->resize(nbStates_);
        for(unsigned int x = 0; x < nbStates_; x++)
        {
          (* likelihoodArray_i_c)[x] = 1.;
        }
      }
    }
  }
  
  unsigned int nbNodes = node->getNumberOfSons();
  
  vector<const VVVdouble *> iLik(nbNodes);
  vector<const VVVdouble *> tProb(nbNodes);
  for(unsigned int n = 0; n < nbNodes; n++)
  {
    const Node * son = node->getSon(n);
    tProb[n] = & pxy_[son->getId()];
    iLik[n] = & (* likelihoods_node)[son->getId()];
  }
  
  if(node->hasFather())
  {
    const Node * father = node->getFather();
    computeLikelihoodFromArrays(iLik, tProb, & (* likelihoods_node)[father->getId()], & pxy_[nodeId], likelihoodArray, nbNodes, nbDistinctSites_, nbClasses_, nbStates_, false);
  }
  else
  {
    computeLikelihoodFromArrays(iLik, tProb, likelihoodArray, nbNodes, nbDistinctSites_, nbClasses_, nbStates_, false);
    
    //We have to account for the equilibrium frequencies:
    for(unsigned int i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble * likelihoodArray_i = & likelihoodArray[i];
      for(unsigned int c = 0; c < nbClasses_; c++)
      {
        Vdouble * likelihoodArray_i_c = & (* likelihoodArray_i)[c];
        for(unsigned int x = 0; x < nbStates_; x++)
        {
          (* likelihoodArray_i_c)[x] *= rootFreqs_[x];
        }
      }
    }
  }
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeLikelihoodFromArrays(
    const vector<const VVVdouble *> & iLik,
    const vector<const VVVdouble *> & tProb,
    VVVdouble & oLik,
    unsigned int nbNodes,
    unsigned int nbDistinctSites,
    unsigned int nbClasses,
    unsigned int nbStates,
    bool reset)
{
  if(reset) resetLikelihoodArray(oLik);

  for(unsigned int n = 0; n < nbNodes; n++)
  {
    const VVVdouble * pxy_n = tProb[n];
    const VVVdouble * iLik_n = iLik[n];

    for(unsigned int i = 0; i < nbDistinctSites; i++)
    {
      //For each site in the sequence,
      const VVdouble * iLik_n_i = & (* iLik_n)[i];
      VVdouble * oLik_i = & (oLik)[i];

      for(unsigned int c = 0; c < nbClasses; c++)
      {
        //For each rate classe,
        const Vdouble * iLik_n_i_c = & (* iLik_n_i)[c];
        Vdouble * oLik_i_c = & (* oLik_i)[c];
        const VVdouble * pxy_n_c = & (* pxy_n)[c];
        for(unsigned int x = 0; x < nbStates; x++)
        {
          //For each initial state,
          const Vdouble * pxy_n_c_x = & (* pxy_n_c)[x];
          double likelihood = 0;
          for(unsigned int y = 0; y < nbStates; y++)
          {
            likelihood += (* pxy_n_c_x)[y] * (* iLik_n_i_c)[y];
          }
          // We store this conditionnal likelihood into the corresponding array:
          (* oLik_i_c)[x] *= likelihood;
        }
      }
    }
  }
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeLikelihoodFromArrays(
    const vector<const VVVdouble *> & iLik,
    const vector<const VVVdouble *> & tProb,
    const VVVdouble * iLikR,
    const VVVdouble * tProbR,
    VVVdouble & oLik,
    unsigned int nbNodes,
    unsigned int nbDistinctSites,
    unsigned int nbClasses,
    unsigned int nbStates,
    bool reset)
{
  if(reset) resetLikelihoodArray(oLik);

  for(unsigned int n = 0; n < nbNodes; n++)
  {
    const VVVdouble * pxy_n = tProb[n];
    const VVVdouble * iLik_n = iLik[n];

    for(unsigned int i = 0; i < nbDistinctSites; i++)
    {
      //For each site in the sequence,
      const VVdouble * iLik_n_i = & (* iLik_n)[i];
      VVdouble * oLik_i = & (oLik)[i];

      for(unsigned int c = 0; c < nbClasses; c++)
      {
        //For each rate classe,
        const Vdouble * iLik_n_i_c = & (* iLik_n_i)[c];
        Vdouble * oLik_i_c = & (* oLik_i)[c];
        const VVdouble * pxy_n_c = & (* pxy_n)[c];
        for(unsigned int x = 0; x < nbStates; x++)
        {
          //For each initial state,
          const Vdouble * pxy_n_c_x = & (* pxy_n_c)[x];
          double likelihood = 0;
          for(unsigned int y = 0; y < nbStates; y++)
          {
            //cout << "1:" << (* pxy_n_c_x)[y]  << endl;
            //cout << "2:" << (* iLik_n_i_c)[y] << endl;
            likelihood += (* pxy_n_c_x)[y] * (* iLik_n_i_c)[y];
            //cout << i << "\t" << c << "\t" << x << "\t" << y << "\t" <<  (* pxy__son_c_x)[y] << "\t" << (* likelihoods_root_son_i_c)[y] << endl;
          }
          // We store this conditionnal likelihood into the corresponding array:
          (* oLik_i_c)[x] *= likelihood;
        }
      }
    }
  }

  // Now deal with the subtree containing the root:
  for(unsigned int i = 0; i < nbDistinctSites; i++)
  {
    //For each site in the sequence,
    const VVdouble * iLikR_i = & (* iLikR)[i];
    VVdouble * oLik_i = & (oLik)[i];

    for(unsigned int c = 0; c < nbClasses; c++)
    {
      //For each rate classe,
      const Vdouble * iLikR_i_c = & (* iLikR_i)[c];
      Vdouble * oLik_i_c = & (* oLik_i)[c];
      const VVdouble * pxyR_c = & (* tProbR)[c];
      for(unsigned int x = 0; x < nbStates; x++)
      {
        double likelihood = 0;
        for(unsigned int y = 0; y < nbStates; y++)
        {
          //For each final state,
          const Vdouble * pxyR_c_y = & (* pxyR_c)[y];
          likelihood += (* pxyR_c_y)[x] * (* iLikR_i_c)[y];
        }
        // We store this conditionnal likelihood into the corresponding array:
        (* oLik_i_c)[x] *= likelihood;
      }
    }
  }
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::displayLikelihood(const Node * node)
{
  cout << "Likelihoods at node " << node->getId() << ": " << endl;
  for(unsigned int n = 0; n < node->getNumberOfSons(); n++)
  {
    const Node * subNode = node->getSon(n);
    cout << "Array for sub-node " << subNode->getId() << endl;
    displayLikelihoodArray(_likelihoodData->getLikelihoodArray(node->getId(), subNode->getId()));
  }
  if(node->hasFather())
  {
    const Node * father = node->getFather();
    cout << "Array for father node " << father->getId() << endl;
    displayLikelihoodArray(_likelihoodData->getLikelihoodArray(node->getId(), father->getId()));
  }
  cout << "                                         ***" << endl;
}

/*******************************************************************************/

