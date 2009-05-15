//
// File: DRNonHomogeneousTreeLikelihood.cpp
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

#include "DRNonHomogeneousTreeLikelihood.h"
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

DRNonHomogeneousTreeLikelihood::DRNonHomogeneousTreeLikelihood(
  const Tree & tree,
  SubstitutionModelSet * modelSet,
  DiscreteDistribution * rDist,
  bool verbose)
throw (Exception):
  AbstractNonHomogeneousTreeLikelihood(tree, modelSet, rDist, verbose),
  _likelihoodData(NULL)
{
  if(!modelSet->isFullySetUpFor(tree))
    throw Exception("DRNonHomogeneousTreeLikelihood(constructor). Model set is not fully specified.");
  _init();
}

/******************************************************************************/

DRNonHomogeneousTreeLikelihood::DRNonHomogeneousTreeLikelihood(
  const Tree & tree,
  const SiteContainer & data,
  SubstitutionModelSet * modelSet,
  DiscreteDistribution * rDist,
  bool verbose)
throw (Exception):
  AbstractNonHomogeneousTreeLikelihood(tree, modelSet, rDist, verbose),
  _likelihoodData(NULL)
{
  if(!modelSet->isFullySetUpFor(tree))
    throw Exception("DRNonHomogeneousTreeLikelihood(constructor). Model set is not fully specified.");
  _init(); 
  setData(data);
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::_init() throw (Exception)
{
  _likelihoodData = new DRASDRTreeLikelihoodData(*_tree, _rateDistribution->getNumberOfCategories());
}

/******************************************************************************/

DRNonHomogeneousTreeLikelihood::DRNonHomogeneousTreeLikelihood(const DRNonHomogeneousTreeLikelihood & lik):
  AbstractNonHomogeneousTreeLikelihood(lik),
  _likelihoodData(NULL)
{
  _likelihoodData = dynamic_cast<DRASDRTreeLikelihoodData *>(lik._likelihoodData->clone());
  _likelihoodData->setTree(*_tree);
}

/******************************************************************************/

DRNonHomogeneousTreeLikelihood & DRNonHomogeneousTreeLikelihood::operator=(const DRNonHomogeneousTreeLikelihood & lik)
{
  AbstractNonHomogeneousTreeLikelihood::operator=(lik);
  if(_likelihoodData) delete _likelihoodData;
  _likelihoodData = dynamic_cast<DRASDRTreeLikelihoodData *>(lik._likelihoodData->clone());
  _likelihoodData->setTree(*_tree);
  return *this;
}

/******************************************************************************/

DRNonHomogeneousTreeLikelihood::~DRNonHomogeneousTreeLikelihood()
{
  delete _likelihoodData;
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::setData(const SiteContainer & sites) throw (Exception)
{
  if(_data) delete _data;
  _data = PatternTools::getSequenceSubset(sites, *_tree->getRootNode());
  if(_verbose) ApplicationTools::displayTask("Initializing data structure");
  _likelihoodData->initLikelihoods(* _data, *_modelSet->getModel(0)); //We assume here that all models have the same number of states, and that they have the same 'init' method,
                                                                      //Which is a reasonable assumption as long as they share the same alphabet.
  if(_verbose) ApplicationTools::displayTaskDone();

  _nbSites = _likelihoodData->getNumberOfSites();
  _nbDistinctSites = _likelihoodData->getNumberOfDistinctSites();
  _nbStates = _likelihoodData->getNumberOfStates();
  
  if(_verbose) ApplicationTools::displayResult("Number of distinct sites",
      TextTools::toString(_nbDistinctSites));
  _initialized = false;
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getLikelihood() const
{
  double l = 1.;
  Vdouble * lik = & _likelihoodData->getRootRateSiteLikelihoodArray(); 
  const vector<unsigned int> * w = & _likelihoodData->getWeights();
  for(unsigned int i = 0; i < _nbDistinctSites; i++)
  {
    l *= std::pow((*lik)[i], (int)(* w)[i]);
  }
  return l;
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getLogLikelihood() const
{
  double ll = 0;
  Vdouble * lik = & _likelihoodData->getRootRateSiteLikelihoodArray(); 
  const vector<unsigned int> * w = & _likelihoodData->getWeights();
  vector<double> la(_nbDistinctSites);
  for(unsigned int i = 0; i < _nbDistinctSites; i++)
  {
    la[i] = (* w)[i] * log((* lik)[i]);
  }
  sort(la.begin(), la.end());
  for(unsigned int i = _nbDistinctSites; i > 0; i--)
    ll += la[i-1];
  return ll;
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getLikelihoodForASite(unsigned int site) const
{
  return _likelihoodData->getRootRateSiteLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)];
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getLogLikelihoodForASite(unsigned int site) const
{
  return log(_likelihoodData->getRootRateSiteLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)]);
}

/******************************************************************************/
double DRNonHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
  return _likelihoodData->getRootSiteLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)][rateClass];
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
  return log(_likelihoodData->getRootSiteLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)][rateClass]);
}

/******************************************************************************/  

double DRNonHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
  return _likelihoodData->getRootLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)][rateClass][state];
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const
{
  return log(_likelihoodData->getRootLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)][rateClass][state]);
}

/******************************************************************************/  

void DRNonHomogeneousTreeLikelihood::setParameters(const ParameterList & parameters)
  throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::fireParameterChanged(const ParameterList & params)
{
  applyParameters();

  if(params.getCommonParametersWith(_rateDistribution->getIndependentParameters()).size() > 0)
  {
    computeAllTransitionProbabilities();
  }
  else
  {
    vector<int> ids;
    vector<string> tmp = params.getCommonParametersWith(_modelSet->getNodeParameters()).getParameterNames();
    for(unsigned int i = 0; i < tmp.size(); i++)
    {
      vector<int> tmpv = _modelSet->getNodesWithParameter(tmp[i]);
      ids = VectorTools::vectorUnion(ids, tmpv);
    }
    tmp = params.getCommonParametersWith(_brLenParameters).getParameterNames();
    vector<const Node *> nodes;
    for(unsigned int i = 0; i < ids.size(); i++)
    {
      nodes.push_back(_idToNode[ids[i]]);
    }
    vector<const Node *> tmpv;
    bool test = false;
    for(unsigned int i = 0; i < tmp.size(); i++)
    {
      if(tmp[i] == "BrLenRoot" || tmp[i] == "RootPosition")
      {
        if(!test)
        {
          tmpv.push_back(_tree->getRootNode()->getSon(0));
          tmpv.push_back(_tree->getRootNode()->getSon(1));
          test = true; //Add only once.
        }
      }
      else
        tmpv.push_back(_nodes[TextTools::to<unsigned int>(tmp[i].substr(5))]);
    }
    nodes = VectorTools::vectorUnion(nodes, tmpv);

    for(unsigned int i = 0; i < nodes.size(); i++)
    {
      computeTransitionProbabilitiesForNode(nodes[i]);
    }
    _rootFreqs = _modelSet->getRootFrequencies();
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

double DRNonHomogeneousTreeLikelihood::getValue() const
throw (Exception)
{
  if(!isInitialized()) throw Exception("DRNonHomogeneousTreeLikelihood::getValue(). Instance is not initialized.");
  return - getLogLikelihood();
}

/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/  

void DRNonHomogeneousTreeLikelihood::computeTreeDLikelihoodAtNode(const Node * node)
{
  const Node * father = node->getFather();
  VVVdouble * _likelihoods_father_node = & _likelihoodData->getLikelihoodArray(father->getId(), node->getId());
  Vdouble * _dLikelihoods_node = & _likelihoodData->getDLikelihoodArray(node->getId());
  VVVdouble *  _pxy_node = &  _pxy[node->getId()];
  VVVdouble * _dpxy_node = & _dpxy[node->getId()];
  VVVdouble larray;
  computeLikelihoodAtNode(father->getId(), larray);
  Vdouble * rootLikelihoodsSR = & _likelihoodData->getRootRateSiteLikelihoodArray();

  double dLi, dLic, dLicx, numerator, denominator;
  for(unsigned int i = 0; i < _nbDistinctSites; i++)
  {
    VVdouble * _likelihoods_father_node_i = & (* _likelihoods_father_node)[i];
    VVdouble * larray_i = & larray[i];
    dLi = 0;
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      Vdouble * _likelihoods_father_node_i_c = & (* _likelihoods_father_node_i)[c];
      Vdouble * larray_i_c = & (* larray_i)[c];
      VVdouble *  _pxy_node_c = & (*  _pxy_node)[c];
      VVdouble * _dpxy_node_c = & (* _dpxy_node)[c];
      dLic = 0;
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        numerator = 0;
        denominator = 0;
        Vdouble *  _pxy_node_c_x = & (*  _pxy_node_c)[x];
        Vdouble * _dpxy_node_c_x = & (* _dpxy_node_c)[x];
        dLicx = 0;
        for(unsigned int y = 0; y < _nbStates; y++)
        {
          numerator   += (* _dpxy_node_c_x)[y] * (* _likelihoods_father_node_i_c)[y];
          denominator += (*  _pxy_node_c_x)[y] * (* _likelihoods_father_node_i_c)[y];
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

void DRNonHomogeneousTreeLikelihood::computeTreeDLikelihoods()
{
  for(unsigned int k = 0; k < _nbNodes; k++)
  {
    computeTreeDLikelihoodAtNode(_nodes[k]);
  }
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getFirstOrderDerivative(const string & variable) const
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
  const vector<unsigned int> * w = & _likelihoodData->getWeights();
  Vdouble * _dLikelihoods_branch;
  if(variable == "BrLenRoot")
  {
    _dLikelihoods_branch = & _likelihoodData->getDLikelihoodArray(_root1);
    double d1 = 0;
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
      d1 -= (* w)[i] * (* _dLikelihoods_branch)[i];
    _dLikelihoods_branch = & _likelihoodData->getDLikelihoodArray(_root2);
    double d2 = 0;
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
      d2 -= (* w)[i] * (* _dLikelihoods_branch)[i];
    double pos = getParameterValue("RootPosition");
    return pos * d1 + (1. - pos) * d2;
  }
  else if(variable == "RootPosition")
  {
    _dLikelihoods_branch = & _likelihoodData->getDLikelihoodArray(_root1);
    double d1 = 0;
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
      d1 -= (* w)[i] * (* _dLikelihoods_branch)[i];
    _dLikelihoods_branch = & _likelihoodData->getDLikelihoodArray(_root2);
    double d2 = 0;
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
      d2 -= (* w)[i] * (* _dLikelihoods_branch)[i];
    double len = getParameterValue("BrLenRoot");
    return len * (d1 - d2);
  }
  else
  {
    // Get the node with the branch whose length must be derivated:
    unsigned int brI = TextTools::to<unsigned int>(variable.substr(5));
    const Node * branch = _nodes[brI];
    _dLikelihoods_branch = & _likelihoodData->getDLikelihoodArray(branch->getId());
    double d = 0;
    const vector<unsigned int> * w = & _likelihoodData->getWeights();
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
      d += (* w)[i] * (* _dLikelihoods_branch)[i];
    return -d;
  }  
}

/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/  

void DRNonHomogeneousTreeLikelihood::computeTreeD2LikelihoodAtNode(const Node * node)
{
  const Node * father = node->getFather();
  VVVdouble * _likelihoods_father_node = & _likelihoodData->getLikelihoodArray(father->getId(), node->getId());
  Vdouble * _d2Likelihoods_node = & _likelihoodData->getD2LikelihoodArray(node->getId());  
  VVVdouble *   _pxy_node = &   _pxy[node->getId()];
  VVVdouble * _d2pxy_node = & _d2pxy[node->getId()];
  VVVdouble larray;
  computeLikelihoodAtNode(father->getId(), larray);
  Vdouble * rootLikelihoodsSR = & _likelihoodData->getRootRateSiteLikelihoodArray();
  
  double d2Li, d2Lic, d2Licx, numerator, denominator;

  for(unsigned int i = 0; i < _nbDistinctSites; i++)
  {
    VVdouble * _likelihoods_father_node_i = & (* _likelihoods_father_node)[i];
    VVdouble * larray_i = & larray[i];
    d2Li = 0;
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      Vdouble * _likelihoods_father_node_i_c = & (* _likelihoods_father_node_i)[c];
      Vdouble * larray_i_c = & (* larray_i)[c];
      VVdouble *   _pxy_node_c = & (*   _pxy_node)[c];
      VVdouble * _d2pxy_node_c = & (* _d2pxy_node)[c];
      d2Lic = 0;
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        numerator = 0;
        denominator = 0;
        Vdouble *   _pxy_node_c_x = & (*   _pxy_node_c)[x];
        Vdouble * _d2pxy_node_c_x = & (* _d2pxy_node_c)[x];
        d2Licx = 0;
        for(unsigned int y = 0; y < _nbStates; y++)
        {
          numerator   += (* _d2pxy_node_c_x)[y] * (* _likelihoods_father_node_i_c)[y];
          denominator += (*   _pxy_node_c_x)[y] * (* _likelihoods_father_node_i_c)[y];
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

void DRNonHomogeneousTreeLikelihood::computeTreeD2Likelihoods()
{
  for(unsigned int k = 0; k < _nbNodes; k++)
  {
    computeTreeD2LikelihoodAtNode(_nodes[k]);
  }
}

/******************************************************************************/

double DRNonHomogeneousTreeLikelihood::getSecondOrderDerivative(const string & variable) const
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
  
  const vector<unsigned int> * w = & _likelihoodData->getWeights();
  //We can't deduce second order derivatives regarding BrLenRoot and RootPosition from the
  //branch length derivatives. We need a bit more calculations...
  //NB: we could save a few calculations here...
  if(variable == "BrLenRoot")
  { 
    const Node* father = _tree->getRootNode();
     
    // Compute dLikelihoods array for the father node.
    // Fist initialize to 1:
    VVVdouble dLikelihoods_father(_nbDistinctSites); 
    VVVdouble d2Likelihoods_father(_nbDistinctSites); 
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      VVdouble * dLikelihoods_father_i = & dLikelihoods_father[i];
      VVdouble * d2Likelihoods_father_i = & d2Likelihoods_father[i];
      dLikelihoods_father_i->resize(_nbClasses);
      d2Likelihoods_father_i->resize(_nbClasses);
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        Vdouble * dLikelihoods_father_i_c = & (* dLikelihoods_father_i)[c];
        Vdouble * d2Likelihoods_father_i_c = & (* d2Likelihoods_father_i)[c];
        dLikelihoods_father_i_c->resize(_nbStates);
        d2Likelihoods_father_i_c->resize(_nbStates);
        for(unsigned int s = 0; s < _nbStates; s++)
        {
          (* dLikelihoods_father_i_c)[s] = 1.;  
          (* d2Likelihoods_father_i_c)[s] = 1.;  
        }
      }
    }

    unsigned int nbNodes = father->getNumberOfSons();
    for(unsigned int l = 0; l < nbNodes; l++)
    {
      const Node * son = father->getSon(l);
    
      if(son->getId() == _root1)
      {
        VVVdouble * _likelihoods_root1 = & _likelihoodData->getLikelihoodArray(father->getId(), _root1);
        VVVdouble * _likelihoods_root2 = & _likelihoodData->getLikelihoodArray(father->getId(), _root2);
        double pos = getParameterValue("RootPosition");

        VVVdouble * _d2pxy_root1 = & _d2pxy[_root1];
        VVVdouble * _d2pxy_root2 = & _d2pxy[_root2];
        VVVdouble * _dpxy_root1  = & _dpxy[_root1];
        VVVdouble * _dpxy_root2  = & _dpxy[_root2];
        VVVdouble * _pxy_root1   = & _pxy[_root1];
        VVVdouble * _pxy_root2   = & _pxy[_root2];
        for(unsigned int i = 0; i < _nbDistinctSites; i++)
        {
          VVdouble * _likelihoods_root1_i = & (* _likelihoods_root1)[i];
          VVdouble * _likelihoods_root2_i = & (* _likelihoods_root2)[i];
          VVdouble * dLikelihoods_father_i = & dLikelihoods_father[i];
          VVdouble * d2Likelihoods_father_i = & d2Likelihoods_father[i];
          for(unsigned int c = 0; c < _nbClasses; c++)
          {
            Vdouble * _likelihoods_root1_i_c = & (* _likelihoods_root1_i)[c];
            Vdouble * _likelihoods_root2_i_c = & (* _likelihoods_root2_i)[c];
            Vdouble * dLikelihoods_father_i_c = & (* dLikelihoods_father_i)[c];
            Vdouble * d2Likelihoods_father_i_c = & (* d2Likelihoods_father_i)[c];
            VVdouble * _d2pxy_root1_c = & (* _d2pxy_root1)[c];
            VVdouble * _d2pxy_root2_c = & (* _d2pxy_root2)[c];
            VVdouble * _dpxy_root1_c  = & (* _dpxy_root1)[c];
            VVdouble * _dpxy_root2_c  = & (* _dpxy_root2)[c];
            VVdouble * _pxy_root1_c   = & (* _pxy_root1)[c];
            VVdouble * _pxy_root2_c   = & (* _pxy_root2)[c];
            for(unsigned int x = 0; x < _nbStates; x++)
            {
              Vdouble * _d2pxy_root1_c_x = & (* _d2pxy_root1_c)[x];
              Vdouble * _d2pxy_root2_c_x = & (* _d2pxy_root2_c)[x];
              Vdouble * _dpxy_root1_c_x  = & (* _dpxy_root1_c)[x];
              Vdouble * _dpxy_root2_c_x  = & (* _dpxy_root2_c)[x];
              Vdouble * _pxy_root1_c_x   = & (* _pxy_root1_c)[x];
              Vdouble * _pxy_root2_c_x   = & (* _pxy_root2_c)[x];
              double d2l1 = 0, d2l2 = 0, dl1 = 0, dl2 = 0, l1 = 0, l2 = 0;
              for(unsigned int y = 0; y < _nbStates; y++)
              {
                d2l1 += (* _d2pxy_root1_c_x)[y] * (* _likelihoods_root1_i_c)[y];
                d2l2 += (* _d2pxy_root2_c_x)[y] * (* _likelihoods_root2_i_c)[y];
                dl1  += (* _dpxy_root1_c_x)[y]  * (* _likelihoods_root1_i_c)[y];
                dl2  += (* _dpxy_root2_c_x)[y]  * (* _likelihoods_root2_i_c)[y];
                l1   += (* _pxy_root1_c_x)[y]   * (* _likelihoods_root1_i_c)[y];
                l2   += (* _pxy_root2_c_x)[y]   * (* _likelihoods_root2_i_c)[y];
              }
              double dl = pos * dl1 * l2 + (1. - pos) * dl2 * l1;
              double d2l = pos * pos * d2l1 * l2 + (1. - pos) * (1. - pos) * d2l2 * l1 + 2 * pos * (1. - pos) * dl1 * dl2;
              (* dLikelihoods_father_i_c)[x] *= dl;
              (* d2Likelihoods_father_i_c)[x] *= d2l;
            }
          }
        }
      }
      else if(son->getId() == _root2)
      {
        //Do nothing, this was accounted when dealing with _root1 
      }
      else
      {
        //Account for a putative multifurcation:
        VVVdouble * _likelihoods_son = & _likelihoodData->getLikelihoodArray(father->getId(), son->getId());

        VVVdouble * _pxy_son = & _pxy[son->getId()];
        for(unsigned int i = 0; i < _nbDistinctSites; i++)
        {
          VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[i];
          VVdouble * dLikelihoods_father_i = & dLikelihoods_father[i];
          VVdouble * d2Likelihoods_father_i = & d2Likelihoods_father[i];
          for(unsigned int c = 0; c < _nbClasses; c++)
          {
            Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
            Vdouble * dLikelihoods_father_i_c = & (* dLikelihoods_father_i)[c];
            Vdouble * d2Likelihoods_father_i_c = & (* d2Likelihoods_father_i)[c];
            VVdouble * _pxy_son_c = & (* _pxy_son)[c];
            for(unsigned int x = 0; x < _nbStates; x++)
            {
              double dl = 0;
              Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
              for(unsigned int y = 0; y < _nbStates; y++)
              {
                dl += (* _pxy_son_c_x)[y] * (* _likelihoods_son_i_c)[y];
              }
              (* dLikelihoods_father_i_c)[x] *= dl;
              (* d2Likelihoods_father_i_c)[x] *= dl;
            }
          }
        }
      }  
    }
    Vdouble * rootLikelihoodsSR = & _likelihoodData->getRootRateSiteLikelihoodArray();
    double d2l = 0, dlx, d2lx;
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      VVdouble * dLikelihoods_father_i = & dLikelihoods_father[i];
      VVdouble * d2Likelihoods_father_i = & d2Likelihoods_father[i];
      dlx = 0, d2lx = 0;
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        Vdouble * dLikelihoods_father_i_c = & (* dLikelihoods_father_i)[c];
        Vdouble * d2Likelihoods_father_i_c = & (* d2Likelihoods_father_i)[c];
        for(unsigned int x = 0; x < _nbStates; x++)
        {
          dlx += _rateDistribution->getProbability(c) * _rootFreqs[x] * (* dLikelihoods_father_i_c)[x];
          d2lx += _rateDistribution->getProbability(c) * _rootFreqs[x] * (* d2Likelihoods_father_i_c)[x];
        }
      }
      d2l += (* w)[i] * (d2lx / (* rootLikelihoodsSR)[i] - pow(dlx / (* rootLikelihoodsSR)[i], 2));
    }
    return - d2l;
  }
  else if(variable == "RootPosition")
  {
    const Node* father = _tree->getRootNode();
     
    // Compute dLikelihoods array for the father node.
    // Fist initialize to 1:
    VVVdouble dLikelihoods_father(_nbDistinctSites); 
    VVVdouble d2Likelihoods_father(_nbDistinctSites); 
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      VVdouble * dLikelihoods_father_i = & dLikelihoods_father[i];
      VVdouble * d2Likelihoods_father_i = & d2Likelihoods_father[i];
      dLikelihoods_father_i->resize(_nbClasses);
      d2Likelihoods_father_i->resize(_nbClasses);
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        Vdouble * dLikelihoods_father_i_c = & (* dLikelihoods_father_i)[c];
        Vdouble * d2Likelihoods_father_i_c = & (* d2Likelihoods_father_i)[c];
        dLikelihoods_father_i_c->resize(_nbStates);
        d2Likelihoods_father_i_c->resize(_nbStates);
        for(unsigned int s = 0; s < _nbStates; s++)
        {
          (* dLikelihoods_father_i_c)[s] = 1.;  
          (* d2Likelihoods_father_i_c)[s] = 1.;  
        }
      }
    }

    unsigned int nbNodes = father->getNumberOfSons();
    for(unsigned int l = 0; l < nbNodes; l++)
    {
      const Node * son = father->getSon(l);
    
      if(son->getId() == _root1)
      {
        VVVdouble * _likelihoods_root1 = & _likelihoodData->getLikelihoodArray(father->getId(), _root1);
        VVVdouble * _likelihoods_root2 = & _likelihoodData->getLikelihoodArray(father->getId(), _root2);
        double len = getParameterValue("BrLenRoot");

        VVVdouble * _d2pxy_root1 = & _d2pxy[_root1];
        VVVdouble * _d2pxy_root2 = & _d2pxy[_root2];
        VVVdouble * _dpxy_root1  = & _dpxy[_root1];
        VVVdouble * _dpxy_root2  = & _dpxy[_root2];
        VVVdouble * _pxy_root1   = & _pxy[_root1];
        VVVdouble * _pxy_root2   = & _pxy[_root2];
        for(unsigned int i = 0; i < _nbDistinctSites; i++)
        {
          VVdouble * _likelihoods_root1_i = & (* _likelihoods_root1)[i];
          VVdouble * _likelihoods_root2_i = & (* _likelihoods_root2)[i];
          VVdouble * dLikelihoods_father_i = & dLikelihoods_father[i];
          VVdouble * d2Likelihoods_father_i = & d2Likelihoods_father[i];
          for(unsigned int c = 0; c < _nbClasses; c++)
          {
            Vdouble * _likelihoods_root1_i_c = & (* _likelihoods_root1_i)[c];
            Vdouble * _likelihoods_root2_i_c = & (* _likelihoods_root2_i)[c];
            Vdouble * dLikelihoods_father_i_c = & (* dLikelihoods_father_i)[c];
            Vdouble * d2Likelihoods_father_i_c = & (* d2Likelihoods_father_i)[c];
            VVdouble * _d2pxy_root1_c = & (* _d2pxy_root1)[c];
            VVdouble * _d2pxy_root2_c = & (* _d2pxy_root2)[c];
            VVdouble * _dpxy_root1_c  = & (* _dpxy_root1)[c];
            VVdouble * _dpxy_root2_c  = & (* _dpxy_root2)[c];
            VVdouble * _pxy_root1_c   = & (* _pxy_root1)[c];
            VVdouble * _pxy_root2_c   = & (* _pxy_root2)[c];
            for(unsigned int x = 0; x < _nbStates; x++)
            {
              Vdouble * _d2pxy_root1_c_x = & (* _d2pxy_root1_c)[x];
              Vdouble * _d2pxy_root2_c_x = & (* _d2pxy_root2_c)[x];
              Vdouble * _dpxy_root1_c_x  = & (* _dpxy_root1_c)[x];
              Vdouble * _dpxy_root2_c_x  = & (* _dpxy_root2_c)[x];
              Vdouble * _pxy_root1_c_x   = & (* _pxy_root1_c)[x];
              Vdouble * _pxy_root2_c_x   = & (* _pxy_root2_c)[x];
              double d2l1 = 0, d2l2 = 0, dl1 = 0, dl2 = 0, l1 = 0, l2 = 0;
              for(unsigned int y = 0; y < _nbStates; y++)
              {
                d2l1 += (* _d2pxy_root1_c_x)[y] * (* _likelihoods_root1_i_c)[y];
                d2l2 += (* _d2pxy_root2_c_x)[y] * (* _likelihoods_root2_i_c)[y];
                dl1  += (* _dpxy_root1_c_x)[y]  * (* _likelihoods_root1_i_c)[y];
                dl2  += (* _dpxy_root2_c_x)[y]  * (* _likelihoods_root2_i_c)[y];
                l1   += (* _pxy_root1_c_x)[y]   * (* _likelihoods_root1_i_c)[y];
                l2   += (* _pxy_root2_c_x)[y]   * (* _likelihoods_root2_i_c)[y];
              }
              double dl = len * (dl1 * l2 - dl2 * l1);
              double d2l = len * len * (d2l1 * l2 + d2l2 * l1 - 2 * dl1 * dl2);
              (* dLikelihoods_father_i_c)[x] *= dl;
              (* d2Likelihoods_father_i_c)[x] *= d2l;
            }
          }
        }
      }
      else if(son->getId() == _root2)
      {
        //Do nothing, this was accounted when dealing with _root1 
      }
      else
      {
        //Account for a putative multifurcation:
        VVVdouble * _likelihoods_son = & _likelihoodData->getLikelihoodArray(father->getId(), son->getId());

        VVVdouble * _pxy_son = & _pxy[son->getId()];
        for(unsigned int i = 0; i < _nbDistinctSites; i++)
        {
          VVdouble * _likelihoods_son_i = & (* _likelihoods_son)[i];
          VVdouble * dLikelihoods_father_i = & dLikelihoods_father[i];
          VVdouble * d2Likelihoods_father_i = & d2Likelihoods_father[i];
          for(unsigned int c = 0; c < _nbClasses; c++)
          {
            Vdouble * _likelihoods_son_i_c = & (* _likelihoods_son_i)[c];
            Vdouble * dLikelihoods_father_i_c = & (* dLikelihoods_father_i)[c];
            Vdouble * d2Likelihoods_father_i_c = & (* d2Likelihoods_father_i)[c];
            VVdouble * _pxy_son_c = & (* _pxy_son)[c];
            for(unsigned int x = 0; x < _nbStates; x++)
            {
              double dl = 0;
              Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
              for(unsigned int y = 0; y < _nbStates; y++)
              {
                dl += (* _pxy_son_c_x)[y] * (* _likelihoods_son_i_c)[y];
              }
              (* dLikelihoods_father_i_c)[x] *= dl;
              (* d2Likelihoods_father_i_c)[x] *= dl;
            }
          }
        }
      }  
    }
    Vdouble * rootLikelihoodsSR = & _likelihoodData->getRootRateSiteLikelihoodArray();
    double d2l = 0, dlx, d2lx;
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      VVdouble * dLikelihoods_father_i = & dLikelihoods_father[i];
      VVdouble * d2Likelihoods_father_i = & d2Likelihoods_father[i];
      dlx = 0, d2lx = 0;
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        Vdouble * dLikelihoods_father_i_c = & (* dLikelihoods_father_i)[c];
        Vdouble * d2Likelihoods_father_i_c = & (* d2Likelihoods_father_i)[c];
        for(unsigned int x = 0; x < _nbStates; x++)
        {
          dlx += _rateDistribution->getProbability(c) * _rootFreqs[x] * (* dLikelihoods_father_i_c)[x];
          d2lx += _rateDistribution->getProbability(c) * _rootFreqs[x] * (* d2Likelihoods_father_i_c)[x];
        }
      }
      d2l += (* w)[i] * (d2lx / (* rootLikelihoodsSR)[i] - pow(dlx / (* rootLikelihoodsSR)[i], 2));
    }
    return - d2l;
  }
  else
  {
    Vdouble * _dLikelihoods_branch;
    Vdouble * _d2Likelihoods_branch;
    // Get the node with the branch whose length must be derivated:
    unsigned int brI = TextTools::to<unsigned int>(variable.substr(5));
    const Node * branch = _nodes[brI];
    _dLikelihoods_branch = & _likelihoodData->getDLikelihoodArray(branch->getId());
    _d2Likelihoods_branch = & _likelihoodData->getD2LikelihoodArray(branch->getId());
    double d2l = 0;
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
      d2l += (* w)[i] * ((* _d2Likelihoods_branch)[i] - pow((* _dLikelihoods_branch)[i], 2));
    return - d2l;
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::resetLikelihoodArrays(const Node * node)
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
  
void DRNonHomogeneousTreeLikelihood::computeTreeLikelihood()
{
  computeSubtreeLikelihoodPostfix(_tree->getRootNode());
  computeSubtreeLikelihoodPrefix(_tree->getRootNode());
  computeRootLikelihood();
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::computeSubtreeLikelihoodPostfix(const Node * node)
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
      for(unsigned int i = 0; i < _nbDistinctSites; i++)
      {
        //For each site in the sequence,
        Vdouble * _likelihoods_leaf_i = & (* _likelihoods_leaf)[i];
        VVdouble * _likelihoods_node_son_i = & (* _likelihoods_node_son)[i];
        for(unsigned int c = 0; c < _nbClasses; c++)
        {
          //For each rate classe,
          Vdouble * _likelihoods_node_son_i_c = & (* _likelihoods_node_son_i)[c];
          for(unsigned int x = 0; x < _nbStates; x++)
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
        tProb[n] = & _pxy[sonSon->getId()];
        iLik[n] = & (* _likelihoods_son)[sonSon->getId()];
      }
      computeLikelihoodFromArrays(iLik, tProb, *_likelihoods_node_son, nbSons, _nbDistinctSites, _nbClasses, _nbStates, false); 
    }
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::computeSubtreeLikelihoodPrefix(const Node * node)
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
      for(unsigned int i = 0; i < _nbDistinctSites; i++)
      {
        //For each site in the sequence,
        Vdouble * _likelihoods_leaf_i = & (* _likelihoods_leaf)[i];
        VVdouble * _likelihoods_node_father_i = & (* _likelihoods_node_father)[i];
        for(unsigned int c = 0; c < _nbClasses; c++)
        {
          //For each rate classe,
          Vdouble * _likelihoods_node_father_i_c = & (* _likelihoods_node_father_i)[c];
          for(unsigned int x = 0; x < _nbStates; x++)
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
  
      unsigned int nbSons = nodes.size(); // In case of a bifurcating tree this is equal to 1.
      
      vector<const VVVdouble *> iLik(nbSons);
      vector<const VVVdouble *> tProb(nbSons);
      for(unsigned int n = 0; n < nbSons; n++)
      {
        const Node * fatherSon = nodes[n];
        tProb[n] = & _pxy[fatherSon->getId()];
        iLik[n] = & (* _likelihoods_father)[fatherSon->getId()];
      }
    
      if(father->hasFather())
      {
        const Node * fatherFather = father->getFather();
        computeLikelihoodFromArrays(iLik, tProb, & (* _likelihoods_father)[fatherFather->getId()], & _pxy[father->getId()], *_likelihoods_node_father, nbSons, _nbDistinctSites, _nbClasses, _nbStates, false); 
      }
      else
      {
        computeLikelihoodFromArrays(iLik, tProb, *_likelihoods_node_father, nbSons, _nbDistinctSites, _nbClasses, _nbStates, false); 
      }
    }

    if(!father->hasFather())
    {
      //We have to account for the root frequencies:
      for(unsigned int i = 0; i < _nbDistinctSites; i++)
      {
        VVdouble * _likelihoods_node_father_i = & (*_likelihoods_node_father)[i];
        for(unsigned int c = 0; c < _nbClasses; c++)
        {
          Vdouble * _likelihoods_node_father_i_c = & (* _likelihoods_node_father_i)[c];
          for(unsigned int x = 0; x < _nbStates; x++)
          {
            (* _likelihoods_node_father_i_c)[x] *= _rootFreqs[x];
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

void DRNonHomogeneousTreeLikelihood::computeRootLikelihood()
{
  const Node * root = _tree->getRootNode();
  VVVdouble * rootLikelihoods = & _likelihoodData->getRootLikelihoodArray();
  // Set all likelihoods to 1 for a start:
  if(root->isLeaf())
  {
    VVdouble * leavesLikelihoods_root = & _likelihoodData->getLeafLikelihoods(root->getId());
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      VVdouble * rootLikelihoods_i = & (* rootLikelihoods)[i];
      Vdouble * leavesLikelihoods_root_i = & (* leavesLikelihoods_root)[i];
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        Vdouble * rootLikelihoods_i_c = & (* rootLikelihoods_i)[c];
        for(unsigned int x = 0; x < _nbStates; x++)
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
    tProb[n] = & _pxy[son->getId()];
    iLik[n] = & (* likelihoods_root)[son->getId()];
  }
  computeLikelihoodFromArrays(iLik, tProb, *rootLikelihoods, nbNodes, _nbDistinctSites, _nbClasses, _nbStates, false);

  Vdouble p = _rateDistribution->getProbabilities();
  VVdouble * rootLikelihoodsS  = & _likelihoodData->getRootSiteLikelihoodArray();
  Vdouble  * rootLikelihoodsSR = & _likelihoodData->getRootRateSiteLikelihoodArray();
  for(unsigned int i = 0; i < _nbDistinctSites; i++)
  {
    //For each site in the sequence,
    VVdouble * rootLikelihoods_i = & (* rootLikelihoods)[i];
    Vdouble * rootLikelihoodsS_i = & (* rootLikelihoodsS)[i];
    (* rootLikelihoodsSR)[i] = 0;
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      //For each rate classe,
      Vdouble * rootLikelihoods_i_c = & (* rootLikelihoods_i)[c];
      double * rootLikelihoodsS_i_c = & (* rootLikelihoodsS_i)[c];
      (* rootLikelihoodsS_i_c) = 0;
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        //For each initial state,
        (* rootLikelihoodsS_i_c) += _rootFreqs[x] * (* rootLikelihoods_i_c)[x];
      }
      (* rootLikelihoodsSR)[i] += p[c] * (* rootLikelihoodsS_i_c);
    }

    //Final checking (for numerical errors):
    if((* rootLikelihoodsSR)[i] < 0) (* rootLikelihoodsSR)[i] = 0.;
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::computeLikelihoodAtNode(int nodeId, VVVdouble& likelihoodArray) const
{
  const Node * node = _tree->getNode(nodeId);

  likelihoodArray.resize(_nbDistinctSites);
  map<int, VVVdouble> * likelihoods_node = & _likelihoodData->getLikelihoodArrays(nodeId);
  
  //Initialize likelihood array:
  if(node->isLeaf())
  {
    VVdouble * leavesLikelihoods_node = & _likelihoodData->getLeafLikelihoods(nodeId);
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      VVdouble * likelihoodArray_i = & likelihoodArray[i];
      Vdouble * leavesLikelihoods_node_i = & (* leavesLikelihoods_node)[i];
      likelihoodArray_i->resize(_nbClasses);
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        Vdouble * likelihoodArray_i_c = & (* likelihoodArray_i)[c];
        likelihoodArray_i_c->resize(_nbStates);
        for(unsigned int x = 0; x < _nbStates; x++)
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
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      VVdouble * likelihoodArray_i = & likelihoodArray[i];
      likelihoodArray_i->resize(_nbClasses);
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        Vdouble * likelihoodArray_i_c = & (* likelihoodArray_i)[c];
        likelihoodArray_i_c->resize(_nbStates);
        for(unsigned int x = 0; x < _nbStates; x++)
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
    tProb[n] = & _pxy[son->getId()];
    iLik[n] = & (* likelihoods_node)[son->getId()];
  }
  
  if(node->hasFather())
  {
    const Node * father = node->getFather();
    computeLikelihoodFromArrays(iLik, tProb, & (* likelihoods_node)[father->getId()], & _pxy[nodeId], likelihoodArray, nbNodes, _nbDistinctSites, _nbClasses, _nbStates, false);
  }
  else
  {
    computeLikelihoodFromArrays(iLik, tProb, likelihoodArray, nbNodes, _nbDistinctSites, _nbClasses, _nbStates, false);
    
    //We have to account for the root frequencies:
     for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      VVdouble * likelihoodArray_i = & likelihoodArray[i];
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        Vdouble * likelihoodArray_i_c = & (* likelihoodArray_i)[c];
        for(unsigned int x = 0; x < _nbStates; x++)
        {
          (* likelihoodArray_i_c)[x] *= _rootFreqs[x];
        }
      }
    }
  }
}

/******************************************************************************/

void DRNonHomogeneousTreeLikelihood::computeLikelihoodFromArrays(
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

void DRNonHomogeneousTreeLikelihood::computeLikelihoodFromArrays(
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
            //cout << i << "\t" << c << "\t" << x << "\t" << y << "\t" <<  (* _pxy_son_c_x)[y] << "\t" << (* likelihoods_root_son_i_c)[y] << endl;
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

void DRNonHomogeneousTreeLikelihood::displayLikelihood(const Node * node)
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

