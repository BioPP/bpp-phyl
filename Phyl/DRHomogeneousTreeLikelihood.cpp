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

// From NumCalc:
using namespace VectorFunctions;

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
  //AbstractParametrizable(),
  //AbstractTreeLikelihood(),
  //AbstractDiscreteRatesAcrossSitesTreeLikelihood(rDist, verbose),
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  _likelihoodData(NULL),
  _brLikFunction(NULL),
  _brentOptimizer(NULL),
  _brLenNNIValues(),
  _brLenNNIParams()
{
  if(verbose) ApplicationTools::message << "Double-Recursive Homogeneous Tree Likelihood" << endl;  
  _likelihoodData = new DRASDRTreeLikelihoodData(*_tree, rDist->getNumberOfCategories());
    
  // Now initializes all parameters:
  initParameters();
  computeAllTransitionProbabilities();
  //fireParameterChanged(_parameters);
  
  _brentOptimizer = new BrentOneDimension();
  _brentOptimizer->setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_AUTO);
  _brentOptimizer->setProfiler(NULL);
  _brentOptimizer->setMessageHandler(NULL);
  _brentOptimizer->setVerbose(0);
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
  //AbstractParametrizable(),
  //AbstractTreeLikelihood(),
  //AbstractDiscreteRatesAcrossSitesTreeLikelihood(rDist, verbose),
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  _likelihoodData(NULL),
  _brLikFunction(NULL),
  _brentOptimizer(NULL),
  _brLenNNIValues(),
  _brLenNNIParams()
{
  if(verbose) ApplicationTools::message << "Double-Recursive Homogeneous Tree Likelihood" << endl;  
  _likelihoodData = new DRASDRTreeLikelihoodData(*_tree, rDist->getNumberOfCategories());
  
  _brentOptimizer = new BrentOneDimension();
  _brentOptimizer->setConstraintPolicy(AbstractOptimizer::CONSTRAINTS_AUTO);
  _brentOptimizer->setProfiler(NULL);
  _brentOptimizer->setMessageHandler(NULL);
  _brentOptimizer->setVerbose(0);
    
  setData(data);
  
  // Now initializes all parameters:
  initParameters();
  fireParameterChanged(_parameters);  
}

/******************************************************************************/

DRHomogeneousTreeLikelihood::DRHomogeneousTreeLikelihood(const DRHomogeneousTreeLikelihood & lik):
  //AbstractParametrizable(lik),
  //AbstractTreeLikelihood(lik),
  //AbstractDiscreteRatesAcrossSitesTreeLikelihood(lik),
  AbstractHomogeneousTreeLikelihood(lik),
  _likelihoodData(NULL)
{
  _likelihoodData = lik._likelihoodData->clone();
  _likelihoodData->setTree(*_tree);
  _brLikFunction  = dynamic_cast<BranchLikelihood *>(lik._brLikFunction->clone());
  _brentOptimizer = dynamic_cast<BrentOneDimension *>(lik._brentOptimizer->clone());
  _brLenNNIValues = lik._brLenNNIValues;
  _brLenNNIParams = lik._brLenNNIParams;
}

/******************************************************************************/

DRHomogeneousTreeLikelihood & DRHomogeneousTreeLikelihood::operator=(const DRHomogeneousTreeLikelihood & lik)
{
  AbstractParametrizable::operator=(lik);
  AbstractTreeLikelihood::operator=(lik);
  AbstractDiscreteRatesAcrossSitesTreeLikelihood::operator=(lik);
  AbstractHomogeneousTreeLikelihood::operator=(lik);
  if(_likelihoodData) delete _likelihoodData;
  _likelihoodData = lik._likelihoodData->clone();
  _likelihoodData->setTree(*_tree);
  if(_brLikFunction) delete _brLikFunction;
  _brLikFunction  = dynamic_cast<BranchLikelihood *>(lik._brLikFunction->clone());
  if(_brentOptimizer) delete _brentOptimizer;
  _brentOptimizer = dynamic_cast<BrentOneDimension *>(lik._brentOptimizer->clone());
  _brLenNNIValues = lik._brLenNNIValues;
  _brLenNNIParams = lik._brLenNNIParams;
  _parameters = lik._parameters;
  return *this;
}

/******************************************************************************/

DRHomogeneousTreeLikelihood::~DRHomogeneousTreeLikelihood()
{
  delete _likelihoodData;
  if(_brLikFunction) delete _brLikFunction;
  delete _brentOptimizer;
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::setData(const SiteContainer & sites) throw (Exception)
{
  _data = PatternTools::getSequenceSubset(sites, *_tree->getRootNode());
   if(_verbose) ApplicationTools::displayTask("Initializing data structure");
  _likelihoodData->initLikelihoods(*_data, *_model);
  if(_verbose) ApplicationTools::displayTaskDone();

  _nbSites = _likelihoodData->getNumberOfSites();
  _nbDistinctSites = _likelihoodData->getNumberOfDistinctSites();
  _nbStates = _likelihoodData->getNumberOfStates();
  
  if(_verbose) ApplicationTools::displayResult("Number of distinct sites",
      TextTools::toString(_nbDistinctSites));
  
  if(_brLikFunction) delete _brLikFunction;
  _brLikFunction = new BranchLikelihood(_likelihoodData->getWeights());
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLikelihood() const
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

double DRHomogeneousTreeLikelihood::getLogLikelihood() const
{
  double ll = 0;
  Vdouble * lik = & _likelihoodData->getRootRateSiteLikelihoodArray(); 
  const vector<unsigned int> * w = & _likelihoodData->getWeights();
  for(unsigned int i = 0; i < _nbDistinctSites; i++)
  {
    ll += (* w)[i] * log((* lik)[i]);
    //cout << i << "\t" << (* w)[i] << "\t" << log((* lik)[i]) << endl;
  }
  return ll;
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLikelihoodForASite(unsigned int site) const
{
  //double l = 0;
  //for(unsigned int i = 0; i < _nbClasses; i++)
  //{
  //  l += getLikelihoodForASiteForARateClass(site, i) * _rateDistribution->getProbability(i);
  //}
  //return l;
  return _likelihoodData->getRootRateSiteLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)];
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihoodForASite(unsigned int site) const
{
  //double l = 0;
  //for(unsigned int i = 0; i < _nbClasses; i++)
  //{
  //  l += getLikelihoodForASiteForARateClass(site, i) * _rateDistribution->getProbability(i);
  //}
  ////if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
  //return log(l);
  return log(_likelihoodData->getRootRateSiteLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)]);
}

/******************************************************************************/
double DRHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
  //double l = 0;
  //for(unsigned int i = 0; i < _nbStates; i++)
  //{
  //  l += _rootLikelihoods[_rootPatternLinks[site]][rateClass][i] * _model->freq(i);
  //}
  //return l;
  return _likelihoodData->getRootSiteLikelihoodArray()[_likelihoodData->getRootArrayPosition(site)][rateClass];
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const
{
  //double l = 0;
  //for(unsigned int i = 0; i < _nbStates; i++)
  //{
  //  l += _rootLikelihoods[_rootPatternLinks[site]][rateClass][i] * _model->freq(i);
  //}
  ////if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
  //return log(l);
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
  || _model->getParameters().getCommonParametersWith(params).size() > 0)
  {
    //Rate parameter changed, need to recompute all probs:
    computeAllTransitionProbabilities();
  }
  else if(params.size() > 0)
  {
    ////We may save some computations:
    //for(unsigned int i = 0; i < params.size(); i++)
    //{
    //  string s = params[i]->getName();
    //  cout << s << endl;
    //  if(s.substr(0,5) == "BrLen")
    //  {
    //    //Branch length parameter:
    //    computeTransitionProbabilitiesForNode(_tree->getNode(TextTools::toInt(s.substr(5))));
    //  }
    //}
  }
  computeAllTransitionProbabilities();

  computeTreeLikelihood();
  if(_computeDerivatives)
  {
    computeTreeDLikelihoods();  
    computeTreeD2Likelihoods();
  }
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getValue() const
throw (Exception)
{
  //double f = - getLogLikelihood(); // For minimization.
  //if(isnan(f)) f = -log(0.); // (+inf if unlikely!)
  //return f;
  return - getLogLikelihood();
}

/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/  

void DRHomogeneousTreeLikelihood::computeTreeDLikelihoodAtNode(const Node * node)
{
  const Node * father = node->getFather();
  VVVdouble * _likelihoods_father_node = & _likelihoodData->getLikelihoodArray(father, node);
  Vdouble * _dLikelihoods_node = & _likelihoodData->getDLikelihoodArray(node);
  VVVdouble *  _pxy_node = &  _pxy[node->getId()];
  VVVdouble * _dpxy_node = & _dpxy[node->getId()];
  VVVdouble larray = computeLikelihoodAtNode(father);
  Vdouble * rootLikelihoodsSR = & _likelihoodData->getRootRateSiteLikelihoodArray();

  for(unsigned int i = 0; i < _nbDistinctSites; i++)
  {
    VVdouble * _likelihoods_father_node_i = & (* _likelihoods_father_node)[i];
    VVdouble * larray_i = & larray[i];
    double dLi = 0;
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      Vdouble * _likelihoods_father_node_i_c = & (* _likelihoods_father_node_i)[c];
      Vdouble * larray_i_c = & (* larray_i)[c];
      VVdouble *  _pxy_node_c = & (*  _pxy_node)[c];
      VVdouble * _dpxy_node_c = & (* _dpxy_node)[c];
      double dLic = 0;
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        double numerator = 0;
        double denominator = 0;
        Vdouble *  _pxy_node_c_x = & (*  _pxy_node_c)[x];
        Vdouble * _dpxy_node_c_x = & (* _dpxy_node_c)[x];
        double dLicx = 0;
        for(unsigned int y = 0; y < _nbStates; y++)
        {
          numerator   += (* _dpxy_node_c_x)[y] * (* _likelihoods_father_node_i_c)[y];
          denominator += (*  _pxy_node_c_x)[y] * (* _likelihoods_father_node_i_c)[y];
        }
        dLicx = (* larray_i_c)[x] * numerator / denominator;
        dLic += _model->freq(x) * dLicx;  
      }
      dLi += _rateDistribution->getProbability(c) * dLic;
    }
    (* _dLikelihoods_node)[i] = dLi / (* rootLikelihoodsSR)[i]; 
  }
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeTreeDLikelihoods()
{
  for(unsigned int k = 0; k < _nbNodes; k++)
  {
    computeTreeDLikelihoodAtNode(_nodes[k]);
  }
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getFirstOrderDerivative(const string & variable) const
throw (Exception)
{ 
  Parameter * p = _parameters.getParameter(variable);
  if(p == NULL) throw ParameterNotFoundException("HomogeneousTreeLikelihood::df", variable);
  if(getRateDistributionParameters().getParameter(variable) != NULL)
  {
    throw Exception("Derivatives respective to rate distribution parameter are not implemented.");
  }
  if(getSubstitutionModelParameters().getParameter(variable) != NULL)
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }
  
  //
  // Computation for branch lengths:
  //
  
  // Get the node with the branch whose length must be derivated:
  int brI = TextTools::toInt(variable.substr(5));
  Node * branch = _nodes[brI];
  Vdouble * _dLikelihoods_branch = & _likelihoodData->getDLikelihoodArray(branch);
  double d = 0;
  //for(unsigned int i = 0; i < _nbSites; i++) d += (* _dLikelihoods_branch)[_rootPatternLinks[i]];
  const vector<unsigned int> * w = & _likelihoodData->getWeights();
  for(unsigned int i = 0; i < _nbDistinctSites; i++) d += (* w)[i] * (* _dLikelihoods_branch)[i];
  return -d;
}

/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/  

void DRHomogeneousTreeLikelihood::computeTreeD2LikelihoodAtNode(const Node * node)
{
  const Node * father = node->getFather();
  VVVdouble * _likelihoods_father_node = & _likelihoodData->getLikelihoodArray(father, node);
  Vdouble * _d2Likelihoods_node = & _likelihoodData->getD2LikelihoodArray(node);  
  VVVdouble *   _pxy_node = &   _pxy[node->getId()];
  VVVdouble * _d2pxy_node = & _d2pxy[node->getId()];
  VVVdouble larray = computeLikelihoodAtNode(father);
  Vdouble * rootLikelihoodsSR = & _likelihoodData->getRootRateSiteLikelihoodArray();
  
  for(unsigned int i = 0; i < _nbDistinctSites; i++)
  {
    VVdouble * _likelihoods_father_node_i = & (* _likelihoods_father_node)[i];
    VVdouble * larray_i = & larray[i];
    double d2Li = 0;
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      Vdouble * _likelihoods_father_node_i_c = & (* _likelihoods_father_node_i)[c];
      Vdouble * larray_i_c = & (* larray_i)[c];
      VVdouble *   _pxy_node_c = & (*   _pxy_node)[c];
      VVdouble * _d2pxy_node_c = & (* _d2pxy_node)[c];
      double d2Lic = 0;
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        double numerator = 0;
        double denominator = 0;
        Vdouble *   _pxy_node_c_x = & (*   _pxy_node_c)[x];
        Vdouble * _d2pxy_node_c_x = & (* _d2pxy_node_c)[x];
        double d2Licx = 0;
        for(unsigned int y = 0; y < _nbStates; y++)
        {
          numerator   += (* _d2pxy_node_c_x)[y] * (* _likelihoods_father_node_i_c)[y];
          denominator += (*   _pxy_node_c_x)[y] * (* _likelihoods_father_node_i_c)[y];
        }
        d2Licx = (* larray_i_c)[x] * numerator / denominator;
        d2Lic += _model->freq(x) * d2Licx;
      }
      d2Li += _rateDistribution->getProbability(c) * d2Lic;
    }
    (* _d2Likelihoods_node)[i] = d2Li / (* rootLikelihoodsSR)[i]; 
  }
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeTreeD2Likelihoods()
{
  for(unsigned int k = 0; k < _nbNodes; k++)
  {
    computeTreeD2LikelihoodAtNode(_nodes[k]);
  }
}

/******************************************************************************/

double DRHomogeneousTreeLikelihood::getSecondOrderDerivative(const string & variable) const
throw (Exception)
{
  Parameter * p = _parameters.getParameter(variable);
  if(p == NULL) throw ParameterNotFoundException("HomogeneousTreeLikelihood::df", variable);
  if(getRateDistributionParameters().getParameter(variable) != NULL)
  {
    throw Exception("Derivatives respective to rate distribution parameter are not implemented.");
  }
  if(getSubstitutionModelParameters().getParameter(variable) != NULL)
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }
  
  //
  // Computation for branch lengths:
  //
  
  // Get the node with the branch whose length must be derivated:
  int brI = TextTools::toInt(variable.substr(5));
  Node * branch = _nodes[brI];
  Vdouble * _dLikelihoods_branch = & _likelihoodData->getDLikelihoodArray(branch);
  Vdouble * _d2Likelihoods_branch = & _likelihoodData->getD2LikelihoodArray(branch);
  double d2 = 0;
  //for(unsigned int i = 0; i < _nbSites; i++) d2 += (* _d2Likelihoods_branch)[_rootPatternLinks[i]] - pow((* _dLikelihoods_branch)[_rootPatternLinks[i]], 2);
  const vector<unsigned int> * w = & _likelihoodData->getWeights();
  for(unsigned int i = 0; i < _nbDistinctSites; i++) d2 += (* w)[i] * ((* _d2Likelihoods_branch)[i] - pow((* _dLikelihoods_branch)[i], 2));
  return -d2;
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::resetLikelihoodArrays(const Node * node)
{
  for(unsigned int n = 0; n < node->getNumberOfSons(); n++)
  {
    const Node * subNode = node->getSon(n);
    resetLikelihoodArray(_likelihoodData->getLikelihoodArray(node, subNode));
  }
  if(node->hasFather())
  {
    const Node * father = node->getFather();
    resetLikelihoodArray(_likelihoodData->getLikelihoodArray(node, father));
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
  
  map<int, VVVdouble> * _likelihoods_node = & _likelihoodData->getLikelihoodArrays(node);
  unsigned int nbNodes = node->getNumberOfSons();
  for(unsigned int l = 0; l < nbNodes; l++)
  {
    //For each son node...  

    const Node * son = node->getSon(l);
    VVVdouble * _likelihoods_node_son = & (* _likelihoods_node)[son->getId()];
    
    if(son->isLeaf())
    {
      VVdouble * _likelihoods_leaf = & _likelihoodData->getLeafLikelihoods(son);
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
      map<int, VVVdouble> * _likelihoods_son = & _likelihoodData->getLikelihoodArrays(son);
      
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

void DRHomogeneousTreeLikelihood::computeSubtreeLikelihoodPrefix(const Node * node)
{
  if(! node->hasFather())
  { 
    // 'node' is the root of the tree.  
    // Just call the method on each son node:
    unsigned int nbSons = node->getNumberOfSons();
    for(unsigned int n = 0; n < nbSons; n++) computeSubtreeLikelihoodPrefix(node->getSon(n));
    return;
  }
  else
  {
    const Node * father = node->getFather();
    map<int, VVVdouble> * _likelihoods_node = & _likelihoodData->getLikelihoodArrays(node);
    map<int, VVVdouble> * _likelihoods_father = & _likelihoodData->getLikelihoodArrays(father);
    VVVdouble * _likelihoods_node_father = & (* _likelihoods_node)[father->getId()];
    if(node->isLeaf())
    {
      resetLikelihoodArray(*_likelihoods_node_father);
    }
  
    if(father->isLeaf())
    { 
      // If the tree is rooted by a leaf
      VVdouble * _likelihoods_leaf = & _likelihoodData->getLeafLikelihoods(father);
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
  
      unsigned int nbSons = nodes.size(); // In case of a bifurcating tree, this is equal to 1, excepted for the root.
      
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
        // Also take grand-father into account:
        const Node * fatherFather = father->getFather();
        tProb.push_back(& _pxy[father->getId()]); //!!! the difference here is that we use
                                                  //!!! _pxy[father] instead of _pxy[fatherFather].
        iLik.push_back(& (* _likelihoods_father)[fatherFather->getId()]);
        nbSons++;
      }
      computeLikelihoodFromArrays(iLik, tProb, *_likelihoods_node_father, nbSons, _nbDistinctSites, _nbClasses, _nbStates, false); 
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
    VVdouble * leavesLikelihoods_root = & _likelihoodData->getLeafLikelihoods(root);
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
  
  map<int, VVVdouble> * likelihoods_root = & _likelihoodData->getLikelihoodArrays(root);
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

  Vdouble f = _model->getFrequencies();
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
        (* rootLikelihoodsS_i_c) += f[x] * (* rootLikelihoods_i_c)[x];
        //cout << i << "\t" << c << "\t" << x << "\t" << f[x] << "\t" << (* rootLikelihoods_i_c)[x] << endl;
      }
      (* rootLikelihoodsSR)[i] += p[c] * (* rootLikelihoodsS_i)[c];
    }

    //Final checking (for numerical errors):
    if((* rootLikelihoodsSR)[i] < 0) (* rootLikelihoodsSR)[i] = 0.;
  }
}

/******************************************************************************/

VVVdouble DRHomogeneousTreeLikelihood::computeLikelihoodAtNode(const Node * node) const
{
  VVVdouble likelihoodArray(_nbDistinctSites);
  map<int, VVVdouble> * likelihoods_node = & _likelihoodData->getLikelihoodArrays(node);
  
  //Initialize likelihood array:
  if(node->isLeaf())
  {
    VVdouble * leavesLikelihoods_node = & _likelihoodData->getLeafLikelihoods(node);
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
    tProb.push_back(& _pxy[node->getId()]); // and not father!!!
    iLik.push_back(& (* likelihoods_node)[father->getId()]);
    nbNodes++;
  }
  computeLikelihoodFromArrays(iLik, tProb, likelihoodArray, nbNodes, _nbDistinctSites, _nbClasses, _nbStates, false);
  return likelihoodArray;
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
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::displayLikelihood(const Node * node)
{
  cout << "Likelihoods at node " << node->getId() << ": " << endl;
  for(unsigned int n = 0; n < node->getNumberOfSons(); n++) {
    const Node * subNode = node->getSon(n);
    cout << "Array for sub-node " << subNode->getId() << endl;
    displayLikelihoodArray(_likelihoodData->getLikelihoodArray(node, subNode));
  }
  if(node->hasFather())
  {
    const Node * father = node->getFather();
    cout << "Array for father node " << father->getId() << endl;
    displayLikelihoodArray(_likelihoodData->getLikelihoodArray(node, father));
  }
  cout << "                                         ***" << endl;
}

/*******************************************************************************/

void BranchLikelihood::computeLogLikelihood()
{
  _arrayTmp = *_array1;
  _lnL = 0;
  for(unsigned int i = 0; i < _array1->size(); i++)
  {
    VVdouble * arrayTmp_i = & _arrayTmp[i];
    const VVdouble * array2_i = & (*_array2)[i];
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      Vdouble * arrayTmp_i_c = & (*arrayTmp_i)[c];
      const Vdouble * array2_i_c = & (*array2_i)[c];
      VVdouble *_pxy_c = & _pxy[c];
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        Vdouble *_pxy_c_x = & (*_pxy_c)[x];
        double likelihood = 0;
        for(unsigned int y = 0; y < _nbStates; y++)
        {
          likelihood += (*_pxy_c_x)[y] * (*array2_i_c)[y];
        }
        (*arrayTmp_i_c)[x] *= likelihood;
      }
    }
  }

  for(unsigned int i = 0; i < _array1->size(); i++)
  {
    VVdouble * arrayTmp_i = & _arrayTmp[i];
    double Li = 0;
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      Vdouble * arrayTmp_i_c = & (*arrayTmp_i)[c];
      double rc = _rDist->getProbability(c);
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        Li += rc * _model->freq(x) * (* arrayTmp_i_c)[x];
      }
    }
    _lnL -= _weights[i] * log(Li);
  }
}

/*******************************************************************************/

void BranchLikelihood::computeAllTransitionProbabilities()
{
  double l = _parameters.getParameter("BrLen")->getValue(); 

  //Computes all pxy once for all:
  for(unsigned int c = 0; c < _nbClasses; c++)
  {
    VVdouble * _pxy_c = & _pxy[c];
    RowMatrix<double> Q = _model->getPij_t(l * _rDist->getCategory(c));
    for(unsigned int x = 0; x < _nbStates; x++)
    {
      Vdouble * _pxy_c_x = & (* _pxy_c)[x];
      for(unsigned int y = 0; y < _nbStates; y++)
      {
        (* _pxy_c_x)[y] = Q(x, y);
      }
    }
  } 
}

/*******************************************************************************/

double DRHomogeneousTreeLikelihood::testNNI(int nodeId) const throw (NodeException)
{
  const Node * son    = _tree->getNode(nodeId);
	if(!son->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::testNNI(). Node 'son' must not be the root node.", son);
  const Node * parent = son->getFather();
	if(!parent->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::testNNI(). Node 'parent' must not be the root node.", parent);
	const Node * grandFather = parent->getFather();
	//From here: Bifurcation assumed.
	//In case of multifurcation, an arbitrary uncle is chosen.
	//If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
	unsigned int parentPosition = grandFather->getSonPosition(*parent);
	//const Node * uncle = grandFather->getSon(parentPosition > 1 ? parentPosition - 1 : 1 - parentPosition);
	const Node * uncle = grandFather->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);
		
	//Retrieving arrays of interest:
	const DRASDRTreeLikelihoodNodeData * parentData = & _likelihoodData->getNodeData(parent);
	const VVVdouble                    * sonArray   = & parentData->getLikelihoodArrayForNeighbor(son);
	vector<const Node *> parentNeighbors = TreeTemplateTools::getRemainingNeighbors(parent, grandFather, son);
	unsigned int nbParentNeighbors = parentNeighbors.size();
	vector<const VVVdouble *> parentArrays(nbParentNeighbors);
	vector<const VVVdouble *> parentTProbs(nbParentNeighbors);
	for(unsigned int k = 0; k < nbParentNeighbors; k++)
  {
		const Node * n = parentNeighbors[k]; // This neighbor
		parentArrays[k] = & parentData->getLikelihoodArrayForNeighbor(n); 
    //if(n != grandFather) parentTProbs[k] = & _pxy[n->getId()];
    //else                 parentTProbs[k] = & _pxy[parent->getId()];
    parentTProbs[k] = & _pxy[n->getId()];
	}
	
	const DRASDRTreeLikelihoodNodeData * grandFatherData = & _likelihoodData->getNodeData(grandFather);
	const VVVdouble                    * uncleArray      = & grandFatherData->getLikelihoodArrayForNeighbor(uncle); 
	vector<const Node *> grandFatherNeighbors = TreeTemplateTools::getRemainingNeighbors(grandFather, parent, uncle);
	unsigned int nbGrandFatherNeighbors = grandFatherNeighbors.size();
	vector<const VVVdouble *> grandFatherArrays(nbGrandFatherNeighbors);
	vector<const VVVdouble *> grandFatherTProbs(nbGrandFatherNeighbors);
	for(unsigned int k = 0; k < nbGrandFatherNeighbors; k++)
  {
		const Node * n = grandFatherNeighbors[k]; // This neighbor
		grandFatherArrays[k] = & grandFatherData->getLikelihoodArrayForNeighbor(n); 
    if(grandFather->getFather() == NULL || n != grandFather->getFather())
    {
      grandFatherTProbs[k] = & _pxy[n->getId()];
    }
    else
    {
      grandFatherTProbs[k] = & _pxy[grandFather->getId()];
    }
	}

  //Compute array 1: grand father array
  VVVdouble array1 = *sonArray;
  resetLikelihoodArray(array1);
	grandFatherArrays.push_back(sonArray);
	grandFatherTProbs.push_back(& _pxy[son->getId()]);
  computeLikelihoodFromArrays(grandFatherArrays, grandFatherTProbs, array1, nbGrandFatherNeighbors + 1, _nbDistinctSites, _nbClasses, _nbStates, false); 
  
  //Compute array 2: parent array
  VVVdouble array2 = *uncleArray;
  resetLikelihoodArray(array2);
	parentArrays.push_back(uncleArray);
	parentTProbs.push_back(& _pxy[uncle->getId()]);
  computeLikelihoodFromArrays(parentArrays, parentTProbs, array2, nbParentNeighbors + 1, _nbDistinctSites, _nbClasses, _nbStates, false); 

  //Initialize BranchLikelihood:
  _brLikFunction->initModel(_model, _rateDistribution);
  _brLikFunction->initLikelihoods(&array1, &array2);
  ParameterList parameters;
  Parameter brLen = *_parameters.getParameter("BrLen" + TextTools::toString(parent->getId()));
  brLen.setName("BrLen");
  parameters.addParameter(brLen);
  _brLikFunction->setParameters(parameters);
  
  //Re-estimate branch length:
  _brentOptimizer->setFunction(_brLikFunction);
  _brentOptimizer->getStopCondition()->setTolerance(0.1);
  _brentOptimizer->setInitialInterval(brLen.getValue(), brLen.getValue()+0.01);
  _brentOptimizer->init(parameters);
  _brentOptimizer->optimize();
  //_brLenNNIValues[nodeId] = _brLikFunction->getParameterValue("BrLen");
  _brLenNNIValues[nodeId] = _brentOptimizer->getParameters().getParameter("BrLen")->getValue();
  _brLikFunction->resetLikelihoods(); //Array1 and Array2 will be destroyed after this function call.
                                      //We should not keep pointers towards them...

  //Return the resulting likelihood:
  return _brLikFunction->getValue() - getValue();
}
    
/*******************************************************************************/

void DRHomogeneousTreeLikelihood::doNNI(int nodeId) throw (NodeException)
{
  //Perform the topological move, the likelihood array will have to be recomputed...
  Node * son    = _tree->getNode(nodeId);
	if(!son->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::testNNI(). Node 'son' must not be the root node.", son);
  Node * parent = son->getFather();
	if(!parent->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::testNNI(). Node 'parent' must not be the root node.", parent);
	Node * grandFather = parent->getFather();
	//From here: Bifurcation assumed.
	//In case of multifurcation, an arbitrary uncle is chosen.
	//If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
	unsigned int parentPosition = grandFather->getSonPosition(*parent);
	Node * uncle = grandFather->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);
	//Swap nodes:
	parent->removeSon(*son);
	grandFather->removeSon(*uncle);
	parent->addSon(*uncle);
	grandFather->addSon(*son);
  string name = "BrLen" + TextTools::toString(parent->getId());
  if(_brLenNNIValues.find(nodeId) != _brLenNNIValues.end())
  {
    double length = _brLenNNIValues[nodeId];
    _brLenParameters.setParameterValue(name, length);
    _parameters.setParameterValue(name, length);
    parent->setDistanceToFather(length);
  }
  else cerr << "ERROR, branch not found: " << nodeId << endl;
  try { _brLenNNIParams.addParameter(*_brLenParameters.getParameter(name)); }
  catch(ParameterException & ex)
  {
    cerr << "DEBUG:" << endl;
    _brLenNNIParams.printParameters(cerr);
    cerr << "DEBUG:" << name << endl;
  }
  //In case of copy of this object, we must remove the constraint associated to this stored parameter:
  //(It should be also possible to update the pointer in the copy constructor,
  //but we do not need the constraint info here...).
  (*_brLenNNIParams.rbegin())->removeConstraint();
}

/*******************************************************************************/

