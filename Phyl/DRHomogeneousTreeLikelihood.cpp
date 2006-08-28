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
#include "SitePatterns.h"

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

void DRASDRTreeLikelihoodData::initLikelihoods(const SiteContainer & sites, const SubstitutionModel & model) throw (Exception)
{  
  if(sites.getNumberOfSequences() == 1) throw Exception("Error, only 1 sequence!");
  if(sites.getNumberOfSequences() == 0) throw Exception("Error, no sequence!");
  if(sites.getAlphabet()->getAlphabetType()
      != model.getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("DRASDRTreeLikelihoodData::initLikelihoods. Data and model must have the same alphabet type.",
        sites.getAlphabet(),
        model.getAlphabet());
  _alphabet = sites.getAlphabet();
  _nbStates = model.getNumberOfStates();
  _nbSites  = sites.getNumberOfSites();
  
  SitePatterns pattern(sites);
  if(_shrunkData != NULL) delete _shrunkData;
  _shrunkData = pattern.getSites();
  _rootWeights = pattern.getWeights();
  _rootPatternLinks = pattern.getIndices();
  _nbDistinctSites = _shrunkData->getNumberOfSites();
  
  //Init data:
  // Clone data for more efficiency on sequences access:
  const SiteContainer * sequences = new AlignedSequenceContainer(* _shrunkData);
  initLikelihoods(_tree->getRootNode(), * sequences, model);
  delete sequences;

  // Now initialize root likelihoods and derivatives:
  _rootLikelihoods.resize(_nbDistinctSites);
  _rootLikelihoodsS.resize(_nbDistinctSites);
  _rootLikelihoodsSR.resize(_nbDistinctSites);
  for(unsigned int i = 0; i < _nbDistinctSites; i++)
  {
    VVdouble * _rootLikelihoods_i = & _rootLikelihoods[i];
    Vdouble * _rootLikelihoodsS_i = & _rootLikelihoodsS[i];
    _rootLikelihoods_i->resize(_nbClasses);
    _rootLikelihoodsS_i->resize(_nbClasses);
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      Vdouble * _rootLikelihoods_i_c = & (* _rootLikelihoods_i)[c];
      _rootLikelihoods_i_c->resize(_nbStates);
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        (* _rootLikelihoods_i_c)[x] = 1.;
      }
    }
  }
}

/******************************************************************************/

void DRASDRTreeLikelihoodData::initLikelihoods(const Node * node, const SiteContainer & sites, const SubstitutionModel & model) throw (Exception)
{
  if(node->isLeaf())
  {
    // Init leaves likelihoods:
    const Sequence * seq;
    try
    {
      seq = sites.getSequence(node->getName());
    }
    catch (SequenceNotFoundException & snfe)
    {
      throw SequenceNotFoundException("DRASDRTreeLikelihoodData::initlikelihoods. Leaf name in tree not found in site container: ", (node->getName()));
    }
    VVdouble * leavesLikelihoods_leaf = & _leafData[node].getLikelihoodArray();
    leavesLikelihoods_leaf->resize(_nbDistinctSites);
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      Vdouble * leavesLikelihoods_leaf_i = & (* leavesLikelihoods_leaf)[i];
      leavesLikelihoods_leaf_i->resize(_nbStates);
      int state = seq->getValue(i);
      double test = 0.;
      for(unsigned int s = 0; s < _nbStates; s++)
      {
        //Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
        //otherwise value set to 0:
        ( * leavesLikelihoods_leaf_i)[s] = model.getInitValue(s, state);
        test += ( * leavesLikelihoods_leaf_i)[s];
      }
      if(test < 0.000001) cerr << "WARNING!!! Likelihood will be 0 for this site." << endl;
    }
  }

  //cout << "Node:\t" << node->getId() << endl;
  // We initialize each son node first:
  unsigned int nbSonNodes = node->getNumberOfSons();
  for(unsigned int l = 0; l < nbSonNodes; l++)
  {
    //For each son node,
    initLikelihoods(node->getSon(l), sites, model);
  }

  //Initialize likelihood vector:
  map<const Node *, VVVdouble> * _likelihoods_node = & _nodeData[node].getLikelihoodArrays();
  
  int nbSons = node->getNumberOfSons();
  
  for(int n = (node->hasFather() ? -1 : 0); n < nbSons; n++)
  {
    const Node * neighbor = (* node)[n];
    VVVdouble * _likelihoods_node_neighbor = & (* _likelihoods_node)[neighbor];
    
    _likelihoods_node_neighbor->resize(_nbDistinctSites);

    if(neighbor->isLeaf())
    {
      VVdouble * _leavesLikelihoods_leaf = & _leafData[neighbor].getLikelihoodArray();
      for(unsigned int i = 0; i < _nbDistinctSites; i++)
      {
        Vdouble  * _leavesLikelihoods_leaf_i = & (* _leavesLikelihoods_leaf)[i];
        VVdouble * _likelihoods_node_neighbor_i = & (* _likelihoods_node_neighbor)[i];
        _likelihoods_node_neighbor_i->resize(_nbClasses);
        for(unsigned int c = 0; c < _nbClasses; c++)
        {
          Vdouble * _likelihoods_node_neighbor_i_c = & (* _likelihoods_node_neighbor_i)[c];
          _likelihoods_node_neighbor_i_c->resize(_nbStates);
          for(unsigned int s = 0; s < _nbStates; s++)
          {
            (* _likelihoods_node_neighbor_i_c)[s] = (* _leavesLikelihoods_leaf_i)[s];
          }
        }
      }
    }
    else
    {
      for(unsigned int i = 0; i < _nbDistinctSites; i++)
      {
        VVdouble * _likelihoods_node_neighbor_i = & (* _likelihoods_node_neighbor)[i];
        _likelihoods_node_neighbor_i->resize(_nbClasses);
        for(unsigned int c = 0; c < _nbClasses; c++)
        {
          Vdouble * _likelihoods_node_neighbor_i_c = & (* _likelihoods_node_neighbor_i)[c];
          _likelihoods_node_neighbor_i_c->resize(_nbStates);
          for(unsigned int s = 0; s < _nbStates; s++)
          {
            (* _likelihoods_node_neighbor_i_c)[s] = 1.; //All likelihoods are initialized to 1.
          }
        }
      }
    }
  }

  // Initialize d and d2 likelihoods:
  Vdouble * _dLikelihoods_node = & _nodeData[node].getDLikelihoodArray();
  Vdouble * _d2Likelihoods_node = & _nodeData[node].getD2LikelihoodArray();
  _dLikelihoods_node->resize(_nbDistinctSites);
  _d2Likelihoods_node->resize(_nbDistinctSites);
    
}

/******************************************************************************/

DRHomogeneousTreeLikelihood::DRHomogeneousTreeLikelihood(
  TreeTemplate<Node> * tree,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose
)  throw (Exception):
  AbstractDiscreteRatesAcrossSitesTreeLikelihood(rDist, verbose), // We must do this since AbstractTreeLikelihood is virtual
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose)
{
  if(verbose) ApplicationTools::message << "Double-Recursive Homogeneous Tree Likelihood" << endl;  
  _likelihoodData = new DRASDRTreeLikelihoodData(*_tree, rDist->getNumberOfCategories());
    
  // Now initializes all parameters:
  initParameters();
  computeAllTransitionProbabilities();
  //fireParameterChanged(_parameters);
}

/******************************************************************************/

DRHomogeneousTreeLikelihood::DRHomogeneousTreeLikelihood(
  TreeTemplate<Node> * tree,
  const SiteContainer & data,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose
)  throw (Exception):
  AbstractDiscreteRatesAcrossSitesTreeLikelihood(rDist, verbose), // We must do this since AbstractTreeLikelihood is virtual
  AbstractHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose)
{
  if(verbose) ApplicationTools::message << "Double-Recursive Homogeneous Tree Likelihood" << endl;  
  _likelihoodData = new DRASDRTreeLikelihoodData(*_tree, rDist->getNumberOfCategories());
    
  setData(data);
  
  // Now initializes all parameters:
  initParameters();
  fireParameterChanged(_parameters);
}

/******************************************************************************/

DRHomogeneousTreeLikelihood::~DRHomogeneousTreeLikelihood()
{
  delete _likelihoodData;
  delete _data;
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::setData(const SiteContainer & sites) throw (Exception)
{
  _data = PatternTools::getSequenceSubset(sites, * _tree->getRootNode());
   if(_verbose) ApplicationTools::displayTask("Initializing data structure");
  _likelihoodData->initLikelihoods(*_data, *_model);
  if(_verbose) ApplicationTools::displayTaskDone();

  _nbSites = _likelihoodData->getNumberOfSites();
  _nbDistinctSites = _likelihoodData->getNumberOfDistinctSites();
  _nbStates = _likelihoodData->getNumberOfStates();
  
  if(_verbose) ApplicationTools::displayResult("Number of distinct sites",
      TextTools::toString(_nbDistinctSites));
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

  // For now we ignore the parameter that changed and we recompute all arrays...
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
  //double f = - geitLogLikelihood(); // For minimization.
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
  VVVdouble *  _pxy_node = &  _pxy[node];
  VVVdouble * _dpxy_node = & _dpxy[node];
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
  VVVdouble *   _pxy_node = &   _pxy[node];
  VVVdouble * _d2pxy_node = & _d2pxy[node];
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
  
  map<const Node *, VVVdouble> * _likelihoods_node = & _likelihoodData->getLikelihoodArrays(node);
  unsigned int nbNodes = node->getNumberOfSons();
  for(unsigned int l = 0; l < nbNodes; l++)
  {
    //For each son node...  

    const Node * son = node->getSon(l);
    VVVdouble * _likelihoods_node_son = & (* _likelihoods_node)[son];
    
    if(son->isLeaf())
    {
      VVdouble * _likelihoods_leaf = & _likelihoodData->getLeafLikelihoods(son);
      for(unsigned int i = 0; i < _nbDistinctSites; i++)
      {
        //For each site in the sequence,
        Vdouble * _likelihoods_leaf_i     = & (* _likelihoods_leaf)[i];
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
      map<const Node *, VVVdouble> * _likelihoods_son = & _likelihoodData->getLikelihoodArrays(son);

      for(unsigned int n = 0; n < nbSons; n++)
      {
        const Node * sonson = son->getSon(n);
        VVVdouble * _pxy_sonson = & _pxy[sonson];
        VVVdouble * _likelihoods_son_son = & (* _likelihoods_son)[sonson];

        for(unsigned int i = 0; i < _nbDistinctSites; i++)
        {
          //For each site in the sequence,
          VVdouble * _likelihoods_son_son_i  = & (* _likelihoods_son_son)[i];
          VVdouble * _likelihoods_node_son_i = & (* _likelihoods_node_son)[i];
          for(unsigned int c = 0; c < _nbClasses; c++)
          {
            //For each rate classe,
            Vdouble * _likelihoods_son_son_i_c  = & (* _likelihoods_son_son_i)[c];
            Vdouble * _likelihoods_node_son_i_c = & (* _likelihoods_node_son_i)[c];
            VVdouble * _pxy_sonson_c = & (* _pxy_sonson)[c];
            for(unsigned int x = 0; x < _nbStates; x++)
            {
              //For each initial state,
              Vdouble * _pxy_sonson_c_x = & (* _pxy_sonson_c)[x];
              double likelihood = 0;
              for(unsigned int y = 0; y < _nbStates; y++)
              {
                likelihood += (* _pxy_sonson_c_x)[y] * (* _likelihoods_son_son_i_c)[y];
              }
              // We store this conditionnal likelihood into the corresponding array:
              (* _likelihoods_node_son_i_c)[x] *= likelihood;
            }
          }
        }
      }
    }
  }
}

/******************************************************************************/

void DRHomogeneousTreeLikelihood::computeSubtreeLikelihoodPrefix(const Node * node)
{
  //if(node->isLeaf())
  //{
  //  return;
  //}
  if(! node->hasFather())
  { // 'node' is the root of the tree.  
    // Just call the method on each son node:
    unsigned int nbSons = node->getNumberOfSons();
    for(unsigned int n = 0; n < nbSons; n++) computeSubtreeLikelihoodPrefix(node->getSon(n));
    return;
  }
  else
  {
    const Node * father = node->getFather();
    map<const Node *, VVVdouble> * _likelihoods_node = & _likelihoodData->getLikelihoodArrays(node);
    map<const Node *, VVVdouble> * _likelihoods_father = & _likelihoodData->getLikelihoodArrays(father);
    VVVdouble * _likelihoods_node_father = & (* _likelihoods_node)[father];
  
    if(father->isLeaf())
    { // If the tree is rooted by a leaf
      VVdouble * _likelihoods_leaf = & _likelihoodData->getLeafLikelihoods(father);
      for(unsigned int i = 0; i < _nbDistinctSites; i++) {
        //For each site in the sequence,
        Vdouble * _likelihoods_leaf_i     = & (* _likelihoods_leaf)[i];
        VVdouble * _likelihoods_node_father_i = & (* _likelihoods_node_father)[i];
        for(unsigned int c = 0; c < _nbClasses; c++) {
          //For each rate classe,
          Vdouble * _likelihoods_node_father_i_c = & (* _likelihoods_node_father_i)[c];
          for(unsigned int x = 0; x < _nbStates; x++) {
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
      for(unsigned int n = 0; n < nbSons; n++)
      {
        const Node * fatherSon = nodes[n];
        VVVdouble * _pxy_fatherSon = & _pxy[fatherSon];
        VVVdouble * _likelihoods_father_son = & (* _likelihoods_father)[fatherSon];

        for(unsigned int i = 0; i < _nbDistinctSites; i++)
        {
          //For each site in the sequence,
          VVdouble * _likelihoods_father_son_i  = & (* _likelihoods_father_son)[i];
          VVdouble * _likelihoods_node_father_i = & (* _likelihoods_node_father)[i];
          for(unsigned int c = 0; c < _nbClasses; c++)
          {
            //For each rate classe,
            Vdouble * _likelihoods_father_son_i_c  = & (* _likelihoods_father_son_i)[c];
            Vdouble * _likelihoods_node_father_i_c = & (* _likelihoods_node_father_i)[c];
            VVdouble * _pxy_fatherSon_c = & (* _pxy_fatherSon)[c];
            for(unsigned int x = 0; x < _nbStates; x++)
            {
              //For each initial state,
              Vdouble * _pxy_fatherSon_c_x = & (* _pxy_fatherSon_c)[x];
              double likelihood = 0;
              for(unsigned int y = 0; y < _nbStates; y++)
              {
                likelihood += (* _pxy_fatherSon_c_x)[y] * (* _likelihoods_father_son_i_c)[y];
              }
              // We store this conditionnal likelihood into the corresponding array:
              (* _likelihoods_node_father_i_c)[x] *= likelihood;
            }
          }
        }
      }
    
      if(father->hasFather())
      {
        // Also take grand-father into account:
        const Node * fatherFather = father->getFather();
        VVVdouble * _pxy_fatherFather = & _pxy[father]; //!!! the difference here is that we use
                                                         //!!! _pxy[father] instead of _pxy[fatherFather].
        VVVdouble * _likelihoods_father_father = & (* _likelihoods_father)[fatherFather];

        for(unsigned int i = 0; i < _nbDistinctSites; i++)
        {
          //For each site in the sequence,
          VVdouble * _likelihoods_father_father_i = & (* _likelihoods_father_father)[i];
          VVdouble * _likelihoods_node_father_i   = & (* _likelihoods_node_father)[i];
          for(unsigned int c = 0; c < _nbClasses; c++)
          {
            //For each rate classe,
            Vdouble * _likelihoods_father_father_i_c = & (* _likelihoods_father_father_i)[c];
            Vdouble * _likelihoods_node_father_i_c   = & (* _likelihoods_node_father_i)[c];
            VVdouble * _pxy_fatherFather_c = & (* _pxy_fatherFather)[c];
            for(unsigned int x = 0; x < _nbStates; x++)
            {
              //For each initial state,
              Vdouble * _pxy_fatherFather_c_x = & (* _pxy_fatherFather_c)[x];
              double likelihood = 0;
              for(unsigned int y = 0; y < _nbStates; y++)
              {
                likelihood += (* _pxy_fatherFather_c_x)[y] * (* _likelihoods_father_father_i_c)[y];
              }
              // We store this conditionnal likelihood into the corresponding array:
              (* _likelihoods_node_father_i_c)[x] *= likelihood;
            }
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
  
  map<const Node *, VVVdouble> * likelihoods_root = & _likelihoodData->getLikelihoodArrays(root);
  unsigned int nbNodes = root->getNumberOfSons();
  for(unsigned int n = 0; n < nbNodes; n++)
  {
    const Node * son = root->getSon(n);
    VVVdouble * _pxy_son = & _pxy[son];
    VVVdouble * likelihoods_root_son = & (* likelihoods_root)[son];

    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      //For each site in the sequence,
      VVdouble * likelihoods_root_son_i  = & (* likelihoods_root_son)[i];
      VVdouble * rootLikelihoods_i = & (* rootLikelihoods)[i];
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        //For each rate classe,
        Vdouble * likelihoods_root_son_i_c  = & (* likelihoods_root_son_i)[c];
        Vdouble * rootLikelihoods_i_c = & (* rootLikelihoods_i)[c];
        VVdouble * _pxy_son_c = & (* _pxy_son)[c];
        for(unsigned int x = 0; x < _nbStates; x++)
        {
          //For each initial state,
          Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
          double likelihood = 0;
          for(unsigned int y = 0; y < _nbStates; y++)
          {
            likelihood += (* _pxy_son_c_x)[y] * (* likelihoods_root_son_i_c)[y];
            //cout << i << "\t" << c << "\t" << x << "\t" << y << "\t" <<  (* _pxy_son_c_x)[y] << "\t" << (* likelihoods_root_son_i_c)[y] << endl;
          }
          // We store this conditionnal likelihood into the corresponding array:
          (* rootLikelihoods_i_c)[x] *= likelihood;
        }
      }
    }
  }

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
  map<const Node *, VVVdouble> * likelihoods_node = & _likelihoodData->getLikelihoodArrays(node);
  
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
        for(unsigned int s = 0; s < _nbStates; s++)
        {
          (* likelihoodArray_i_c)[s] = 1.;
        }
      }
    }
  }
  
  unsigned int nbNodes = node->getNumberOfSons();
  for(unsigned int n = 0; n < nbNodes; n++)
  {
    const Node * son = node->getSon(n);
    VVVdouble * _pxy_son = & _pxy[son];
    VVVdouble * likelihoods_node_son = & (* likelihoods_node)[son];

    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      //For each site in the sequence,
      VVdouble * likelihoods_node_son_i  = & (* likelihoods_node_son)[i];
      VVdouble * likelihoodArray_i = & likelihoodArray[i];
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        //For each rate classe,
        Vdouble * likelihoods_node_son_i_c  = & (* likelihoods_node_son_i)[c];
        Vdouble * likelihoodArray_i_c = & (* likelihoodArray_i)[c];
        VVdouble * _pxy_son_c = & (* _pxy_son)[c];
        for(unsigned int x = 0; x < _nbStates; x++)
        {
          //For each initial state,
          Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
          double likelihood = 0;
          for(unsigned int y = 0; y < _nbStates; y++)
          {
            likelihood += (* _pxy_son_c_x)[y] * (* likelihoods_node_son_i_c)[y];
          }
          // We store this conditionnal likelihood into the corresponding array:
          (* likelihoodArray_i_c)[x] *= likelihood;
        }
      }
    }
  }
  
  if(node->hasFather())
  {
    const Node * son = node->getFather();
    VVVdouble * _pxy_son = & _pxy[node]; // and not son!!!
    VVVdouble * likelihoods_node_son = & (* likelihoods_node)[son];

    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      //For each site in the sequence,
      VVdouble * likelihoods_node_son_i  = & (* likelihoods_node_son)[i];
      VVdouble * likelihoodArray_i = & likelihoodArray[i];
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        //For each rate classe,
        Vdouble * likelihoods_node_son_i_c  = & (* likelihoods_node_son_i)[c];
        Vdouble * likelihoodArray_i_c = & (* likelihoodArray_i)[c];
        VVdouble * _pxy_son_c = & (* _pxy_son)[c];
        for(unsigned int x = 0; x < _nbStates; x++)
        {
          //For each initial state,
          Vdouble * _pxy_son_c_x = & (* _pxy_son_c)[x];
          double likelihood = 0;
          for(unsigned int y = 0; y < _nbStates; y++)
          {
            likelihood += (* _pxy_son_c_x)[y] * (* likelihoods_node_son_i_c)[y];
          }
          // We store this conditionnal likelihood into the corresponding array:
          (* likelihoodArray_i_c)[x] *= likelihood;
        }
      }
    }
  }
  return likelihoodArray;
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

