//
// File: ClockTreeLikelihood.cpp
// Created by: Benoît Nabholz
// Created on: Fri Apr 06 14:11 2007
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

#include "ClockTreeLikelihood.h"
#include "TreeTemplateTools.h"
#include <iostream>
using namespace std;

/******************************************************************************/

IncludingInterval ClockTreeLikelihood::PERCENT_CONSTRAINT(0.0,1.0);

/******************************************************************************/

ClockTreeLikelihood::ClockTreeLikelihood(
  const Tree & tree,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose)
throw (Exception):
  DRHomogeneousTreeLikelihood(tree, model, rDist, false, verbose)
{
  init();
}

/******************************************************************************/

ClockTreeLikelihood::ClockTreeLikelihood(
  const Tree & tree,
  const SiteContainer & data,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose)
throw (Exception):
  DRHomogeneousTreeLikelihood(tree, data, model, rDist, false, verbose)
{
  init();
}

/******************************************************************************/

void ClockTreeLikelihood::init()
{
  //Check is the tree is rooted:
  if(!_tree->isRooted()) throw Exception("ClockTreeLikelihood::init(). Tree is unrooted!");
  if(TreeTemplateTools::isMultifurcating(*_tree->getRootNode())) throw Exception("ClockTreeLikelihood::init(). Tree is multifurcating.");
  TreeTools::convertToClockTree(*_tree, _tree->getRootId(), true);
}

/******************************************************************************/

void ClockTreeLikelihood::applyParameters() throw (Exception)
{
  if(!_initialized) throw Exception("ClockTreeLikelihood::applyParameters(). Object not initialized.");
   //Apply branch lengths:
  TreeTemplateTools::getHeights(* _tree->getRootNode(), _previousHeights);
  _previousTotalHeight = _previousHeights[_tree->getRootNode()];
  computeBranchLengthsFromHeights(_tree->getRootNode());
  //Apply substitution model parameters:
  _model->matchParametersValues(_parameters);
  //Apply rate distribution parameters:
  _rateDistribution->matchParametersValues(_parameters);
}

/******************************************************************************/

void ClockTreeLikelihood::initBranchLengthsParameters()
{
  //Check branch lengths first:
  for(unsigned int i = 0; i < _nbNodes; i++)
  {
    double d = _minimumBrLen;
    if(!_nodes[i]->hasDistanceToFather())
    {
      cout << "WARNING!!! Missing branch length " << i << ". Value is set to " << _minimumBrLen << endl;
      _nodes[i]->setDistanceToFather(_minimumBrLen);
    }
    else
    {
      d = _nodes[i]->getDistanceToFather();
      if (d < _minimumBrLen)
      {
        cout << "WARNING!!! Branch length " << i << " is too small: " << d << ". Value is set to " << _minimumBrLen << endl;
        _nodes[i]->setDistanceToFather(_minimumBrLen);
        d = _minimumBrLen;
      }
    }
  }

  _brLenParameters.reset();
  _totalHeightParameter.reset();
  map<const Node *, double> heights;
  TreeTemplateTools::getHeights(*_tree->getRootNode(), heights);
  double totalHeight = heights[_tree->getRootNode()];
  _totalHeightParameter.addParameter(Parameter("TotalHeight", totalHeight, _brLenConstraint)); 
   for(map<const Node *, double>::iterator it = heights.begin(); it != heights.end(); it++)
  {
    if(!it->first->isLeaf() && it->first->hasFather())
    {
      _brLenParameters.addParameter(Parameter("HeightP" + TextTools::toString(it->first->getId()), it->second / totalHeight, &PERCENT_CONSTRAINT));
    }
  }
}

/******************************************************************************/

double ClockTreeLikelihood::computeBranchLengthsFromHeights(Node * node) throw (Exception)
{
  double minimumBrLenP = _minimumBrLen / _previousTotalHeight;
  if(node->isLeaf()) return minimumBrLenP;
  double totalHeight = _parameters.getParameter("TotalHeight")->getValue();
  //Compute bounds:
  double previousFatherHeightP = 1.;
  if(node->hasFather())
  {
    Node * father = node->getFather();
    previousFatherHeightP = _previousHeights[father] / _previousTotalHeight;
  }
  double nodeHeightP = 1.;

  //This parameter:
  if(node->hasFather())
  {
    Parameter * p = _parameters.getParameter(string("HeightP") + TextTools::toString(node->getId()));
    if(p == NULL) throw Exception("ClockTreeLikelihood::computeBranchLengthsFromHeights(). Parameter " + string("HeightP") + TextTools::toString(node->getId()) + " was not found."); 
    nodeHeightP = p->getValue();

    if(nodeHeightP > previousFatherHeightP - minimumBrLenP)
    {
      p->setValue(previousFatherHeightP - minimumBrLenP);
      //cout << "TOO LARGE: " << node->getId() << "\t" << nodeHeight << "\t" << previousFatherHeight << endl;
      nodeHeightP = previousFatherHeightP - minimumBrLenP;
    }
    double maxSonHeightP = minimumBrLenP;
    for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
    {
      Node * son = node->getSon(i);
      if(!son->isLeaf())
      {
        double previousSonHeightP = _previousHeights[son] / _previousTotalHeight;
        if(previousSonHeightP > maxSonHeightP) maxSonHeightP = previousSonHeightP;
      }
    }
    if(nodeHeightP < maxSonHeightP + minimumBrLenP)
    {
      p->setValue(maxSonHeightP + minimumBrLenP);
      //cout << "TOO SMALL: " << node->getId() << "\t" << nodeHeight << "\t" << maxSonHeight << endl;
      nodeHeightP = maxSonHeightP + minimumBrLenP;
    }
  }

  //Recursive call:
  for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    Node * son = node->getSon(i);
    double sonHeightP = computeBranchLengthsFromHeights(son);
    son->setDistanceToFather(std::max(_minimumBrLen, (nodeHeightP - sonHeightP) * totalHeight));
  }
  return nodeHeightP;
}

/******************************************************************************/

void ClockTreeLikelihood::fireParameterChanged(const ParameterList & params)
{
  applyParameters();

  computeAllTransitionProbabilities();

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

ParameterList ClockTreeLikelihood::getNonDerivableParameters() const throw (Exception)
{
  if(!_initialized) throw Exception("ClockTreeLikelihood::getNonDerivableParameters(). Object is not initialized.");
  ParameterList tmp = DRHomogeneousTreeLikelihood::getNonDerivableParameters();
  tmp.addParameters(getTotalHeightParameter());
  return tmp;
}

/******************************************************************************/

ParameterList ClockTreeLikelihood::getTotalHeightParameter() const throw (Exception)
{
  if(!_initialized) throw Exception("ClockTreeLikelihood::getTotalHeightParameter(). Object is not initialized.");
  return _parameters.getCommonParametersWith(_totalHeightParameter);
}

/******************************************************************************/

void ClockTreeLikelihood::initParameters()
{
  // Reset parameters:
  _parameters.reset();
  
  // Branch lengths:
  initBranchLengthsParameters();
  _parameters.addParameters(_brLenParameters);

  // Total height:
  _parameters.addParameters(_totalHeightParameter);
  
  // Substitution model:
  _parameters.addParameters(_model->getParameters());
  
  // Rate distribution:
  _parameters.addParameters(_rateDistribution->getParameters());
}

/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/  

void ClockTreeLikelihood::computeTreeDLikelihoodAtNode(const Node * node)
{
  if(node->isLeaf()) return; //The height of a leaf is always 0!
  double totalHeight  = _parameters.getParameter("TotalHeight")->getValue();
  const Node * father = node->getFather();
  const Node * son1   = node->getSon(0);
  const Node * son2   = node->getSon(1);
  Vdouble * _dLikelihoods_node = & _likelihoodData->getDLikelihoodArray(node);
  //Branch length toward father:
  VVVdouble * _likelihoods_node_father = & _likelihoodData->getLikelihoodArray(node, father);
  VVVdouble * _pxy_node0 = & _pxy[node->getId()];
  VVVdouble * _dpxy_node0 = & _dpxy[node->getId()];
  //Branch length toward son 1:
  VVVdouble * _likelihoods_node_son1   = & _likelihoodData->getLikelihoodArray(node, son1);
  VVVdouble * _pxy_node1 = & _pxy[son1->getId()];
  VVVdouble * _dpxy_node1 = & _dpxy[son1->getId()];
  //Branch length toward son 2:
  VVVdouble * _likelihoods_node_son2   = & _likelihoodData->getLikelihoodArray(node, son2);
  VVVdouble * _pxy_node2 = & _pxy[son2->getId()];
  VVVdouble * _dpxy_node2 = & _dpxy[son2->getId()];
  Vdouble * rootLikelihoodsSR = & _likelihoodData->getRootRateSiteLikelihoodArray();
 
  double dLi, dLic, dLicx;
  double l0, dl0; //Father
  double l1, dl1; //Son 1
  double l2, dl2; //Son 2

  for(unsigned int i = 0; i < _nbDistinctSites; i++)
  {
    VVdouble * _likelihoods_node_father_i = & (* _likelihoods_node_father)[i];
    VVdouble * _likelihoods_node_son1_i = & (* _likelihoods_node_son1)[i];
    VVdouble * _likelihoods_node_son2_i = & (* _likelihoods_node_son2)[i];
    dLi = 0;
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      Vdouble * _likelihoods_node_father_i_c = & (* _likelihoods_node_father_i)[c];
      Vdouble * _likelihoods_node_son1_i_c = & (* _likelihoods_node_son1_i)[c];
      Vdouble * _likelihoods_node_son2_i_c = & (* _likelihoods_node_son2_i)[c];
      VVdouble * _pxy_node0_c = & (* _pxy_node0)[c];
      VVdouble * _pxy_node1_c = & (* _pxy_node1)[c];
      VVdouble * _pxy_node2_c = & (* _pxy_node2)[c];
      VVdouble * _dpxy_node0_c = & (* _dpxy_node0)[c];
      VVdouble * _dpxy_node1_c = & (* _dpxy_node1)[c];
      VVdouble * _dpxy_node2_c = & (* _dpxy_node2)[c];
      dLic = 0;
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        l0 = 0; dl0 = 0;
        l1 = 0; dl1 = 0;
        l2 = 0; dl2 = 0;
        Vdouble * _pxy_node0_c_x = & (* _pxy_node0_c)[x];
        Vdouble * _pxy_node1_c_x = & (* _pxy_node1_c)[x];
        Vdouble * _pxy_node2_c_x = & (* _pxy_node2_c)[x];
        Vdouble * _dpxy_node0_c_x = & (* _dpxy_node0_c)[x];
        Vdouble * _dpxy_node1_c_x = & (* _dpxy_node1_c)[x];
        Vdouble * _dpxy_node2_c_x = & (* _dpxy_node2_c)[x];
        dLicx = 0;
        for(unsigned int y = 0; y < _nbStates; y++)
        {
          l0 += (* _pxy_node0_c_x)[y] * (* _likelihoods_node_father_i_c)[y];
          l1 += (* _pxy_node1_c_x)[y] * (* _likelihoods_node_son1_i_c)[y];
          l2 += (* _pxy_node2_c_x)[y] * (* _likelihoods_node_son2_i_c)[y];
          dl0 += - totalHeight * (* _dpxy_node0_c_x)[y] * (* _likelihoods_node_father_i_c)[y];
          dl1 += totalHeight * (* _dpxy_node1_c_x)[y] * (* _likelihoods_node_son1_i_c)[y];
          dl2 += totalHeight * (* _dpxy_node2_c_x)[y] * (* _likelihoods_node_son2_i_c)[y];
        }
        dLicx = dl0 * l1 * l2 + l0 * dl1 * l2 + l0 * l1 * dl2;
        dLic += _model->freq(x) * dLicx;  
      }
      dLi += _rateDistribution->getProbability(c) * dLic;
    }
    (* _dLikelihoods_node)[i] = dLi / (* rootLikelihoodsSR)[i]; 
  }
}

/******************************************************************************/

void ClockTreeLikelihood::computeTreeDLikelihoods()
{
  for(unsigned int k = 0; k < _nbNodes; k++)
  {
    computeTreeDLikelihoodAtNode(_nodes[k]);
  }
}

/******************************************************************************/

double ClockTreeLikelihood::getFirstOrderDerivative(const string & variable) const
throw (Exception)
{ 
  Parameter * p = _parameters.getParameter(variable);
  if(p == NULL) throw ParameterNotFoundException("ClockTreeLikelihood::getFirstOrderDerivative", variable);
  if(getRateDistributionParameters().getParameter(variable) != NULL)
  {
    throw Exception("Derivatives respective to rate distribution parameters are not implemented.");
  }
  if(getSubstitutionModelParameters().getParameter(variable) != NULL)
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }
  if(variable == "TotalHeight")
  {
    throw Exception("Derivative respective to 'TotalHeight' not implemented.");
  }
  
  //
  // Computation for branch lengths:
  //
  
  // Get the node with the branch whose length must be derivated:
  int brI = TextTools::toInt(variable.substr(7));
  Node * branch = _nodes[brI];
  Vdouble * _dLikelihoods_branch = & _likelihoodData->getDLikelihoodArray(branch);
  double d = 0;
  const vector<unsigned int> * w = & _likelihoodData->getWeights();
  for(unsigned int i = 0; i < _nbDistinctSites; i++)
    d += (* w)[i] * (* _dLikelihoods_branch)[i];
  return -d;
}

/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/  

void ClockTreeLikelihood::computeTreeD2LikelihoodAtNode(const Node * node)
{
  if(node->isLeaf()) return; //The height of a leaf is always 0!
  double totalHeight  = _parameters.getParameter("TotalHeight")->getValue();
  const Node * father = node->getFather();
  const Node * son1   = node->getSon(0);
  const Node * son2   = node->getSon(1);
  Vdouble * _d2Likelihoods_node = & _likelihoodData->getD2LikelihoodArray(node);
  //Branch length toward father:
  VVVdouble * _likelihoods_node_father = & _likelihoodData->getLikelihoodArray(node, father);
  VVVdouble * _pxy_node0 = & _pxy[node->getId()];
  VVVdouble * _dpxy_node0 = & _dpxy[node->getId()];
  VVVdouble * _d2pxy_node0 = & _d2pxy[node->getId()];
  //Branch length toward son 1:
  VVVdouble * _likelihoods_node_son1   = & _likelihoodData->getLikelihoodArray(node, son1);
  VVVdouble * _pxy_node1 = & _pxy[son1->getId()];
  VVVdouble * _dpxy_node1 = & _dpxy[son1->getId()];
  VVVdouble * _d2pxy_node1 = & _d2pxy[son1->getId()];
  //Branch length toward son 2:
  VVVdouble * _likelihoods_node_son2   = & _likelihoodData->getLikelihoodArray(node, son2);
  VVVdouble * _pxy_node2 = & _pxy[son2->getId()];
  VVVdouble * _dpxy_node2 = & _dpxy[son2->getId()];
  VVVdouble * _d2pxy_node2 = & _d2pxy[son2->getId()];
  Vdouble * rootLikelihoodsSR = & _likelihoodData->getRootRateSiteLikelihoodArray();
 
  double d2Li, d2Lic, d2Licx;
  double l0, dl0, d2l0; //Father
  double l1, dl1, d2l1; //Son 1
  double l2, dl2, d2l2; //Son 2

  for(unsigned int i = 0; i < _nbDistinctSites; i++)
  {
    VVdouble * _likelihoods_node_father_i = & (* _likelihoods_node_father)[i];
    VVdouble * _likelihoods_node_son1_i = & (* _likelihoods_node_son1)[i];
    VVdouble * _likelihoods_node_son2_i = & (* _likelihoods_node_son2)[i];
    d2Li = 0;
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      Vdouble * _likelihoods_node_father_i_c = & (* _likelihoods_node_father_i)[c];
      Vdouble * _likelihoods_node_son1_i_c = & (* _likelihoods_node_son1_i)[c];
      Vdouble * _likelihoods_node_son2_i_c = & (* _likelihoods_node_son2_i)[c];
      VVdouble * _pxy_node0_c = & (* _pxy_node0)[c];
      VVdouble * _pxy_node1_c = & (* _pxy_node1)[c];
      VVdouble * _pxy_node2_c = & (* _pxy_node2)[c];
      VVdouble * _dpxy_node0_c = & (* _dpxy_node0)[c];
      VVdouble * _dpxy_node1_c = & (* _dpxy_node1)[c];
      VVdouble * _dpxy_node2_c = & (* _dpxy_node2)[c];
      VVdouble * _d2pxy_node0_c = & (* _d2pxy_node0)[c];
      VVdouble * _d2pxy_node1_c = & (* _d2pxy_node1)[c];
      VVdouble * _d2pxy_node2_c = & (* _d2pxy_node2)[c];
      d2Lic = 0;
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        l0 = 0; dl0 = 0; d2l0 = 0;
        l1 = 0; dl1 = 0; d2l1 = 0;
        l2 = 0; dl2 = 0; d2l2 = 0;
        Vdouble * _pxy_node0_c_x = & (* _pxy_node0_c)[x];
        Vdouble * _pxy_node1_c_x = & (* _pxy_node1_c)[x];
        Vdouble * _pxy_node2_c_x = & (* _pxy_node2_c)[x];
        Vdouble * _dpxy_node0_c_x = & (* _dpxy_node0_c)[x];
        Vdouble * _dpxy_node1_c_x = & (* _dpxy_node1_c)[x];
        Vdouble * _dpxy_node2_c_x = & (* _dpxy_node2_c)[x];
        Vdouble * _d2pxy_node0_c_x = & (* _d2pxy_node0_c)[x];
        Vdouble * _d2pxy_node1_c_x = & (* _d2pxy_node1_c)[x];
        Vdouble * _d2pxy_node2_c_x = & (* _d2pxy_node2_c)[x];
        d2Licx = 0;
        for(unsigned int y = 0; y < _nbStates; y++)
        {
          l0 += (* _pxy_node0_c_x)[y] * (* _likelihoods_node_father_i_c)[y];
          l1 += (* _pxy_node1_c_x)[y] * (* _likelihoods_node_son1_i_c)[y];
          l2 += (* _pxy_node2_c_x)[y] * (* _likelihoods_node_son2_i_c)[y];
          dl0 += - totalHeight * (* _dpxy_node0_c_x)[y] * (* _likelihoods_node_father_i_c)[y];
          dl1 += totalHeight * (* _dpxy_node1_c_x)[y] * (* _likelihoods_node_son1_i_c)[y];
          dl2 += totalHeight * (* _dpxy_node2_c_x)[y] * (* _likelihoods_node_son2_i_c)[y];
          d2l0 += totalHeight * totalHeight * (* _d2pxy_node0_c_x)[y] * (* _likelihoods_node_father_i_c)[y];
          d2l1 += totalHeight * totalHeight * (* _d2pxy_node1_c_x)[y] * (* _likelihoods_node_son1_i_c)[y];
          d2l2 += totalHeight * totalHeight * (* _d2pxy_node2_c_x)[y] * (* _likelihoods_node_son2_i_c)[y];
        }
        d2Licx = d2l0 * l1 * l2 
               + l0 * d2l1 * l2
               + l0 * l1 * d2l2
               + 2 * (dl0 * dl1 * l2
                    + dl0 * l1 * dl2
                    + l0 * dl1 * dl2);
        d2Lic += _model->freq(x) * d2Licx;  
      }
      d2Li += _rateDistribution->getProbability(c) * d2Lic;
    }
    (* _d2Likelihoods_node)[i] = d2Li / (* rootLikelihoodsSR)[i]; 
  }
}

/******************************************************************************/

void ClockTreeLikelihood::computeTreeD2Likelihoods()
{
  for(unsigned int k = 0; k < _nbNodes; k++)
  {
    computeTreeD2LikelihoodAtNode(_nodes[k]);
  }
}

/******************************************************************************/

double ClockTreeLikelihood::getSecondOrderDerivative(const string & variable) const
throw (Exception)
{
  Parameter * p = _parameters.getParameter(variable);
  if(p == NULL) throw ParameterNotFoundException("ClockTreeLikelihood::getSecondOrderDerivative", variable);
  if(getRateDistributionParameters().getParameter(variable) != NULL)
  {
    throw Exception("Derivatives respective to rate distribution parameters are not implemented.");
  }
  if(getSubstitutionModelParameters().getParameter(variable) != NULL)
  {
    throw Exception("Derivatives respective to substitution model parameters are not implemented.");
  }
  if(variable == "TotalHeight")
  {
    throw Exception("Derivative respective to 'TotalHeight' not implemented.");
  }
  
  //
  // Computation for branch lengths:
  //
  
  // Get the node with the branch whose length must be derivated:
  int brI = TextTools::toInt(variable.substr(7));
  Node * branch = _nodes[brI];
  Vdouble * _dLikelihoods_branch = & _likelihoodData->getDLikelihoodArray(branch);
  Vdouble * _d2Likelihoods_branch = & _likelihoodData->getD2LikelihoodArray(branch);
  double d2 = 0;
  const vector<unsigned int> * w = & _likelihoodData->getWeights();
  for(unsigned int i = 0; i < _nbDistinctSites; i++) d2 += (* w)[i] * ((* _d2Likelihoods_branch)[i] - pow((* _dLikelihoods_branch)[i], 2));
  return -d2;
}

/******************************************************************************/

