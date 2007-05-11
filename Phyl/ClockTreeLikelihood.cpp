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

ClockTreeLikelihood::~ClockTreeLikelihood()
{
  //Delete previous contraints:
  for(unsigned int i = 0; i < _heightConstraints.size(); i++) delete _heightConstraints[i];
}

/******************************************************************************/

void ClockTreeLikelihood::init()
{
  //Check is ithe tree is rooted:
  if(!_tree->isRooted()) throw Exception("ClockTreeLikelihood::init(). Tree is unrooted!");
  if(TreeTemplateTools::isMultifurcating(*_tree->getRootNode())) throw Exception("ClockTreeLikelihood::init(). Tree is multifurcating.");
}

/******************************************************************************/

void ClockTreeLikelihood::applyParameters() throw (Exception)
{
  if(!_initialized) throw Exception("ClockTreeLikelihood::applyParameters(). Object not initialized.");
   //Apply branch lengths:
  _conflictingParameters.reset();
  _brLenParameters.matchParametersValues(_parameters);
  _totalHeightParameter.matchParametersValues(_parameters);
  computeBranchLengthsFromHeights(_tree->getRootNode());
  //Update parameters in cas of global constraints:
  _parameters.matchParametersValues(_brLenParameters);

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
      }
    }
  }

  _brLenParameters.reset();
  _totalHeightParameter.reset();

  //Delete previous contraints:
  for(unsigned int i = 0; i < _heightConstraints.size(); i++) delete _heightConstraints[i];
  _heightConstraints.clear();

  map<const Node *, double> heights;
  TreeTemplateTools::getHeights(*_tree->getRootNode(), heights);
  double totalHeight = heights[_tree->getRootNode()];
  _totalHeightParameter.addParameter(Parameter("TotalHeight", totalHeight, _brLenConstraint)); 
  for(map<const Node *, double>::iterator it = heights.begin(); it != heights.end(); it++)
  {
    if(!it->first->isLeaf() && it->first->hasFather())
    {
      Interval * constraint =  new IncludingInterval(0.,1.);
      _heightConstraints.push_back(constraint);
      _brLenParameters.addParameter(Parameter("HeightP" + TextTools::toString(it->first->getId()), it->second / totalHeight, constraint, false));
    }
  }
}

/******************************************************************************/

double ClockTreeLikelihood::computeBranchLengthsFromHeights(Node * node) throw (Exception)
{
  if(node->isLeaf()) return 0.;
  
  double totalHeight = _totalHeightParameter[0]->getValue();
  double nodeHeightP = 1.;
  if(node->hasFather()) //Not the root node
  {
    double fatherHeightP = 1.;
    if(node->getFather()->hasFather())
    {
      Node * father = node->getFather();
      Parameter * pp = _brLenParameters.getParameter(string("HeightP") + TextTools::toString(father->getId()));
      if(pp == NULL) throw Exception("ClockTreeLikelihood::computeBranchLengthsFromHeights(). Parameter HeightP" + TextTools::toString(father->getId()) + " was not found."); 
      fatherHeightP = pp->getValue();
    }

    //This parameter:
    Parameter * p = _brLenParameters.getParameter(string("HeightP") + TextTools::toString(node->getId()));
    if(p == NULL) throw Exception("ClockTreeLikelihood::computeBranchLengthsFromHeights(). Parameter HeightP" + TextTools::toString(node->getId()) + " was not found."); 
    nodeHeightP = p->getValue();

    double maxSonHeightP = 0.;
    for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
    {
      Node * son = node->getSon(i);
      if(!son->isLeaf())
      {
        Parameter * ps = _brLenParameters.getParameter(string("HeightP") + TextTools::toString(son->getId()));
        if(ps == NULL) throw Exception("ClockTreeLikelihood::computeBranchLengthsFromHeights(). Parameter " + string("HeightP") + TextTools::toString(son->getId()) + " was not found."); 
        double sonHeightP = ps->getValue();
        if(sonHeightP > maxSonHeightP)
          maxSonHeightP = sonHeightP;
      }
    }
    if(nodeHeightP < maxSonHeightP)
    {
      //throw Exception("ClockTreeLikelihood::computeBranchLengthsFromHeights(). Height " + TextTools::toString(node->getId()) + " is to small: " + TextTools::toString(nodeHeightP));
      //Conflict: setting height in between...
      //Height of son will be adjusted after recursive call.
      nodeHeightP = (maxSonHeightP + nodeHeightP) / 2.;
      p->setValue(nodeHeightP);
      if(_conflictingParameters.getParameter(p->getName()) == NULL)
        _conflictingParameters.addParameter(*p);
    }
    if(nodeHeightP > fatherHeightP)
    {
      //throw Exception("ClockTreeLikelihood::computeBranchLengthsFromHeights(). Height " + TextTools::toString(node->getId()) + " is to large:" + TextTools::toString(nodeHeightP));
      nodeHeightP = fatherHeightP;
      p->setValue(nodeHeightP);
      if(_conflictingParameters.getParameter(p->getName()) == NULL)
        _conflictingParameters.addParameter(*p);
    }
  }
  
  for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    Node * son = node->getSon(i);
    //Recursive call
    double sonHeightP = computeBranchLengthsFromHeights(son);
    double d = (nodeHeightP - sonHeightP) * totalHeight;
    if(d < _minimumBrLen) 
    {
      //cerr << "DEBUG: node=" << son->getId() << "\t" << d << endl;
      d = _minimumBrLen;
    }
    son->setDistanceToFather(d);
  }

  //Return the height for this node
  return nodeHeightP;
}

/******************************************************************************/

void ClockTreeLikelihood::adjustHeightsUp(const Node *node, double height)
{
  if(!node->hasFather()) return;
  if(!node->getFather()->hasFather()) return;
  const Node *father = node->getFather();

  //Father parameter:
  Parameter * pp = _brLenParameters.getParameter(string("HeightP") + TextTools::toString(father->getId()));
  if(pp == NULL) throw Exception("ClockTreeLikelihood::adjustHeightsUp(). Parameter HeightP" + TextTools::toString(father->getId()) + " was not found."); 
  double fatherHeightP = pp->getValue();
  if(fatherHeightP < height)
  {
    fatherHeightP = height;
    pp->setValue(fatherHeightP);
  }
  adjustHeightsUp(father, fatherHeightP);
}

void ClockTreeLikelihood::adjustHeightsDown(const Node *node, double height)
{
  for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    const Node *son = node->getSon(i);
    if(!son->isLeaf())
    {
      Parameter * ps = _brLenParameters.getParameter(string("HeightP") + TextTools::toString(son->getId()));
      if(ps == NULL) throw Exception("ClockTreeLikelihood::adjustHeightsUp(). Parameter HeightP" + TextTools::toString(son->getId()) + " was not found."); 
      double sonHeightP = ps->getValue();

      if(sonHeightP > height)
      {
        sonHeightP = height;
        ps->setValue(sonHeightP);
      }
      adjustHeightsDown(son, sonHeightP);
    }
  }
}

void ClockTreeLikelihood::adjustHeightsUp2(const Node *node, double ratio)
{
  if(!node->hasFather()) return;
  if(!node->getFather()->hasFather()) return;
  const Node *father = node->getFather();

  //Father parameter:
  Parameter * pp = _brLenParameters.getParameter(string("HeightP") + TextTools::toString(father->getId()));
  if(pp == NULL) throw Exception("ClockTreeLikelihood::adjustHeightsUp(). Parameter HeightP" + TextTools::toString(father->getId()) + " was not found."); 
  double fatherHeightP = pp->getValue();
  fatherHeightP = 1. - (1. - fatherHeightP) * ratio; 
  pp->setValue(fatherHeightP);
  adjustHeightsUp2(father, ratio);
}

void ClockTreeLikelihood::adjustHeightsDown2(const Node *node, double ratio)
{
  for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    const Node *son = node->getSon(i);
    if(!son->isLeaf())
    {
      //Son parameter:
      Parameter * ps = _brLenParameters.getParameter(string("HeightP") + TextTools::toString(son->getId()));
      if(ps == NULL) throw Exception("ClockTreeLikelihood::adjustHeightsUp(). Parameter HeightP" + TextTools::toString(son->getId()) + " was not found."); 
      double sonHeightP = ps->getValue();
      sonHeightP = sonHeightP * ratio; 
      ps->setValue(sonHeightP);
      adjustHeightsDown2(son, ratio);
    }
  }
}

/******************************************************************************/

void ClockTreeLikelihood::resetHeightsConstraints()
{
  resetHeightsConstraints(_tree->getRootNode());
}

/******************************************************************************/

void ClockTreeLikelihood::updateHeightsConstraints()
{
  updateHeightsConstraints(_tree->getRootNode());
}

/******************************************************************************/

void ClockTreeLikelihood::resetHeightsConstraints(const Node * node)
{
  if(node->isLeaf()) return;

  for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    //Recursive call
    const Node * son = node->getSon(i);
    resetHeightsConstraints(son);
  }

  if(node->hasFather()) //Not the root node
  {
    Parameter * p = _brLenParameters.getParameter(string("HeightP") + TextTools::toString(node->getId()));
    if(p == NULL) throw Exception("ClockTreeLikelihood::updateHeightsConstraints(). Parameter HeightP" + TextTools::toString(node->getId()) + " was not found."); 
    //Height is bound by the height of the father node from one side,
    //and the height of the closest descendant from the other side.
    //For optimization reason however, we allow them a minimum range of 2%.
    //Additionally, heights must be within range [0,1] (proportion of total height).
    Interval * constraint = dynamic_cast<Interval *>(p->getConstraint());
    constraint->setLowerBound(0.);
    constraint->setUpperBound(1.);
  }
}

/******************************************************************************/

double ClockTreeLikelihood::updateHeightsConstraints(const Node * node)
{
  if(node->isLeaf()) return _minimumBrLen;

  double maxSonHeightP = 0.;
  for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    //Recursive call
    const Node * son = node->getSon(i);
    double sonHeightP = updateHeightsConstraints(son);
    if(sonHeightP > maxSonHeightP) maxSonHeightP = sonHeightP;
  }

  double nodeHeightP = 1.;
  if(node->hasFather()) //Not the root node
  {
    Parameter * p = _brLenParameters.getParameter(string("HeightP") + TextTools::toString(node->getId()));
    if(p == NULL) throw Exception("ClockTreeLikelihood::updateHeightsConstraints(). Parameter HeightP" + TextTools::toString(node->getId()) + " was not found."); 
    nodeHeightP = p->getValue();
    double fatherHeightP = 1.;
    if(node->getFather()->hasFather())
    {
      const Node * grandFather = node->getFather();
      Parameter * pp = _brLenParameters.getParameter(string("HeightP") + TextTools::toString(grandFather->getId()));
      if(pp == NULL) throw Exception("ClockTreeLikelihood::updateHeightsConstraints(). Parameter HeightP" + TextTools::toString(grandFather->getId()) + " was not found."); 
      fatherHeightP = pp->getValue();
    }
    //Height is bound by the height of the father node from one side,
    //and the height of the closest descendant from the other side.
    //For optimization reason however, we allow them a minimum range of 2%.
    //Additionally, heights must be within range [0,1] (proportion of total height).
    double minBound = maxSonHeightP;
    double maxBound = fatherHeightP;
    if(nodeHeightP - maxSonHeightP < 0.01) minBound = nodeHeightP - 0.01; 
    if(fatherHeightP - nodeHeightP < 0.01) maxBound = nodeHeightP + 0.01;
    if(minBound < 0.) minBound = 0.;
    if(maxBound > 1.) maxBound = 1.;
    Interval * constraint = dynamic_cast<Interval *>(p->getConstraint());
    constraint->setLowerBound(minBound);
    constraint->setUpperBound(maxBound);
  }

  //Return the height of this node
  return nodeHeightP;
}

/******************************************************************************/

void ClockTreeLikelihood::fireParameterChanged(const ParameterList & params)
{
  if(params.size() == 1 && params[0]->getName().substr(0,7) == "HeightP")
  {
    //Only one branch length:
    
    //First method:
    //string focusHeight = params[0]->getName();
    //double oldHeightP = _brLenParameters.getParameter(focusHeight)->getValue();
    //_brLenParameters.matchParametersValues(_parameters);
    //double newHeightP = _brLenParameters.getParameter(focusHeight)->getValue();
    //const Node *focusNode = _tree->getNode(TextTools::toInt(focusHeight.substr(7)));
    //if(newHeightP > oldHeightP)
    //{
    //  adjustHeightsUp(focusNode, newHeightP);
    //}
    //else
    //{
    //  adjustHeightsDown(focusNode, newHeightP);
    //}
    
    //Second method:
    string focusHeight = params[0]->getName();
    double oldHeightP = _brLenParameters.getParameter(focusHeight)->getValue();
    _brLenParameters.matchParametersValues(_parameters);
    double newHeightP = _brLenParameters.getParameter(focusHeight)->getValue();
    const Node *focusNode = _tree->getNode(TextTools::toInt(focusHeight.substr(7)));
    if(newHeightP > oldHeightP)
    {
      double ratio = (1. - newHeightP) / (1. - oldHeightP);
      adjustHeightsUp2(focusNode, ratio);
    }
    else
    {
      double ratio = newHeightP / oldHeightP;
      adjustHeightsDown2(focusNode, ratio);
    }
    _parameters.matchParametersValues(_brLenParameters);
  }
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
  //updateHeightsConstraints();

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
  double totalHeight  = _totalHeightParameter[0]->getValue();
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
  double totalHeight  = _totalHeightParameter[0]->getValue();
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

