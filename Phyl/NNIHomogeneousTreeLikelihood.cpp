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

#include "NNIHomogeneousTreeLikelihood.h"

// From Utils:
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>

// From NumCalc:
#include <NumCalc/AutoParameter.h>

// From the STL:
#include <iostream>
using namespace std;

/*******************************************************************************/

void BranchLikelihood::initModel(const SubstitutionModel *model, const DiscreteDistribution *rDist)
{
  this->_model = model;
  _rDist = rDist;
  _nbStates = model->getNumberOfStates();
  _nbClasses  = rDist->getNumberOfCategories();
  _pxy.resize(_nbClasses);
  for(unsigned int i = 0; i < _nbClasses; i++)
  {
    _pxy[i].resize(_nbStates);
    for(unsigned int j = 0; j < _nbStates; j++)
      _pxy[i][j].resize(_nbStates);
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
        //Li += rc * _model->freq(x) * (* arrayTmp_i_c)[x];
        //freq is already accounted in the array
        Li += rc * (* arrayTmp_i_c)[x];
      }
    }
    _lnL -= _weights[i] * log(Li);
  }
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood::NNIHomogeneousTreeLikelihood(
  const Tree & tree,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose)
throw (Exception):
  DRHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  _brLikFunction(NULL),
  _brentOptimizer(NULL),
  _brLenNNIValues(),
  _brLenNNIParams()
{
  _brentOptimizer = new BrentOneDimension();
  _brentOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  _brentOptimizer->setProfiler(NULL);
  _brentOptimizer->setMessageHandler(NULL);
  _brentOptimizer->setVerbose(0);
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood::NNIHomogeneousTreeLikelihood(
  const Tree & tree,
  const SiteContainer & data,
  SubstitutionModel * model,
  DiscreteDistribution * rDist,
  bool checkRooted,
  bool verbose)
throw (Exception):
  DRHomogeneousTreeLikelihood(tree, data, model, rDist, checkRooted, verbose),
  _brLikFunction(NULL),
  _brentOptimizer(NULL),
  _brLenNNIValues(),
  _brLenNNIParams()
{
  _brentOptimizer = new BrentOneDimension();
  _brentOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  _brentOptimizer->setProfiler(NULL);
  _brentOptimizer->setMessageHandler(NULL);
  _brentOptimizer->setVerbose(0);
  //We have to do this since the DRHomogeneousTreeLikelihood constructor will not call the overloaded setData method:
  _brLikFunction = new BranchLikelihood(_likelihoodData->getWeights());
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood::NNIHomogeneousTreeLikelihood(const NNIHomogeneousTreeLikelihood & lik):
  DRHomogeneousTreeLikelihood(lik),
  _brLikFunction(NULL),
  _brentOptimizer(NULL),
  _brLenNNIValues(),
  _brLenNNIParams()
{
  _brLikFunction  = dynamic_cast<BranchLikelihood *>(lik._brLikFunction->clone());
  _brentOptimizer = dynamic_cast<BrentOneDimension *>(lik._brentOptimizer->clone());
  _brLenNNIValues = lik._brLenNNIValues;
  _brLenNNIParams = lik._brLenNNIParams;
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood & NNIHomogeneousTreeLikelihood::operator=(const NNIHomogeneousTreeLikelihood & lik)
{
  DRHomogeneousTreeLikelihood::operator=(lik);
  if(_brLikFunction) delete _brLikFunction;
  _brLikFunction  = dynamic_cast<BranchLikelihood *>(lik._brLikFunction->clone());
  if(_brentOptimizer) delete _brentOptimizer;
  _brentOptimizer = dynamic_cast<BrentOneDimension *>(lik._brentOptimizer->clone());
  _brLenNNIValues = lik._brLenNNIValues;
  _brLenNNIParams = lik._brLenNNIParams;
  return *this;
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood::~NNIHomogeneousTreeLikelihood()
{
  if(_brLikFunction) delete _brLikFunction;
  delete _brentOptimizer;
}

/******************************************************************************/

double NNIHomogeneousTreeLikelihood::testNNI(int nodeId) const throw (NodeException)
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
  if(grandFather->getFather() == NULL)
  {
    //This is the root node, we have to account for the ancestral frequencies:
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      for(unsigned int j = 0; j < _nbClasses; j++)
      {
        for(unsigned int x = 0; x < _nbStates; x++)
          array1[i][j][x] *= _model->freq(x);
      }
    }
  }
  
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

void NNIHomogeneousTreeLikelihood::doNNI(int nodeId) throw (NodeException)
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
    Parameter * p = _parameters.getParameter(name);
    if(p) p->setValue(length);
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

