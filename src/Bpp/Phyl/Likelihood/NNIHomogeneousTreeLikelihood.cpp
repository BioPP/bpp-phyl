//
// File: DRHomogeneousTreeLikelihood.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 17 18:14:51 2003
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

#include "NNIHomogeneousTreeLikelihood.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/AutoParameter.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/*******************************************************************************/
void BranchLikelihood::initModel(const SubstitutionModel* model, const DiscreteDistribution* rDist)
{
  model_ = model;
  rDist_ = rDist;
  nbStates_ = model->getNumberOfStates();
  nbClasses_  = rDist->getNumberOfCategories();
  pxy_.resize(nbClasses_);
  for (size_t i = 0; i < nbClasses_; i++)
  {
    pxy_[i].resize(nbStates_);
    for (size_t j = 0; j < nbStates_; j++)
    {
      pxy_[i][j].resize(nbStates_);
    }
  }
}

/*******************************************************************************/
void BranchLikelihood::computeAllTransitionProbabilities()
{
  double l = getParameterValue("BrLen");

  // Computes all pxy once for all:
  for (size_t c = 0; c < nbClasses_; c++)
  {
    VVdouble* pxy__c = &pxy_[c];
    RowMatrix<double> Q = model_->getPij_t(l * rDist_->getCategory(c));
    for (size_t x = 0; x < nbStates_; x++)
    {
      Vdouble* pxy__c_x = &(*pxy__c)[x];
      for (size_t y = 0; y < nbStates_; y++)
      {
        (*pxy__c_x)[y] = Q(x, y);
      }
    }
  }
}

/*******************************************************************************/
void BranchLikelihood::computeLogLikelihood()
{
  lnL_ = 0;

  vector<double> la(array1_->size());
  for (size_t i = 0; i < array1_->size(); i++)
  {
    double Li = 0;
    for (size_t c = 0; c < nbClasses_; c++)
    {
      double rc = rDist_->getProbability(c);
      for (size_t x = 0; x < nbStates_; x++)
      {
        for (size_t y = 0; y < nbStates_; y++)
        {
          Li += rc * (*array1_)[i][c][x] * pxy_[c][x][y] * (*array2_)[i][c][y];
        }
      }
    }
    la[i] = weights_[i] * log(Li);
  }

  sort(la.begin(), la.end());
  for (size_t i = array1_->size(); i > 0; i--)
  {
    lnL_ -= la[i - 1];
  }
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood::NNIHomogeneousTreeLikelihood(
  const Tree& tree,
  SubstitutionModel* model,
  DiscreteDistribution* rDist,
  bool checkRooted,
  bool verbose)
throw (Exception) :
  DRHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  brLikFunction_(0),
  brentOptimizer_(0),
  brLenNNIValues_(),
  brLenNNIParams_()
{
  brentOptimizer_ = new BrentOneDimension();
  brentOptimizer_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  brentOptimizer_->setProfiler(0);
  brentOptimizer_->setMessageHandler(0);
  brentOptimizer_->setVerbose(0);
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood::NNIHomogeneousTreeLikelihood(
  const Tree& tree,
  const SiteContainer& data,
  SubstitutionModel* model,
  DiscreteDistribution* rDist,
  bool checkRooted,
  bool verbose)
throw (Exception) :
  DRHomogeneousTreeLikelihood(tree, data, model, rDist, checkRooted, verbose),
  brLikFunction_(0),
  brentOptimizer_(0),
  brLenNNIValues_(),
  brLenNNIParams_()
{
  brentOptimizer_ = new BrentOneDimension();
  brentOptimizer_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  brentOptimizer_->setProfiler(0);
  brentOptimizer_->setMessageHandler(0);
  brentOptimizer_->setVerbose(0);
  // We have to do this since the DRHomogeneousTreeLikelihood constructor will not call the overloaded setData method:
  brLikFunction_ = new BranchLikelihood(getLikelihoodData()->getWeights());
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood::NNIHomogeneousTreeLikelihood(const NNIHomogeneousTreeLikelihood& lik) :
  DRHomogeneousTreeLikelihood(lik),
  brLikFunction_(0),
  brentOptimizer_(0),
  brLenNNIValues_(),
  brLenNNIParams_()
{
  brLikFunction_  = dynamic_cast<BranchLikelihood*>(lik.brLikFunction_->clone());
  brentOptimizer_ = dynamic_cast<BrentOneDimension*>(lik.brentOptimizer_->clone());
  brLenNNIValues_ = lik.brLenNNIValues_;
  brLenNNIParams_ = lik.brLenNNIParams_;
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood& NNIHomogeneousTreeLikelihood::operator=(const NNIHomogeneousTreeLikelihood& lik)
{
  DRHomogeneousTreeLikelihood::operator=(lik);
  if (brLikFunction_) delete brLikFunction_;
  brLikFunction_  = dynamic_cast<BranchLikelihood*>(lik.brLikFunction_->clone());
  if (brentOptimizer_) delete brentOptimizer_;
  brentOptimizer_ = dynamic_cast<BrentOneDimension*>(lik.brentOptimizer_->clone());
  brLenNNIValues_ = lik.brLenNNIValues_;
  brLenNNIParams_ = lik.brLenNNIParams_;
  return *this;
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood::~NNIHomogeneousTreeLikelihood()
{
  if (brLikFunction_) delete brLikFunction_;
  delete brentOptimizer_;
}

/******************************************************************************/
double NNIHomogeneousTreeLikelihood::testNNI(int nodeId) const throw (NodeException)
{
  const Node* son    = tree_->getNode(nodeId);
  if (!son->hasFather()) throw NodePException("DRHomogeneousTreeLikelihood::testNNI(). Node 'son' must not be the root node.", son);
  const Node* parent = son->getFather();
  if (!parent->hasFather()) throw NodePException("DRHomogeneousTreeLikelihood::testNNI(). Node 'parent' must not be the root node.", parent);
  const Node* grandFather = parent->getFather();
  // From here: Bifurcation assumed.
  // In case of multifurcation, an arbitrary uncle is chosen.
  // If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
  size_t parentPosition = grandFather->getSonPosition(parent);
  // const Node * uncle = grandFather->getSon(parentPosition > 1 ? parentPosition - 1 : 1 - parentPosition);
  const Node* uncle = grandFather->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);

  // Retrieving arrays of interest:
  const DRASDRTreeLikelihoodNodeData* parentData = &getLikelihoodData()->getNodeData(parent->getId());
  const VVVdouble* sonArray   = &parentData->getLikelihoodArrayForNeighbor(son->getId());
  vector<const Node*> parentNeighbors = TreeTemplateTools::getRemainingNeighbors(parent, grandFather, son);
  size_t nbParentNeighbors = parentNeighbors.size();
  vector<const VVVdouble*> parentArrays(nbParentNeighbors);
  vector<const VVVdouble*> parentTProbs(nbParentNeighbors);
  for (size_t k = 0; k < nbParentNeighbors; k++)
  {
    const Node* n = parentNeighbors[k]; // This neighbor
    parentArrays[k] = &parentData->getLikelihoodArrayForNeighbor(n->getId());
    // if(n != grandFather) parentTProbs[k] = & pxy_[n->getId()];
    // else                 parentTProbs[k] = & pxy_[parent->getId()];
    parentTProbs[k] = &pxy_[n->getId()];
  }

  const DRASDRTreeLikelihoodNodeData* grandFatherData = &getLikelihoodData()->getNodeData(grandFather->getId());
  const VVVdouble* uncleArray      = &grandFatherData->getLikelihoodArrayForNeighbor(uncle->getId());
  vector<const Node*> grandFatherNeighbors = TreeTemplateTools::getRemainingNeighbors(grandFather, parent, uncle);
  size_t nbGrandFatherNeighbors = grandFatherNeighbors.size();
  vector<const VVVdouble*> grandFatherArrays;
  vector<const VVVdouble*> grandFatherTProbs;
  for (size_t k = 0; k < nbGrandFatherNeighbors; k++)
  {
    const Node* n = grandFatherNeighbors[k]; // This neighbor
    if (grandFather->getFather() == NULL || n != grandFather->getFather())
    {
      grandFatherArrays.push_back(&grandFatherData->getLikelihoodArrayForNeighbor(n->getId()));
      grandFatherTProbs.push_back(&pxy_[n->getId()]);
    }
  }

  // Compute array 1: grand father array
  VVVdouble array1 = *sonArray;
  resetLikelihoodArray(array1);
  grandFatherArrays.push_back(sonArray);
  grandFatherTProbs.push_back(&pxy_[son->getId()]);
  if (grandFather->hasFather())
  {
    computeLikelihoodFromArrays(grandFatherArrays, grandFatherTProbs, &grandFatherData->getLikelihoodArrayForNeighbor(grandFather->getFather()->getId()), &pxy_[grandFather->getId()], array1, nbGrandFatherNeighbors, nbDistinctSites_, nbClasses_, nbStates_, false);
  }
  else
  {
    computeLikelihoodFromArrays(grandFatherArrays, grandFatherTProbs, array1, nbGrandFatherNeighbors + 1, nbDistinctSites_, nbClasses_, nbStates_, false);

    // This is the root node, we have to account for the ancestral frequencies:
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      for (size_t j = 0; j < nbClasses_; j++)
      {
        for (size_t x = 0; x < nbStates_; x++)
        {
          array1[i][j][x] *= rootFreqs_[x];
        }
      }
    }
  }

  // Compute array 2: parent array
  VVVdouble array2 = *uncleArray;
  resetLikelihoodArray(array2);
  parentArrays.push_back(uncleArray);
  parentTProbs.push_back(&pxy_[uncle->getId()]);
  computeLikelihoodFromArrays(parentArrays, parentTProbs, array2, nbParentNeighbors + 1, nbDistinctSites_, nbClasses_, nbStates_, false);

  // Initialize BranchLikelihood:
  brLikFunction_->initModel(model_, rateDistribution_);
  brLikFunction_->initLikelihoods(&array1, &array2);
  ParameterList parameters;
  size_t pos = 0;
  while (pos < nodes_.size() && nodes_[pos]->getId() != parent->getId()) pos++;
  if (pos == nodes_.size()) throw Exception("NNIHomogeneousTreeLikelihood::testNNI. Unvalid node id.");
  Parameter brLen = getParameter("BrLen" + TextTools::toString(pos));
  brLen.setName("BrLen");
  parameters.addParameter(brLen);
  brLikFunction_->setParameters(parameters);

  // Re-estimate branch length:
  brentOptimizer_->setFunction(brLikFunction_);
  brentOptimizer_->getStopCondition()->setTolerance(0.1);
  brentOptimizer_->setInitialInterval(brLen.getValue(), brLen.getValue() + 0.01);
  brentOptimizer_->init(parameters);
  brentOptimizer_->optimize();
  // brLenNNIValues_[nodeId] = brLikFunction_->getParameterValue("BrLen");
  brLenNNIValues_[nodeId] = brentOptimizer_->getParameters().getParameter("BrLen").getValue();
  brLikFunction_->resetLikelihoods(); // Array1 and Array2 will be destroyed after this function call.
                                      // We should not keep pointers towards them...

  // Return the resulting likelihood:
  return brLikFunction_->getValue() - getValue();
}

/*******************************************************************************/
void NNIHomogeneousTreeLikelihood::doNNI(int nodeId) throw (NodeException)
{
  // Perform the topological move, the likelihood array will have to be recomputed...
  Node* son    = tree_->getNode(nodeId);
  if (!son->hasFather()) throw NodePException("DRHomogeneousTreeLikelihood::testNNI(). Node 'son' must not be the root node.", son);
  Node* parent = son->getFather();
  if (!parent->hasFather()) throw NodePException("DRHomogeneousTreeLikelihood::testNNI(). Node 'parent' must not be the root node.", parent);
  Node* grandFather = parent->getFather();
  // From here: Bifurcation assumed.
  // In case of multifurcation, an arbitrary uncle is chosen.
  // If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
  size_t parentPosition = grandFather->getSonPosition(parent);
  Node* uncle = grandFather->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);
  // Swap nodes:
  parent->removeSon(son);
  grandFather->removeSon(uncle);
  parent->addSon(uncle);
  grandFather->addSon(son);
  size_t pos = 0;
  while (pos < nodes_.size() && nodes_[pos]->getId() != parent->getId()) pos++;
  if (pos == nodes_.size()) throw Exception("NNIHomogeneousTreeLikelihood::doNNI. Unvalid node id.");

  string name = "BrLen" + TextTools::toString(pos);
  if (brLenNNIValues_.find(nodeId) != brLenNNIValues_.end())
  {
    double length = brLenNNIValues_[nodeId];
    brLenParameters_.setParameterValue(name, length);
    getParameter_(name).setValue(length);
    parent->setDistanceToFather(length);
  }
  else cerr << "ERROR, branch not found: " << nodeId << endl;
  try
  {
    brLenNNIParams_.addParameter(brLenParameters_.getParameter(name));
  }
  catch (ParameterException& ex)
  {
    StdErr errout;
    cerr << "DEBUG:" << endl;
    brLenNNIParams_.printParameters(errout);
    cerr << "DEBUG:" << name << endl;
  }
  // In case of copy of this object, we must remove the constraint associated to this stored parameter:
  // (It should be also possible to update the pointer in the copy constructor,
  // but we do not need the constraint info here...).
  brLenNNIParams_[brLenNNIParams_.size() - 1].removeConstraint();
}

/*******************************************************************************/

