// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Text/TextTools.h>

#include "NNIHomogeneousTreeLikelihood.h"

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/*******************************************************************************/
void BranchLikelihood::initModel(
    std::shared_ptr<const TransitionModelInterface> model,
    std::shared_ptr<const DiscreteDistributionInterface> rDist)
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
    std::shared_ptr<TransitionModelInterface> model,
    std::shared_ptr<DiscreteDistributionInterface> rDist,
    bool checkRooted,
    bool verbose) :
  DRHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose),
  brLikFunction_(),
  brentOptimizer_(),
  brLenNNIValues_(),
  brLenNNIParams_()
{
  brentOptimizer_ = make_unique<BrentOneDimension>();
  brentOptimizer_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  brentOptimizer_->setProfiler(0);
  brentOptimizer_->setMessageHandler(0);
  brentOptimizer_->setVerbose(0);
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood::NNIHomogeneousTreeLikelihood(
    const Tree& tree,
    const AlignmentDataInterface& data,
    std::shared_ptr<TransitionModelInterface> model,
    std::shared_ptr<DiscreteDistributionInterface> rDist,
    bool checkRooted,
    bool verbose) :
  DRHomogeneousTreeLikelihood(tree, data, model, rDist, checkRooted, verbose),
  brLikFunction_(),
  brentOptimizer_(),
  brLenNNIValues_(),
  brLenNNIParams_()
{
  brentOptimizer_ = make_unique<BrentOneDimension>();
  brentOptimizer_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  brentOptimizer_->setProfiler(0);
  brentOptimizer_->setMessageHandler(0);
  brentOptimizer_->setVerbose(0);
  // We have to do this since the DRHomogeneousTreeLikelihood constructor will not call the overloaded setData method:
  brLikFunction_ = make_shared<BranchLikelihood>(likelihoodData().getWeights());
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood::NNIHomogeneousTreeLikelihood(const NNIHomogeneousTreeLikelihood& lik) :
  DRHomogeneousTreeLikelihood(lik),
  brLikFunction_(),
  brentOptimizer_(),
  brLenNNIValues_(),
  brLenNNIParams_()
{
  brLikFunction_  = shared_ptr<BranchLikelihood>(lik.brLikFunction_->clone());
  brentOptimizer_ = unique_ptr<BrentOneDimension>(lik.brentOptimizer_->clone());
  brLenNNIValues_ = lik.brLenNNIValues_;
  brLenNNIParams_ = lik.brLenNNIParams_;
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood& NNIHomogeneousTreeLikelihood::operator=(const NNIHomogeneousTreeLikelihood& lik)
{
  DRHomogeneousTreeLikelihood::operator=(lik);
  brLikFunction_  = shared_ptr<BranchLikelihood>(lik.brLikFunction_->clone());
  brentOptimizer_ = unique_ptr<BrentOneDimension>(lik.brentOptimizer_->clone());
  brLenNNIValues_ = lik.brLenNNIValues_;
  brLenNNIParams_ = lik.brLenNNIParams_;
  return *this;
}

/******************************************************************************/

NNIHomogeneousTreeLikelihood::~NNIHomogeneousTreeLikelihood() {}

/******************************************************************************/

double NNIHomogeneousTreeLikelihood::testNNI(int nodeId) const
{
  const Node* son    = tree_->getNode(nodeId);
  if (!son->hasFather())
    throw NodePException("DRHomogeneousTreeLikelihood::testNNI(). Node 'son' must not be the root node.", son);
  const Node* parent = son->getFather();
  if (!parent->hasFather())
    throw NodePException("DRHomogeneousTreeLikelihood::testNNI(). Node 'parent' must not be the root node.", parent);
  const Node* grandFather = parent->getFather();
  // From here: Bifurcation assumed.
  // In case of multifurcation, an arbitrary uncle is chosen.
  // If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
  size_t parentPosition = grandFather->getSonPosition(parent);
  // const Node * uncle = grandFather->getSon(parentPosition > 1 ? parentPosition - 1 : 1 - parentPosition);
  const Node* uncle = grandFather->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);

  // Retrieving arrays of interest:
  const DRASDRTreeLikelihoodNodeData* parentData = &likelihoodData().getNodeData(parent->getId());
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

  const DRASDRTreeLikelihoodNodeData* grandFatherData = &likelihoodData().getNodeData(grandFather->getId());
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
  while (pos < nodes_.size() && nodes_[pos]->getId() != parent->getId())
    pos++;
  if (pos == nodes_.size())
    throw Exception("NNIHomogeneousTreeLikelihood::testNNI. Invalid node id.");
  Parameter brLen = parameter("BrLen" + TextTools::toString(pos));
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
  brLenNNIValues_[nodeId] = brentOptimizer_->getParameters().parameter("BrLen").getValue();
  brLikFunction_->resetLikelihoods(); // Array1 and Array2 will be destroyed after this function call.
                                      // We should not keep pointers towards them...

  // Return the resulting likelihood:
  return brLikFunction_->getValue() - getValue();
}

/*******************************************************************************/
void NNIHomogeneousTreeLikelihood::doNNI(int nodeId)
{
  // Perform the topological move, the likelihood array will have to be recomputed...
  Node* son    = tree_->getNode(nodeId);
  if (!son->hasFather())
    throw NodePException("DRHomogeneousTreeLikelihood::testNNI(). Node 'son' must not be the root node.", son);
  Node* parent = son->getFather();
  if (!parent->hasFather())
    throw NodePException("DRHomogeneousTreeLikelihood::testNNI(). Node 'parent' must not be the root node.", parent);
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
  while (pos < nodes_.size() && nodes_[pos]->getId() != parent->getId())
    pos++;
  if (pos == nodes_.size())
    throw Exception("NNIHomogeneousTreeLikelihood::doNNI. Invalid node id.");

  string name = "BrLen" + TextTools::toString(pos);
  if (brLenNNIValues_.find(nodeId) != brLenNNIValues_.end())
  {
    double length = brLenNNIValues_[nodeId];
    brLenParameters_.setParameterValue(name, length);
    getParameter_(name).setValue(length);
    parent->setDistanceToFather(length);
  }
  else
    cerr << "ERROR, branch not found: " << nodeId << endl;
  try
  {
    brLenNNIParams_.addParameter(brLenParameters_.parameter(name));
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
