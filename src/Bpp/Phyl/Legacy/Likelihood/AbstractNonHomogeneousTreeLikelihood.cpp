// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../../PatternTools.h"
#include "AbstractNonHomogeneousTreeLikelihood.h"

// From SeqLib:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

AbstractNonHomogeneousTreeLikelihood::AbstractNonHomogeneousTreeLikelihood(
  const Tree& tree,
  shared_ptr<SubstitutionModelSet> modelSet,
  shared_ptr<DiscreteDistributionInterface> rDist,
  bool verbose,
  bool reparametrizeRoot) :
  AbstractDiscreteRatesAcrossSitesTreeLikelihood(rDist, verbose),
  modelSet_(0),
  brLenParameters_(),
  pxy_(),
  dpxy_(),
  d2pxy_(),
  rootFreqs_(),
  nodes_(),
  idToNode_(),
  nbSites_(),
  nbDistinctSites_(),
  nbClasses_(),
  nbStates_(),
  nbNodes_(),
  verbose_(),
  minimumBrLen_(),
  maximumBrLen_(),
  brLenConstraint_(),
  reparametrizeRoot_(reparametrizeRoot),
  root1_(),
  root2_()
{
  init_(tree, modelSet, rDist, verbose);
}

/******************************************************************************/

AbstractNonHomogeneousTreeLikelihood::AbstractNonHomogeneousTreeLikelihood(
  const AbstractNonHomogeneousTreeLikelihood& lik) :
  AbstractDiscreteRatesAcrossSitesTreeLikelihood(lik),
  modelSet_(lik.modelSet_),
  brLenParameters_(lik.brLenParameters_),
  pxy_(lik.pxy_),
  dpxy_(lik.dpxy_),
  d2pxy_(lik.d2pxy_),
  rootFreqs_(lik.rootFreqs_),
  nodes_(),
  idToNode_(),
  nbSites_(lik.nbSites_),
  nbDistinctSites_(lik.nbDistinctSites_),
  nbClasses_(lik.nbClasses_),
  nbStates_(lik.nbStates_),
  nbNodes_(lik.nbNodes_),
  verbose_(lik.verbose_),
  minimumBrLen_(lik.minimumBrLen_),
  maximumBrLen_(lik.maximumBrLen_),
  brLenConstraint_(lik.brLenConstraint_->clone()),
  reparametrizeRoot_(lik.reparametrizeRoot_),
  root1_(lik.root1_),
  root2_(lik.root2_)
{
  nodes_ = tree_->getNodes();
  nodes_.pop_back(); // Remove the root node (the last added!).
  // Rebuild nodes index:
  for (unsigned int i = 0; i < nodes_.size(); i++)
  {
    const Node* node = nodes_[i];
    idToNode_[node->getId()] = node;
  }
}

/******************************************************************************/

AbstractNonHomogeneousTreeLikelihood& AbstractNonHomogeneousTreeLikelihood::operator=(
  const AbstractNonHomogeneousTreeLikelihood& lik)
{
  AbstractDiscreteRatesAcrossSitesTreeLikelihood::operator=(lik);
  modelSet_          = lik.modelSet_;
  brLenParameters_   = lik.brLenParameters_;
  pxy_               = lik.pxy_;
  dpxy_              = lik.dpxy_;
  d2pxy_             = lik.d2pxy_;
  rootFreqs_         = lik.rootFreqs_;
  nodes_             = tree_->getNodes();
  nodes_.pop_back(); // Remove the root node (the last added!).
  nbSites_           = lik.nbSites_;
  nbDistinctSites_   = lik.nbDistinctSites_;
  nbClasses_         = lik.nbClasses_;
  nbStates_          = lik.nbStates_;
  nbNodes_           = lik.nbNodes_;
  verbose_           = lik.verbose_;
  minimumBrLen_      = lik.minimumBrLen_;
  maximumBrLen_      = lik.maximumBrLen_;
  brLenConstraint_   = std::shared_ptr<ConstraintInterface>(lik.brLenConstraint_->clone());
  reparametrizeRoot_ = lik.reparametrizeRoot_;
  root1_             = lik.root1_;
  root2_             = lik.root2_;
  // Rebuild nodes index:
  for (unsigned int i = 0; i < nodes_.size(); i++)
  {
    const Node* node = nodes_[i];
    idToNode_[node->getId()] = node;
  }
  return *this;
}

/******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::init_(
  const Tree& tree,
  shared_ptr<SubstitutionModelSet> modelSet,
  shared_ptr<DiscreteDistributionInterface> rDist,
  bool verbose)
{
  TreeTools::checkIds(tree, true);
  tree_ = make_unique< TreeTemplate<Node> >(tree);
  root1_ = tree_->getRootNode()->getSon(0)->getId();
  root2_ = tree_->getRootNode()->getSon(1)->getId();
  nodes_ = tree_->getNodes();
  nodes_.pop_back(); // Remove the root node (the last added!).
  nbNodes_ = nodes_.size();
  // Build nodes index:
  for (size_t i = 0; i < nodes_.size(); ++i)
  {
    const Node* node = nodes_[i];
    idToNode_[node->getId()] = node;
  }
  nbClasses_ = rateDistribution_->getNumberOfCategories();

  verbose_ = verbose;

  minimumBrLen_ = 0.000001;
  maximumBrLen_ = 10000;
  brLenConstraint_ = std::make_shared<IntervalConstraint>(minimumBrLen_, maximumBrLen_, true, true);
  setSubstitutionModelSet(modelSet);
}

/******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::setSubstitutionModelSet(
     shared_ptr<SubstitutionModelSet> modelSet
     )
{
  // Check:
  if (data_)
  {
    if (modelSet->getAlphabet()->getAlphabetType() != data_->getAlphabet()->getAlphabetType())
      throw Exception("AbstractBranchNonHomogeneousTreeLikelihood::setSubstitutionModelSet(). Model alphabet do not match existing data.");
  }

  modelSet_ = modelSet;

  if (data_)
  {
    if (modelSet->getNumberOfStates() != modelSet_->getNumberOfStates())
      setData(*data_); // Have to reinitialize the whole data structure.
  }

  nbStates_ = modelSet->getNumberOfStates();

  // Allocate transition probabilities arrays:
  for (unsigned int l = 0; l < nbNodes_; l++)
  {
    // For each son node,i
    Node* son = nodes_[l];

    VVVdouble* pxy__son = &pxy_[son->getId()];
    pxy__son->resize(nbClasses_);
    for (unsigned int c = 0; c < nbClasses_; c++)
    {
      VVdouble* pxy__son_c = &(*pxy__son)[c];
      pxy__son_c->resize(nbStates_);
      for (unsigned int x = 0; x < nbStates_; x++)
      {
        (*pxy__son_c)[x].resize(nbStates_);
      }
    }

    VVVdouble* dpxy__son = &dpxy_[son->getId()];
    dpxy__son->resize(nbClasses_);
    for (unsigned int c = 0; c < nbClasses_; c++)
    {
      VVdouble* dpxy__son_c = &(*dpxy__son)[c];
      dpxy__son_c->resize(nbStates_);
      for (unsigned int x = 0; x < nbStates_; x++)
      {
        (*dpxy__son_c)[x].resize(nbStates_);
      }
    }

    VVVdouble* d2pxy__son = &d2pxy_[son->getId()];
    d2pxy__son->resize(nbClasses_);
    for (unsigned int c = 0; c < nbClasses_; c++)
    {
      VVdouble* d2pxy__son_c = &(*d2pxy__son)[c];
      d2pxy__son_c->resize(nbStates_);
      for (unsigned int x = 0; x < nbStates_; x++)
      {
        (*d2pxy__son_c)[x].resize(nbStates_);
      }
    }
  }

  // We have to reset parameters. If the instance is not initialized, this will be done by the initialize method.
  if (initialized_)
  {
    initParameters();
    computeAllTransitionProbabilities();
    fireParameterChanged(getParameters());
  }
}

/******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::initialize()
{
  if (initialized_)
    throw Exception("AbstractBranchNonHomogeneousTreeLikelihood::initialize(). Object is already initialized.");
  if (!data_)
    throw Exception("AbstractBranchNonHomogeneousTreeLikelihood::initialize(). Data are no set.");
  initParameters();
  initialized_ = true;

  computeAllTransitionProbabilities();
  fireParameterChanged(getParameters());
}

/******************************************************************************/

ParameterList AbstractNonHomogeneousTreeLikelihood::getBranchLengthsParameters() const
{
  if (!initialized_)
    throw Exception("AbstractBranchNonHomogeneousTreeLikelihood::getBranchLengthsParameters(). Object is not initialized.");
  return brLenParameters_.getCommonParametersWith(getParameters());
}

/******************************************************************************/

ParameterList AbstractNonHomogeneousTreeLikelihood::getSubstitutionModelParameters() const
{
  if (!initialized_)
    throw Exception("AbstractBranchNonHomogeneousTreeLikelihood::getSubstitutionModelParameters(). Object is not initialized.");
  return modelSet_->getParameters().getCommonParametersWith(getParameters());
}

/******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::initParameters()
{
  // Reset parameters:
  resetParameters_();

  // Branch lengths:
  initBranchLengthsParameters();
  addParameters_(brLenParameters_);

  // Substitution model:
  addParameters_(modelSet_->getIndependentParameters());

  // Rate distribution:
  addParameters_(rateDistribution_->getIndependentParameters());
}

/******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::applyParameters()
{
  if (!initialized_)
    throw Exception("AbstractBranchNonHomogeneousTreeLikelihood::applyParameters(). Object not initialized.");
  // Apply branch lengths:
  for (unsigned int i = 0; i < nbNodes_; i++)
  {
    int id = nodes_[i]->getId();
    if (reparametrizeRoot_ && id == root1_)
    {
      const Parameter* rootBrLen = &parameter("BrLenRoot");
      const Parameter* rootPos = &parameter("RootPosition");
      nodes_[i]->setDistanceToFather(rootBrLen->getValue() * rootPos->getValue());
    }
    else if (reparametrizeRoot_ && id == root2_)
    {
      const Parameter* rootBrLen = &parameter("BrLenRoot");
      const Parameter* rootPos = &parameter("RootPosition");
      nodes_[i]->setDistanceToFather(rootBrLen->getValue() * (1. - rootPos->getValue()));
    }
    else
    {
      const Parameter* brLen = &parameter(string("BrLen") + TextTools::toString(i));
      if (brLen)
        nodes_[i]->setDistanceToFather(brLen->getValue());
    }
  }
  // Apply substitution model parameters:

  modelSet_->matchParametersValues(getParameters());
  // Apply rate distribution parameters:
  rateDistribution_->matchParametersValues(getParameters());
}

/******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::initBranchLengthsParameters(bool verbose)
{
  brLenParameters_.reset();
  double l1 = 0, l2 = 0;
  for (unsigned int i = 0; i < nbNodes_; i++)
  {
    double d = minimumBrLen_;
    if (!nodes_[i]->hasDistanceToFather())
    {
      if (verbose)
        ApplicationTools::displayWarning("Missing branch length " + TextTools::toString(i) + ". Value is set to " + TextTools::toString(minimumBrLen_));
      nodes_[i]->setDistanceToFather(minimumBrLen_);
    }
    else
    {
      d = nodes_[i]->getDistanceToFather();
      if (d < minimumBrLen_)
      {
        if (verbose)
          ApplicationTools::displayWarning("Branch length " + TextTools::toString(i) + " is too small: " + TextTools::toString(d) + ". Value is set to " + TextTools::toString(minimumBrLen_));
        nodes_[i]->setDistanceToFather(minimumBrLen_);
        d = minimumBrLen_;
      }
      if (d > maximumBrLen_)
      {
        if (verbose)
          ApplicationTools::displayWarning("Branch length " + TextTools::toString(i) + " is too big: " + TextTools::toString(d) + ". Value is set to " + TextTools::toString(maximumBrLen_));
        nodes_[i]->setDistanceToFather(maximumBrLen_);
        d = maximumBrLen_;
      }
    }
    if (reparametrizeRoot_ && nodes_[i]->getId() == root1_)
      l1 = d;
    else if (reparametrizeRoot_ && nodes_[i]->getId() == root2_)
      l2 = d;
    else
    {
      brLenParameters_.addParameter(Parameter("BrLen" + TextTools::toString(i), d, brLenConstraint_));
    }
  }
  if (reparametrizeRoot_)
  {
    brLenParameters_.addParameter(Parameter("BrLenRoot", l1 + l2, brLenConstraint_));
    brLenParameters_.addParameter(Parameter("RootPosition", l1 / (l1 + l2), Parameter::PROP_CONSTRAINT_EX));
  }
}

/*******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::computeAllTransitionProbabilities()
{
  for (unsigned int l = 0; l < nbNodes_; l++)
  {
    // For each node,
    Node* node = nodes_[l];
    computeTransitionProbabilitiesForNode(node);
  }
  rootFreqs_ = modelSet_->getRootFrequencies();
}

/*******************************************************************************/

void AbstractNonHomogeneousTreeLikelihood::computeTransitionProbabilitiesForNode(const Node* node)
{
  auto model = modelSet_->getModelForNode(node->getId());
  double l = node->getDistanceToFather();

  // Computes all pxy and pyx once for all:
  VVVdouble* pxy__node = &pxy_[node->getId()];
  for (unsigned int c = 0; c < nbClasses_; c++)
  {
    VVdouble* pxy__node_c = &(*pxy__node)[c];
    RowMatrix<double> Q = model->getPij_t(l * rateDistribution_->getCategory(c));
    for (unsigned int x = 0; x < nbStates_; x++)
    {
      Vdouble* pxy__node_c_x = &(*pxy__node_c)[x];
      for (unsigned int y = 0; y < nbStates_; y++)
      {
        (*pxy__node_c_x)[y] = Q(x, y);
      }
    }
  }

  if (computeFirstOrderDerivatives_)
  {
    // Computes all dpxy/dt once for all:
    VVVdouble* dpxy__node = &dpxy_[node->getId()];

    for (unsigned int c = 0; c < nbClasses_; c++)
    {
      VVdouble* dpxy__node_c = &(*dpxy__node)[c];
      double rc = rateDistribution_->getCategory(c);

      RowMatrix<double> dQ = model->getdPij_dt(l * rc);

      for (unsigned int x = 0; x < nbStates_; x++)
      {
        Vdouble* dpxy__node_c_x = &(*dpxy__node_c)[x];
        for (unsigned int y = 0; y < nbStates_; y++)
        {
          (*dpxy__node_c_x)[y] = rc * dQ(x, y);
        }
      }
    }
  }

  if (computeSecondOrderDerivatives_)
  {
    // Computes all d2pxy/dt2 once for all:
    VVVdouble* d2pxy__node = &d2pxy_[node->getId()];
    for (unsigned int c = 0; c < nbClasses_; c++)
    {
      VVdouble* d2pxy__node_c = &(*d2pxy__node)[c];
      double rc =  rateDistribution_->getCategory(c);
      RowMatrix<double> d2Q = model->getd2Pij_dt2(l * rc);
      for (unsigned int x = 0; x < nbStates_; x++)
      {
        Vdouble* d2pxy__node_c_x = &(*d2pxy__node_c)[x];
        for (unsigned int y = 0; y < nbStates_; y++)
        {
          (*d2pxy__node_c_x)[y] = rc * rc * d2Q(x, y);
        }
      }
    }
  }
}

/*******************************************************************************/
