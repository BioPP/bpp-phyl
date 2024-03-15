// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/TextTools.h>

#include "../../PatternTools.h"
#include "AbstractHomogeneousTreeLikelihood.h"

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

/******************************************************************************/

AbstractHomogeneousTreeLikelihood::AbstractHomogeneousTreeLikelihood(
  const Tree& tree,
  shared_ptr<TransitionModelInterface> model,
  shared_ptr<DiscreteDistributionInterface> rDist,
  bool checkRooted,
  bool verbose) :
  AbstractDiscreteRatesAcrossSitesTreeLikelihood(rDist, verbose),
  model_(0),
  brLenParameters_(),
  pxy_(),
  dpxy_(),
  d2pxy_(),
  rootFreqs_(),
  nodes_(),
  nbSites_(),
  nbDistinctSites_(),
  nbClasses_(),
  nbStates_(),
  nbNodes_(),
  verbose_(),
  minimumBrLen_(),
  maximumBrLen_(),
  brLenConstraint_()
{
  init_(tree, model, rDist, checkRooted, verbose);
}

/******************************************************************************/

AbstractHomogeneousTreeLikelihood::AbstractHomogeneousTreeLikelihood(
  const AbstractHomogeneousTreeLikelihood& lik) :
  AbstractDiscreteRatesAcrossSitesTreeLikelihood(lik),
  model_(lik.model_),
  brLenParameters_(lik.brLenParameters_),
  pxy_(lik.pxy_),
  dpxy_(lik.dpxy_),
  d2pxy_(lik.d2pxy_),
  rootFreqs_(lik.rootFreqs_),
  nodes_(),
  nbSites_(lik.nbSites_),
  nbDistinctSites_(lik.nbDistinctSites_),
  nbClasses_(lik.nbClasses_),
  nbStates_(lik.nbStates_),
  nbNodes_(lik.nbNodes_),
  verbose_(lik.verbose_),
  minimumBrLen_(lik.minimumBrLen_),
  maximumBrLen_(lik.maximumBrLen_),
  brLenConstraint_(lik.brLenConstraint_->clone())
{
  nodes_ = tree_->getNodes();
  nodes_.pop_back(); // Remove the root node (the last added!).
}

/******************************************************************************/

AbstractHomogeneousTreeLikelihood& AbstractHomogeneousTreeLikelihood::operator=(
  const AbstractHomogeneousTreeLikelihood& lik)
{
  AbstractDiscreteRatesAcrossSitesTreeLikelihood::operator=(lik);
  model_           = lik.model_;
  brLenParameters_ = lik.brLenParameters_;
  pxy_             = lik.pxy_;
  dpxy_            = lik.dpxy_;
  d2pxy_           = lik.d2pxy_;
  rootFreqs_       = lik.rootFreqs_;
  nodes_ = tree_->getNodes();
  nodes_.pop_back(); // Remove the root node (the last added!).
  nbSites_         = lik.nbSites_;
  nbDistinctSites_ = lik.nbDistinctSites_;
  nbClasses_       = lik.nbClasses_;
  nbStates_        = lik.nbStates_;
  nbNodes_         = lik.nbNodes_;
  verbose_         = lik.verbose_;
  minimumBrLen_    = lik.minimumBrLen_;
  maximumBrLen_    = lik.maximumBrLen_;
  brLenConstraint_ = std::shared_ptr<ConstraintInterface>(lik.brLenConstraint_->clone());
  return *this;
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::init_(
  const Tree& tree,
  std::shared_ptr<TransitionModelInterface> model,
  std::shared_ptr<DiscreteDistributionInterface> rDist,
  bool checkRooted,
  bool verbose)
{
  TreeTools::checkIds(tree, true);
  tree_ = make_unique<TreeTemplate<Node>>(tree);
  if (checkRooted && tree_->isRooted())
  {
    if (verbose)
      ApplicationTools::displayWarning("Tree has been unrooted.");
    tree_->unroot();
  }
  nodes_ = tree_->getNodes();
  nodes_.pop_back(); // Remove the root node (the last added!).
  nbNodes_ = nodes_.size();
  nbClasses_ = rateDistribution_->getNumberOfCategories();
  setModel(model);

  verbose_ = verbose;

  minimumBrLen_ = 0.000001;
  maximumBrLen_ = 10000;
  brLenConstraint_ = std::make_shared<IntervalConstraint>(minimumBrLen_, maximumBrLen_, true, true);
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::setModel(std::shared_ptr<TransitionModelInterface> model)
{
  // Check:
  if (data_)
  {
    if (model->getAlphabet()->getAlphabetType() != data_->getAlphabet()->getAlphabetType())
      throw Exception("AbstractHomogeneousTreeLikelihood::setSubstitutionModel(). Model alphabet do not match existing data.");
  }

  model_ = model;

  if (data_)
  {
    if (model->getNumberOfStates() != model_->getNumberOfStates())
      setData(*data_);                                                  // Have to reinitialize the whole data structure.
  }

  nbStates_ = model->getNumberOfStates();

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
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::initialize()
{
  if (initialized_)
    throw Exception("AbstractHomogeneousTreeLikelihood::initialize(). Object is already initialized.");
  if (!data_)
    throw Exception("AbstractHomogeneousTreeLikelihood::initialize(). Data are no set.");
  initParameters();
  initialized_ = true;
  fireParameterChanged(getParameters());
}

/******************************************************************************/

ParameterList AbstractHomogeneousTreeLikelihood::getBranchLengthsParameters() const
{
  if (!initialized_)
    throw Exception("AbstractHomogeneousTreeLikelihood::getBranchLengthsParameters(). Object is not initialized.");
  return brLenParameters_.getCommonParametersWith(getParameters());
}

/******************************************************************************/

ParameterList AbstractHomogeneousTreeLikelihood::getSubstitutionModelParameters() const
{
  if (!initialized_)
    throw Exception("AbstractHomogeneousTreeLikelihood::getSubstitutionModelParameters(). Object is not initialized.");
  return model_->getParameters().getCommonParametersWith(getParameters());
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::initParameters()
{
  // Reset parameters:
  resetParameters_();

  // Branch lengths:
  initBranchLengthsParameters();
  addParameters_(brLenParameters_);

  // Substitution model:
  addParameters_(model_->getIndependentParameters());

  // Rate distribution:
  addParameters_(rateDistribution_->getIndependentParameters());
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::applyParameters()
{
  if (!initialized_)
    throw Exception("AbstractHomogeneousTreeLikelihood::applyParameters(). Object not initialized.");
  // Apply branch lengths:
  // brLenParameters_.matchParametersValues(parameters_); Not necessary!
  for (unsigned int i = 0; i < nbNodes_; i++)
  {
    const Parameter* brLen = &parameter(string("BrLen") + TextTools::toString(i));
    if (brLen)
      nodes_[i]->setDistanceToFather(brLen->getValue());
  }
  // Apply substitution model parameters:
  model_->matchParametersValues(getParameters());
  rootFreqs_ = model_->getFrequencies();
  // Apply rate distribution parameters:
  rateDistribution_->matchParametersValues(getParameters());
}

/******************************************************************************/

void AbstractHomogeneousTreeLikelihood::initBranchLengthsParameters(bool verbose)
{
  brLenParameters_.reset();
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
    brLenParameters_.addParameter(Parameter("BrLen" + TextTools::toString(i), d, brLenConstraint_));
  }
}

/*******************************************************************************/

void AbstractHomogeneousTreeLikelihood::computeAllTransitionProbabilities()
{
  for (unsigned int l = 0; l < nbNodes_; l++)
  {
    // For each node,
    Node* node = nodes_[l];
    computeTransitionProbabilitiesForNode(node);
  }
  rootFreqs_ = model_->getFrequencies();
}

/*******************************************************************************/

void AbstractHomogeneousTreeLikelihood::computeTransitionProbabilitiesForNode(const Node* node)
{
  double l = node->getDistanceToFather();

  // Computes all pxy and pyx once for all:
  VVVdouble* pxy__node = &pxy_[node->getId()];
  for (unsigned int c = 0; c < nbClasses_; c++)
  {
    VVdouble* pxy__node_c = &(*pxy__node)[c];
    RowMatrix<double> Q = model_->getPij_t(l * rateDistribution_->getCategory(c));

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
      RowMatrix<double> dQ = model_->getdPij_dt(l * rc);
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
      RowMatrix<double> d2Q = model_->getd2Pij_dt2(l * rc);
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
