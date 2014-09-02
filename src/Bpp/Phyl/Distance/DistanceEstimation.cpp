//
// File: DistanceEstimation.cpp
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Wed jun 08 10:39 2005
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

#include "DistanceEstimation.h"
#include "../Tree/Tree.h"
#include "../PatternTools.h"
#include "../SitePatterns.h"

// From bpp-core:
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/AutoParameter.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/DistanceMatrix.h>

using namespace bpp;

// From the STL:
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

/******************************************************************************/

TwoTreeLikelihood::TwoTreeLikelihood(
    const std::string& seq1, const std::string& seq2,
    const SiteContainer& data,
    SubstitutionModel* model,
    DiscreteDistribution* rDist,
    bool verbose) throw (Exception) :
  AbstractDiscreteRatesAcrossSitesTreeLikelihood(rDist, verbose),
  shrunkData_(0), seqnames_(2), model_(model), brLenParameters_(), pxy_(), dpxy_(), d2pxy_(),
  rootPatternLinks_(), rootWeights_(), nbSites_(0), nbClasses_(0), nbStates_(0), nbDistinctSites_(0),
  rootLikelihoods_(), rootLikelihoodsS_(), rootLikelihoodsSR_(), dLikelihoods_(), d2Likelihoods_(),
  leafLikelihoods1_(), leafLikelihoods2_(),
  minimumBrLen_(0.000001), brLenConstraint_(0), brLen_(0)
{
  seqnames_[0] = seq1;
  seqnames_[1] = seq2;
  data_ = PatternTools::getSequenceSubset(data, seqnames_);
  if (data_->getAlphabet()->getAlphabetType()
      != model_->getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("TwoTreeTreeLikelihood::TwoTreeTreeLikelihood. Data and model must have the same alphabet type.",
                                    data_->getAlphabet(),
                                    model_->getAlphabet());

  nbSites_   = data_->getNumberOfSites();
  nbClasses_ = rateDistribution_->getNumberOfCategories();
  nbStates_  = model_->getNumberOfStates();
  if (verbose)
    ApplicationTools::displayMessage("Double-Recursive Homogeneous Tree Likelihood");

  // Initialize root patterns:
  SitePatterns pattern(data_);
  shrunkData_       = pattern.getSites();
  rootWeights_      = pattern.getWeights();
  rootPatternLinks_ = pattern.getIndices();
  nbDistinctSites_  = shrunkData_->getNumberOfSites();
  if (verbose)
    ApplicationTools::displayResult("Number of distinct sites", TextTools::toString(nbDistinctSites_));

  // Init _likelihoods:
  if (verbose) ApplicationTools::displayTask("Init likelihoods arrays recursively");
  // Clone data for more efficiency on sequences access:
  const SiteContainer* sequences = new AlignedSequenceContainer(*shrunkData_);
  initTreeLikelihoods(*sequences);
  delete sequences;

  brLen_ = minimumBrLen_;
  brLenConstraint_ = new IntervalConstraint(1, minimumBrLen_, true);

  if (verbose) ApplicationTools::displayTaskDone();
}

/******************************************************************************/

TwoTreeLikelihood::TwoTreeLikelihood(const TwoTreeLikelihood& lik) :
  AbstractDiscreteRatesAcrossSitesTreeLikelihood(lik),
  shrunkData_        (dynamic_cast<SiteContainer*>(lik.shrunkData_->clone())),
  seqnames_          (lik.seqnames_),
  model_             (lik.model_),
  brLenParameters_   (lik.brLenParameters_),
  pxy_               (lik.pxy_),
  dpxy_              (lik.dpxy_),
  d2pxy_             (lik.d2pxy_),
  rootPatternLinks_  (lik.rootPatternLinks_),
  rootWeights_       (lik.rootWeights_),
  nbSites_           (lik.nbSites_),
  nbClasses_         (lik.nbClasses_),
  nbStates_          (lik.nbStates_),
  nbDistinctSites_   (lik.nbDistinctSites_),
  rootLikelihoods_   (lik.rootLikelihoods_),
  rootLikelihoodsS_  (lik.rootLikelihoodsS_),
  rootLikelihoodsSR_ (lik.rootLikelihoodsSR_),
  dLikelihoods_      (lik.dLikelihoods_),
  d2Likelihoods_     (lik.d2Likelihoods_),
  leafLikelihoods1_  (lik.leafLikelihoods1_),
  leafLikelihoods2_  (lik.leafLikelihoods2_),
  minimumBrLen_      (lik.minimumBrLen_),
  brLenConstraint_   (dynamic_cast<Constraint*>(lik.brLenConstraint_->clone())),
  brLen_             (lik.brLen_)
{}

/******************************************************************************/

TwoTreeLikelihood& TwoTreeLikelihood::operator=(const TwoTreeLikelihood& lik)
{
  AbstractDiscreteRatesAcrossSitesTreeLikelihood::operator=(lik);
  shrunkData_        = dynamic_cast<SiteContainer*>(lik.shrunkData_->clone());
  seqnames_          = lik.seqnames_;
  model_             = lik.model_;
  brLenParameters_   = lik.brLenParameters_;
  pxy_               = lik.pxy_;
  dpxy_              = lik.dpxy_;
  d2pxy_             = lik.d2pxy_;
  rootPatternLinks_  = lik.rootPatternLinks_;
  rootWeights_       = lik.rootWeights_;
  nbSites_           = lik.nbSites_;
  nbClasses_         = lik.nbClasses_;
  nbStates_          = lik.nbStates_;
  nbDistinctSites_   = lik.nbDistinctSites_;
  rootLikelihoods_   = lik.rootLikelihoods_;
  rootLikelihoodsS_  = lik.rootLikelihoodsS_;
  rootLikelihoodsSR_ = lik.rootLikelihoodsSR_;
  dLikelihoods_      = lik.dLikelihoods_;
  d2Likelihoods_     = lik.d2Likelihoods_;
  leafLikelihoods1_  = lik.leafLikelihoods1_;
  leafLikelihoods2_  = lik.leafLikelihoods2_;
  minimumBrLen_      = lik.minimumBrLen_;
  brLenConstraint_   = dynamic_cast<Constraint*>(brLenConstraint_->clone());
  brLen_             = lik.brLen_;
  return *this;
}

/******************************************************************************/

TwoTreeLikelihood::~TwoTreeLikelihood()
{
  delete shrunkData_;
  if (brLenConstraint_) delete brLenConstraint_;
}

/******************************************************************************/

void TwoTreeLikelihood::initialize() throw (Exception)
{
  initParameters();
  initialized_ = true;
  fireParameterChanged(getParameters());
}

/******************************************************************************/

ParameterList TwoTreeLikelihood::getBranchLengthsParameters() const
{
  if (!initialized_) throw Exception("TwoTreeLikelihood::getBranchLengthsParameters(). Object is not initialized.");
  return brLenParameters_.getCommonParametersWith(getParameters());
}

/******************************************************************************/

ParameterList TwoTreeLikelihood::getSubstitutionModelParameters() const
{
  if (!initialized_) throw Exception("TwoTreeLikelihood::getSubstitutionModelParameters(). Object is not initialized.");
  return model_->getParameters().getCommonParametersWith(getParameters());
}

/******************************************************************************/

double TwoTreeLikelihood::getLikelihood() const
{
  double l = 1.;
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    l *= std::pow(rootLikelihoodsSR_[i], (int)rootWeights_[i]);
  }
  return l;
}

/******************************************************************************/

double TwoTreeLikelihood::getLogLikelihood() const
{
  double ll = 0;
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    ll += rootWeights_[i] * log(rootLikelihoodsSR_[i]);
  }
  return ll;
}

/******************************************************************************/

double TwoTreeLikelihood::getLikelihoodForASite(size_t site) const
{
  return rootLikelihoodsSR_[rootPatternLinks_[site]];
}

/******************************************************************************/

double TwoTreeLikelihood::getLogLikelihoodForASite(size_t site) const
{
  return log(rootLikelihoodsSR_[rootPatternLinks_[site]]);
}

/******************************************************************************/

double TwoTreeLikelihood::getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  return rootLikelihoodsS_[rootPatternLinks_[site]][rateClass];
}

/******************************************************************************/

double TwoTreeLikelihood::getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
  return log(rootLikelihoodsS_[rootPatternLinks_[site]][rateClass]);
}

/******************************************************************************/

double TwoTreeLikelihood::getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  return rootLikelihoods_[rootPatternLinks_[site]][rateClass][static_cast<size_t>(state)];
}

/******************************************************************************/

double TwoTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
  return log(rootLikelihoods_[rootPatternLinks_[site]][rateClass][static_cast<size_t>(state)]);
}

/******************************************************************************/

void TwoTreeLikelihood::initParameters()
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

void TwoTreeLikelihood::applyParameters() throw (Exception)
{
  // Apply branch length:
  brLen_ = getParameterValue("BrLen");
  // Apply substitution model parameters:
  model_->matchParametersValues(getParameters());
  // Apply rate distribution parameters:
  rateDistribution_->matchParametersValues(getParameters());
}

/******************************************************************************/

void TwoTreeLikelihood::initBranchLengthsParameters()
{
  if (brLen_ < minimumBrLen_)
  {
   ApplicationTools::displayWarning("Branch length is too small: " + TextTools::toString(brLen_) + ". Value is set to " + TextTools::toString(minimumBrLen_));
    brLen_ = minimumBrLen_;
  }
  brLenParameters_.reset();
  brLenParameters_.addParameter(Parameter("BrLen", brLen_, brLenConstraint_));
}

/******************************************************************************/

void TwoTreeLikelihood::setParameters(const ParameterList& parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

void TwoTreeLikelihood::fireParameterChanged(const ParameterList& params)
{
  applyParameters();

  // For now we ignore the parameter that changed and we recompute all arrays...

  // Computes all pxy and pyx once for all:
  pxy_.resize(nbClasses_);
  for (size_t c = 0; c < nbClasses_; c++)
  {
    VVdouble* pxy_c = &pxy_[c];
    pxy_c->resize(nbStates_);
    RowMatrix<double> Q = model_->getPij_t(brLen_ * rateDistribution_->getCategory(c));
    for (size_t x = 0; x < nbStates_; x++)
    {
      Vdouble* pxy_c_x = &(*pxy_c)[x];
      pxy_c_x->resize(nbStates_);
      for (size_t y = 0; y < nbStates_; y++)
      {
        (*pxy_c_x)[y] = Q(x, y);
      }
    }
  }

  if (computeFirstOrderDerivatives_)
  {
    // Computes all dpxy/dt once for all:
    dpxy_.resize(nbClasses_);
    for (size_t c = 0; c < nbClasses_; c++)
    {
      VVdouble* dpxy_c = &dpxy_[c];
      dpxy_c->resize(nbStates_);
      double rc = rateDistribution_->getCategory(c);
      RowMatrix<double> dQ = model_->getdPij_dt(brLen_ * rc);
      for (size_t x = 0; x < nbStates_; x++)
      {
        Vdouble* dpxy_c_x = &(*dpxy_c)[x];
        dpxy_c_x->resize(nbStates_);
        for (size_t y = 0; y < nbStates_; y++)
        {
          (*dpxy_c_x)[y] = rc * dQ(x, y);
        }
      }
    }
  }

  if (computeSecondOrderDerivatives_)
  {
    // Computes all d2pxy/dt2 once for all:
    d2pxy_.resize(nbClasses_);
    for (size_t c = 0; c < nbClasses_; c++)
    {
      VVdouble* d2pxy_c = &d2pxy_[c];
      d2pxy_c->resize(nbStates_);
      double rc = rateDistribution_->getCategory(c);
      RowMatrix<double> d2Q = model_->getd2Pij_dt2(brLen_ * rc);
      for (size_t x = 0; x < nbStates_; x++)
      {
        Vdouble* d2pxy_c_x = &(*d2pxy_c)[x];
        d2pxy_c_x->resize(nbStates_);
        for (size_t y = 0; y < nbStates_; y++)
        {
          (*d2pxy_c_x)[y] = rc * rc * d2Q(x, y);
        }
      }
    }
  }

  computeTreeLikelihood();
  if (computeFirstOrderDerivatives_)
  {
    computeTreeDLikelihood();
  }
  if (computeSecondOrderDerivatives_)
  {
    computeTreeD2Likelihood();
  }
}

/******************************************************************************/

double TwoTreeLikelihood::getValue() const
throw (Exception)
{
  return -getLogLikelihood();
}

/******************************************************************************/

void TwoTreeLikelihood::initTreeLikelihoods(const SequenceContainer& sequences) throw (Exception)
{
  const Sequence* seq1 = &sequences.getSequence(seqnames_[0]);
  const Sequence* seq2 = &sequences.getSequence(seqnames_[1]);
  leafLikelihoods1_.resize(nbDistinctSites_);
  leafLikelihoods2_.resize(nbDistinctSites_);
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
   Vdouble* leafLikelihoods1_i = &leafLikelihoods1_[i];
   Vdouble* leafLikelihoods2_i = &leafLikelihoods2_[i];
   leafLikelihoods1_i->resize(nbStates_);
   leafLikelihoods2_i->resize(nbStates_);
   int state1 = seq1->getValue(i);
   int state2 = seq2->getValue(i);
    for (size_t s = 0; s < nbStates_; s++)
    {
      // Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
      // otherwise value set to 0:
      try
      {
        (*leafLikelihoods1_i)[s] = model_->getInitValue(s, state1);
        (*leafLikelihoods2_i)[s] = model_->getInitValue(s, state2);
      }
      catch (SequenceNotFoundException& snfe)
      {
        throw SequenceNotFoundException("TwoTreeLikelihood::initTreelikelihoods. Leaf name in tree not found in site conainer: ", snfe.getSequenceId());
      }
    }
  }

  // Initialize likelihood vector:
  rootLikelihoods_.resize(nbDistinctSites_);
  rootLikelihoodsS_.resize(nbDistinctSites_);
  rootLikelihoodsSR_.resize(nbDistinctSites_);
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    VVdouble* rootLikelihoods_i = &rootLikelihoods_[i];
    Vdouble* rootLikelihoodsS_i = &rootLikelihoodsS_[i];
    rootLikelihoods_i->resize(nbClasses_);
    rootLikelihoodsS_i->resize(nbClasses_);
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* rootLikelihoods_i_c = &(*rootLikelihoods_i)[c];
      rootLikelihoods_i_c->resize(nbStates_);
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*rootLikelihoods_i_c)[s] = 1.; // All likelihoods are initialized to 1.
      }
    }
  }

  // Initialize d and d2 likelihoods:
  dLikelihoods_.resize(nbDistinctSites_);
  d2Likelihoods_.resize(nbDistinctSites_);
}

/******************************************************************************/

void TwoTreeLikelihood::computeTreeLikelihood()
{
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    VVdouble* rootLikelihoods_i = &rootLikelihoods_[i];
    Vdouble* leafLikelihoods1_i = &leafLikelihoods1_[i];
    Vdouble* leafLikelihoods2_i = &leafLikelihoods2_[i];
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* rootLikelihoods_i_c = &(*rootLikelihoods_i)[c];
      VVdouble* pxy_c = &pxy_[c];
      for (size_t x = 0; x < nbStates_; x++)
      {
        Vdouble* pxy_c_x = &(*pxy_c)[x];
        double l = 0;
        double l1 = (*leafLikelihoods1_i)[x];
        for (size_t y = 0; y < nbStates_; y++)
        {
          double l2 = (*leafLikelihoods2_i)[y];
          l += l1 * l2 * (*pxy_c_x)[y];
        }
        (*rootLikelihoods_i_c)[x] = l;
      }
    }
  }

  Vdouble fr = model_->getFrequencies();
  Vdouble p = rateDistribution_->getProbabilities();
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    // For each site in the sequence,
    VVdouble* rootLikelihoods_i = &rootLikelihoods_[i];
    Vdouble* rootLikelihoodsS_i = &rootLikelihoodsS_[i];
    rootLikelihoodsSR_[i] = 0;
    for (size_t c = 0; c < nbClasses_; c++)
    {
      (*rootLikelihoodsS_i)[c] = 0;
      // For each rate classe,
      Vdouble* rootLikelihoods_i_c = &(*rootLikelihoods_i)[c];
      for (size_t x = 0; x < nbStates_; x++)
      {
        // For each initial state,
        (*rootLikelihoodsS_i)[c] += fr[x] * (*rootLikelihoods_i_c)[x];
      }
      rootLikelihoodsSR_[i] += p[c] * (*rootLikelihoodsS_i)[c];
    }
  }
}

/******************************************************************************/

void TwoTreeLikelihood::computeTreeDLikelihood()
{
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    Vdouble* leafLikelihoods1_i = &leafLikelihoods1_[i];
    Vdouble* leafLikelihoods2_i = &leafLikelihoods2_[i];
    double dli = 0;
    for (size_t c = 0; c < nbClasses_; c++)
    {
      VVdouble* dpxy_c = &dpxy_[c];
      double dlic = 0;
      for (size_t x = 0; x < nbStates_; x++)
      {
        Vdouble* dpxy_c_x = &(*dpxy_c)[x];
        double l1 = (*leafLikelihoods1_i)[x];
        double dlicx = 0;
        for (size_t y = 0; y < nbStates_; y++)
        {
          double l2 = (*leafLikelihoods2_i)[y];
          dlicx += l1 * l2 * (*dpxy_c_x)[y];
        }
        dlic += dlicx * model_->freq(x);
      }
      dli += dlic * rateDistribution_->getProbability(c);
    }
    dLikelihoods_[i] = dli / rootLikelihoodsSR_[i];
  }
}

/******************************************************************************/

void TwoTreeLikelihood::computeTreeD2Likelihood()
{
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    Vdouble* leafLikelihoods1_i = &leafLikelihoods1_[i];
    Vdouble* leafLikelihoods2_i = &leafLikelihoods2_[i];
    double d2li = 0;
    for (size_t c = 0; c < nbClasses_; c++)
    {
      VVdouble* d2pxy_c = &d2pxy_[c];
      double d2lic = 0;
      for (size_t x = 0; x < nbStates_; x++)
      {
        Vdouble* d2pxy_c_x = &(*d2pxy_c)[x];
        double l1 = (*leafLikelihoods1_i)[x];
        double d2licx = 0;
        for (size_t y = 0; y < nbStates_; y++)
        {
          double l2 = (*leafLikelihoods2_i)[y];
          d2licx += l1 * l2 * (*d2pxy_c_x)[y];
        }
        d2lic += d2licx * model_->freq(x);
      }
      d2li += d2lic * rateDistribution_->getProbability(c);
    }
    d2Likelihoods_[i] = d2li / rootLikelihoodsSR_[i];
  }
}

/******************************************************************************/

double TwoTreeLikelihood::getFirstOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("TwoTreeLikelihood::getFirstOrderDerivative().", variable);
  if (getRateDistributionParameters().hasParameter(variable))
  {
    cout << "DEBUB: WARNING!!! Derivatives respective to rate distribution parameter are not implemented." << endl;
    return log(-1.);
  }
  if (getSubstitutionModelParameters().hasParameter(variable))
  {
    cout << "DEBUB: WARNING!!! Derivatives respective to substitution model parameters are not implemented." << endl;
    return log(-1.);
  }

  //
  // Computation for branch lengths:
  //

  // Get the node with the branch whose length must be derivated:
  double d = 0;
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    d += rootWeights_[i] * dLikelihoods_[i];
  }
  return -d;
}

/******************************************************************************/

double TwoTreeLikelihood::getSecondOrderDerivative(const string& variable) const
throw (Exception)
{
  if (!hasParameter(variable))
    throw ParameterNotFoundException("TwoTreeLikelihood::getSecondOrderDerivative().", variable);
  if (getRateDistributionParameters().hasParameter(variable))
  {
    cout << "DEBUB: WARNING!!! Derivatives respective to rate distribution parameter are not implemented." << endl;
    return log(-1.);
  }
  if (getSubstitutionModelParameters().hasParameter(variable))
  {
    cout << "DEBUB: WARNING!!! Derivatives respective to substitution model parameters are not implemented." << endl;
    return log(-1.);
  }

  //
  // Computation for branch lengths:
  //

  // Get the node with the branch whose length must be derivated:
  double d2 = 0;
  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    d2 += rootWeights_[i] * (d2Likelihoods_[i] - pow(dLikelihoods_[i], 2));
  }
  return -d2;
}

/******************************************************************************/

void DistanceEstimation::computeMatrix() throw (NullPointerException)
{
  size_t n = sites_->getNumberOfSequences();
  vector<string> names = sites_->getSequencesNames();
  if (dist_ != 0) delete dist_;
  dist_ = new DistanceMatrix(names);
  optimizer_->setVerbose(static_cast<unsigned int>(max(static_cast<int>(verbose_) - 2, 0)));
  for (size_t i = 0; i < n; ++i)
  {
    (*dist_)(i, i) = 0;
    if (verbose_ == 1)
    {
      ApplicationTools::displayGauge(i, n - 1, '=');
    }
    for (size_t j = i + 1; j < n; j++)
    {
      if (verbose_ > 1)
      {
        ApplicationTools::displayGauge(j - i - 1, n - i - 2, '=');
      }
      TwoTreeLikelihood* lik =
        new TwoTreeLikelihood(names[i], names[j], *sites_, model_.get(), rateDist_.get(), verbose_ > 3);
      lik->initialize();
      lik->enableDerivatives(true);
      size_t d = SymbolListTools::getNumberOfDistinctPositions(sites_->getSequence(i), sites_->getSequence(j));
      size_t g = SymbolListTools::getNumberOfPositionsWithoutGap(sites_->getSequence(i), sites_->getSequence(j));
      lik->setParameterValue("BrLen", g == 0 ? lik->getMinimumBranchLength() : std::max(lik->getMinimumBranchLength(), static_cast<double>(d) / static_cast<double>(g)));
      // Optimization:
      optimizer_->setFunction(lik);
      optimizer_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
      ParameterList params = lik->getBranchLengthsParameters();
      params.addParameters(parameters_);
      optimizer_->init(params);
      optimizer_->optimize();
      // Store results:
      (*dist_)(i, j) = (*dist_)(j, i) = lik->getParameterValue("BrLen");
      delete lik;
    }
    if (verbose_ > 1 && ApplicationTools::message) ApplicationTools::message->endLine();
  }
}

/******************************************************************************/

