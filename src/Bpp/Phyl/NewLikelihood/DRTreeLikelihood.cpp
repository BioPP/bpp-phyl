// File: DRTreeLikelihood.cpp
// Created by: Julien Dutheil
// Created on: Fri Apr 26 20:07 2013
// From file: DRNonHomogeneousTreeLikelihood.cpp
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

#include "DRTreeLikelihood.h"

using namespace std;
using namespace bpp;
using namespace newlik;

/******************************************************************************/

DRTreeLikelihood::DRTreeLikelihood(
  SubstitutionProcess* process,
  bool verbose)
throw (Exception) :
  AbstractTreeLikelihood(process),
  likelihoodData_(0),
  minusLogLik_(-1.)
{
  init_();
}

/******************************************************************************/

DRTreeLikelihood::DRTreeLikelihood(
  const SiteContainer& data,
  SubstitutionProcess* process,
  bool verbose)
throw (Exception) :
  AbstractTreeLikelihood(process),
  likelihoodData_(0),
  minusLogLik_(-1.)
{
  init_();
  setData(data);
}

/******************************************************************************/

void DRTreeLikelihood::init_() throw (Exception)
{
  likelihoodData_.reset(new DRTreeLikelihoodData(
    process_->getNumberOfClasses()));
}

/******************************************************************************/

DRTreeLikelihood::DRTreeLikelihood(const DRTreeLikelihood& lik) :
  AbstractTreeLikelihood(lik),
  likelihoodData_(0),
  minusLogLik_(lik.minusLogLik_)
{
  likelihoodData_.reset(lik.likelihoodData_->clone());
}

/******************************************************************************/

DRTreeLikelihood& DRTreeLikelihood::operator=(const DRTreeLikelihood& lik)
{
  AbstractTreeLikelihood::operator=(lik);
  likelihoodData_.reset(lik.likelihoodData_->clone());
  minusLogLik_ = lik.minusLogLik_;
  return *this;
}

/******************************************************************************/

void DRTreeLikelihood::setData(const SiteContainer& sites) throw (Exception)
{
  try {
    const TreeTemplate<Node>& tt = dynamic_cast<const TreeTemplate<Node>&>(process_->getTree());
    data_.reset(PatternTools::getSequenceSubset(sites, *tt.getRootNode()));
  } catch (exception& e) {
    throw Exception("DEBUG. DRTreeLikelihood::setData. The SubstitutionProcess does not use a TreeTemplate object.");
  }
  if (verbose_)
    ApplicationTools::displayTask("Initializing data structure");
  likelihoodData_->initLikelihoods(*data_, *process_); //We assume here that all models have the same number of states, and that they have the same 'init' method,
                                                       //Which is a reasonable assumption as long as they share the same alphabet.
  if (verbose_)
    ApplicationTools::displayTaskDone();

  nbSites_         = likelihoodData_->getNumberOfSites();
  nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();
  nbStates_        = likelihoodData_->getNumberOfStates();

  if (verbose_) ApplicationTools::displayResult("Number of distinct sites", TextTools::toString(nbDistinctSites_));
  
  initialized_ = true;

  //Recompute likelihood:
  computeTreeLikelihood();
  minusLogLik_ = -getLogLikelihood();
}

/******************************************************************************/

double DRTreeLikelihood::getLogLikelihood() const
{
  vector<double> la(nbSites_);
  for (size_t i = 0; i < nbSites_; ++i)
  {
    la[i] = log(getLikelihoodForASite(i));
  }
  sort(la.begin(), la.end());
  double ll = 0;
  for (size_t i = nbSites_; i > 0; i--)
  {
    ll += la[i - 1];
  }
  return ll;
}

/******************************************************************************/

double RTreeLikelihood::getLikelihoodForASite(size_t site) const
{
  double l = 0;
  for (size_t i = 0; i < nbClasses_; ++i)
  {
    l += getLikelihoodForASiteForAClass(site, i) * process_->getProbabilityForModel(i);
  }
  //if(l <= 0.) cerr << "WARNING!!! Negative likelihood." << endl;
  if (l < 0) l = 0; //May happen because of numerical errors.
  return l;
}

/******************************************************************************/

double RTreeLikelihood::getLikelihoodForASiteForAClass(size_t site, size_t classIndex) const
{
  double l = 0;
  Vdouble* la = &likelihoodData_->getLikelihoodArray(
      process_->getTree().getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][classIndex];
  for (size_t i = 0; i < nbStates_; ++i)
  {
    l += (*la)[i] * process_->getRootFrequencies()[i];
  }
  return l;
}

/******************************************************************************/

double RTreeLikelihood::getLikelihoodForASiteForAState(size_t site, int state) const
{
  double l = 0;
  VVdouble* la = &likelihoodData_->getLikelihoodArray(
      process_->getTree().getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)];
  for (size_t c = 0; c < nbClasses_; ++c)
  {
    l += (*la)[c][state] * process_->getProbabilityForModel(c);
  }
  return l;
}

/******************************************************************************/

double RTreeLikelihood::getLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state) const
{
  return likelihoodData_->getLikelihoodArray(
      process_->getTree().getRootNode()->getId())[likelihoodData_->getRootArrayPosition(site)][classIndex][state];
}

/******************************************************************************/

void RTreeLikelihood::setParameters(const ParameterList& parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

/******************************************************************************/

void RTreeLikelihood::fireParameterChanged(const ParameterList& params)
{
  //TODO: check if that is still needed (jdutheil 03/12/12)
  //applyParameters();

  if (params.size() > 0) {
    computeTreeLikelihood();
    minusLogLik_ = -getLogLikelihood();
  }
}

/******************************************************************************/

double RTreeLikelihood::getValue() const
throw (Exception)
{
  if (!isInitialized()) throw Exception("RTreeLikelihood::getValue(). Instance is not initialized.");
  return minusLogLik_;
}

/******************************************************************************
 *                           Likelihood computation                           *
 ******************************************************************************/

void RTreeLikelihood::computeTreeLikelihood()
{
  computeSubtreeLikelihood_(process_->getTree().getRootNode());
}

/******************************************************************************/


