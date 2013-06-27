//
// File: RateAcrossSitesSubstitutionProcess.cpp
// Created by: Julien Dutheil
// Created on: Tue May 15 13:11 2012
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

#include "RateAcrossSitesSubstitutionProcess.h"

using namespace bpp;

RateAcrossSitesSubstitutionProcess::RateAcrossSitesSubstitutionProcess(SubstitutionModel* model, DiscreteDistribution* rdist, ParametrizableTree* tree) :
  AbstractParameterAliasable(model ? model->getNamespace() : ""),
  AbstractSubstitutionProcess(tree, rdist ? rdist->getNumberOfCategories() : 0),
  model_(model),
  rDist_(rdist)
{
  if (!model)
    throw Exception("RateAcrossSitesSubstitutionProcess. A model instance must be provided.");
  if (!rdist)
    throw Exception("RateAcrossSitesSubstitutionProcess. A rate distribution instance must be provided.");

  // Add parameters:
  addParameters_(tree->getParameters());  //Branch lengths
  addParameters_(model->getParameters()); //Substitution model
  addParameters_(rdist->getParameters()); //Rate distribution
}
    
RateAcrossSitesSubstitutionProcess::RateAcrossSitesSubstitutionProcess(const RateAcrossSitesSubstitutionProcess& rassp) :
  AbstractParameterAliasable(rassp),
  AbstractSubstitutionProcess(rassp),
  model_(rassp.model_->clone()),
  rDist_(rassp.rDist_->clone())
{}

RateAcrossSitesSubstitutionProcess& RateAcrossSitesSubstitutionProcess::operator=(const RateAcrossSitesSubstitutionProcess& rassp)
{
  AbstractParameterAliasable::operator=(rassp);
  AbstractSubstitutionProcess::operator=(rassp);
  model_.reset(rassp.model_->clone());
  rDist_.reset(rassp.rDist_->clone());
  return *this;
}

void RateAcrossSitesSubstitutionProcess::fireParameterChanged(const ParameterList& pl)
{
  //Forward parameters:
  pTree_->matchParametersValues(pl);
  //Update substitution model:
  model_->matchParametersValues(pl);
  //Update rate distribution:
  rDist_->matchParametersValues(pl);
  //Transition probabilities have changed and need to be recomputed:
  AbstractSubstitutionProcess::fireParameterChanged(pl);
}

const Matrix<double>& RateAcrossSitesSubstitutionProcess::getTransitionProbabilities(int nodeId, size_t classIndex) const
{
  size_t i = getModelIndex_(nodeId, classIndex);
  if (!computeProbability_[i]) {
    computeProbability_[i] = true; //We record that we did this computation.
    //The transition matrix was never computed before. We therefore have to compute it first:
    double l = pTree_->getBranchLengthParameter(nodeId).getValue();
    double r = rDist_->getCategory(classIndex);
    probabilities_[i] = model_->getPij_t(l * r);
  }
  return probabilities_[i];
}

const Matrix<double>& RateAcrossSitesSubstitutionProcess::getTransitionProbabilitiesD1(int nodeId, size_t classIndex) const
{
  size_t i = getModelIndex_(nodeId, classIndex);
  if (!computeProbabilityD1_[i]) {
    computeProbabilityD1_[i] = true; //We record that we did this computation.
    //The transition matrix was never computed before. We therefore have to compute it first:
    double l = pTree_->getBranchLengthParameter(nodeId).getValue();
    double r = rDist_->getCategory(classIndex);
    probabilitiesD1_[i] = model_->getdPij_dt(l * r);
  }
  return probabilitiesD1_[i];
}

const Matrix<double>& RateAcrossSitesSubstitutionProcess::getTransitionProbabilitiesD2(int nodeId, size_t classIndex) const
{
  size_t i = getModelIndex_(nodeId, classIndex);
  if (!computeProbabilityD2_[i]) {
    computeProbabilityD2_[i] = true; //We record that we did this computation.
    //The transition matrix was never computed before. We therefore have to compute it first:
    double l = pTree_->getBranchLengthParameter(nodeId).getValue();
    double r = rDist_->getCategory(classIndex);
    probabilitiesD2_[i] = model_->getd2Pij_dt2(l * r);
  }
  return probabilitiesD2_[i];
}


