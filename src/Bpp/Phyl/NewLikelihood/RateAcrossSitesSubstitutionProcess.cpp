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

#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace bpp;
using namespace std;

RateAcrossSitesSubstitutionProcess::RateAcrossSitesSubstitutionProcess(
    SubstitutionModel* model,
    DiscreteDistribution* rdist,
    ParametrizableTree* tree,
    bool checkRooted) :
  AbstractParameterAliasable(model ? model->getNamespace() : ""),
  AbstractSubstitutionProcess(tree, rdist ? rdist->getNumberOfCategories() : 0, model ? model->getNamespace() : ""),
  model_(model),
  rDist_(rdist),
  computingTree_()
{
  if (!model)
    throw Exception("RateAcrossSitesSubstitutionProcess. A model instance must be provided.");
  if (!rdist)
    throw Exception("RateAcrossSitesSubstitutionProcess. A rate distribution instance must be provided.");
  if (checkRooted && pTree_->isRooted())
    throw Exception("RateAcrossSitesSubstitutionProcess (constructor). Tree is rooted.");

  computingTree_.reset(new ComputingTree(*pTree_.get(), *rDist_.get()));
  computingTree_->addModel(model_.get());
  computingTree_->checkModelOnEachNode();  

  // Add parameters:
  addParameters_(tree->getParameters());  //Branch lengths
  addParameters_(model->getIndependentParameters()); //Substitution model
  addParameters_(rdist->getIndependentParameters()); //Rate distribution
}
    
RateAcrossSitesSubstitutionProcess::RateAcrossSitesSubstitutionProcess(const RateAcrossSitesSubstitutionProcess& rassp) :
  AbstractParameterAliasable(rassp),
  AbstractSubstitutionProcess(rassp),
  model_(rassp.model_->clone()),
  rDist_(rassp.rDist_->clone()),
  computingTree_()
{
  computingTree_.reset(new ComputingTree(*pTree_.get(), *rDist_.get()));
  computingTree_->addModel(model_.get());
  computingTree_->checkModelOnEachNode();  
}


RateAcrossSitesSubstitutionProcess& RateAcrossSitesSubstitutionProcess::operator=(const RateAcrossSitesSubstitutionProcess& rassp)
{
  AbstractParameterAliasable::operator=(rassp);
  AbstractSubstitutionProcess::operator=(rassp);
  model_.reset(rassp.model_->clone());
  rDist_.reset(rassp.rDist_->clone());

  computingTree_.reset(new ComputingTree(*pTree_.get(), *rDist_.get()));
  computingTree_->addModel(model_.get());
  computingTree_->checkModelOnEachNode();  
  
  return *this;
}

void RateAcrossSitesSubstitutionProcess::fireParameterChanged(const ParameterList& pl)
{
  //Update substitution model:
  if (model_->matchParametersValues(pl))
    computingTree_->updateAll();
  
  //Update rate distribution:
  rDist_->matchParametersValues(pl);  

  //Transition probabilities have changed and need to be recomputed:
  AbstractSubstitutionProcess::fireParameterChanged(pl);
}

