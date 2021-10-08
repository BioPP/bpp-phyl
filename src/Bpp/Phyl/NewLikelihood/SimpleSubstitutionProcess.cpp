//
// File: SimpleSubstitutionProcess.cpp
// Authors:
//   Julien Dutheil
// Created: 2012-05-15 13:11:00
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/


#include "SimpleSubstitutionProcess.h"

using namespace bpp;
using namespace std;

SimpleSubstitutionProcess::SimpleSubstitutionProcess(std::shared_ptr<BranchModel> model, ParametrizablePhyloTree* tree) :
  AbstractParameterAliasable(""),
  AbstractSubstitutionProcess(tree, 1, model ? model->getNamespace() : ""),
  model_(model)
{
  if (!model)
    throw Exception("SimpleSubstitutionProcess. A model instance must be provided.");

  // Add parameters:
  addParameters_(tree->getParameters());  // Branch lengths
  addParameters_(model->getIndependentParameters()); // Substitution model
}

SimpleSubstitutionProcess::SimpleSubstitutionProcess(const SimpleSubstitutionProcess& ssp) :
  AbstractParameterAliasable(ssp),
  AbstractSubstitutionProcess(ssp),
  model_(ssp.model_->clone())
{
  if (modelScenario_)
    modelScenario_->changeModel(std::dynamic_pointer_cast<MixedTransitionModel>(ssp.model_), std::dynamic_pointer_cast<MixedTransitionModel>(model_));
}

SimpleSubstitutionProcess& SimpleSubstitutionProcess::operator=(const SimpleSubstitutionProcess& ssp)
{
  AbstractParameterAliasable::operator=(ssp);
  AbstractSubstitutionProcess::operator=(ssp);
  model_.reset(ssp.model_->clone());

  if (modelScenario_)
    modelScenario_->changeModel(std::dynamic_pointer_cast<MixedTransitionModel>(ssp.model_), std::dynamic_pointer_cast<MixedTransitionModel>(model_));

  return *this;
}

void SimpleSubstitutionProcess::setModelScenario(std::shared_ptr<ModelScenario> modelpath)
{
  auto vmod = modelpath->getModels();

  if (vmod.size() == 0) // as if no scenario
    return;

  if (vmod.size() != 1)
    throw Exception("SimpleSubstitutionProcess::setModelPath: model path must have exactly one model.");

  if (vmod[0] != model_)
    throw Exception("SimpleSubstitutionProcess::setModelPath: models are different " + vmod[0]->getName() + " != " + model_->getName());

  modelScenario_ = modelpath;
}

void SimpleSubstitutionProcess::fireParameterChanged(const ParameterList& pl)
{
  // Transition probabilities have changed and need to be recomputed:
  AbstractSubstitutionProcess::fireParameterChanged(pl);

  model_->matchParametersValues(pl);
}
