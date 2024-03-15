// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Seq/Alphabet/AlphabetTools.h>

#include "AbstractCodonAAFitnessSubstitutionModel.h"

using namespace bpp;
using namespace std;

/****************************************************************************************/

AbstractCodonAAFitnessSubstitutionModel::AbstractCodonAAFitnessSubstitutionModel(
    shared_ptr<FrequencySetInterface> pfitset,
    shared_ptr<const GeneticCode> pgencode,
    const string& prefix) :
  AbstractParameterAliasable(prefix), pfitset_(pfitset), pgencode_(pgencode), fitName_(""), stateMap_(new CanonicalStateMap(pgencode->getSourceAlphabet(), false)),
  protStateMap_(pfitset->getStateMap()), Ns_(1)
{
  if (!AlphabetTools::isProteicAlphabet(pfitset_->getAlphabet().get()))
    throw Exception("AbstractCodonAAFitnessSubstitutionModel::AbstractCodonAAFitnessSubstitutionModel need Proteic Fitness.");

  fitName_ = "fit_" + pfitset_->getNamespace();
  pfitset_->setNamespace(prefix + fitName_);

  addParameters_(pfitset_->getParameters());
}

AbstractCodonAAFitnessSubstitutionModel::~AbstractCodonAAFitnessSubstitutionModel()
{}

void AbstractCodonAAFitnessSubstitutionModel::fireParameterChanged (const ParameterList& parameters)
{
  AbstractParameterAliasable::fireParameterChanged(parameters);
  if (hasParameter("Ns"))
    Ns_ = getParameterValue("Ns");

  pfitset_->matchParametersValues(parameters);
}

void AbstractCodonAAFitnessSubstitutionModel::setFreq(map<int, double>& frequencies)
{
  pfitset_->setFrequenciesFromAlphabetStatesFrequencies(frequencies);
  matchParametersValues(pfitset_->getParameters() );
}

double AbstractCodonAAFitnessSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  double mu;

  int aai = pgencode_->translate(stateMap_->getAlphabetStateAsInt(i));
  int aaj = pgencode_->translate(stateMap_->getAlphabetStateAsInt(j));

  double phi_i = pfitset_->getFrequencies() [protStateMap_->getModelStates(aai)[0]];
  double phi_j = pfitset_->getFrequencies() [protStateMap_->getModelStates(aaj)[0]];

  if (phi_i == phi_j || Ns_ == 0)
    mu = 1;
  else if (phi_i == 0)
    mu = 100;
  else if (phi_j == 0)
    mu = 0;
  else
  {
    if (Ns_ == 1)
      mu = -(log(phi_i / phi_j) / (1 - (phi_i / phi_j)));
    else
    {
      double x = pow(phi_i / phi_j, Ns_);
      mu = -log(x) / (1 - x);
    }
  }

  return mu;
}
