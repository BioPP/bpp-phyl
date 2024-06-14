// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

# include "AbstractCodonFitnessSubstitutionModel.h"
using namespace bpp;
using namespace std;

/****************************************************************************************/

AbstractCodonFitnessSubstitutionModel::AbstractCodonFitnessSubstitutionModel(
    unique_ptr<FrequencySetInterface> pfitset,
    shared_ptr<const GeneticCode> pgencode,
    const string& prefix) :
  AbstractParameterAliasable(prefix),
  pfitset_(std::move(pfitset)),
  pgencode_(pgencode),
  fitName_("")
{
  if (!dynamic_cast<CodonFrequencySetInterface*>(pfitset_.get()))
    throw Exception("Bad type for fitness parameters" + pfitset_->getName());
  fitName_ = "fit_" + pfitset_->getNamespace();
  pfitset_->setNamespace(prefix + fitName_);
  addParameters_(pfitset_->getParameters());
}
/****************************************************************************************/

AbstractCodonFitnessSubstitutionModel::~AbstractCodonFitnessSubstitutionModel() {}

/****************************************************************************************/

void AbstractCodonFitnessSubstitutionModel::fireParameterChanged (const ParameterList& parameters)
{
  pfitset_->matchParametersValues(parameters);
}

/****************************************************************************************/

void AbstractCodonFitnessSubstitutionModel::setFreq(map<int, double>& frequencies)
{
  pfitset_->setFrequenciesFromAlphabetStatesFrequencies(frequencies);
  matchParametersValues(pfitset_->getParameters() );
}

/****************************************************************************************/

double AbstractCodonFitnessSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  double mu;

  double phi_j = pfitset_->getFrequencies() [j];
  double phi_i = pfitset_->getFrequencies() [i];

  if (phi_i == phi_j)
    mu = 1;
  else if (phi_i == 0)
    mu = 100;
  else if (phi_j == 0)
    mu = 0;
  else
    mu = -(log(phi_i / phi_j) / (1 - (phi_i / phi_j)));
  return mu;
}

/****************************************************************************************/
