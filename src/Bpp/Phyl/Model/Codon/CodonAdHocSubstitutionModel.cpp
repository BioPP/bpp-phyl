// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "CodonAdHocSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel(
    shared_ptr<const GeneticCode> gCode,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    vector<unique_ptr<CoreCodonSubstitutionModelInterface>>& vpmodel,
    const string& name) :
  AbstractParameterAliasable(name + "."),
  AbstractCodonSubstitutionModel(gCode, std::move(pmod), name + "."),
  vModel_(),
  name_(name),
  freqSet_()
{
  for (auto& model : vpmodel)
  {
    if (model)
    {
      model->setNamespace(name + ".");
      if (model->hasCodonFrequencySet())
      {
        if (freqSet_)
          throw Exception("CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel : two sub models with FrequencySet");

        freqSet_ = unique_ptr<CodonFrequencySetInterface>(model->codonFrequencySet().clone());
      }
      addParameters_(model->getParameters());
      vModel_.push_back(std::move(model));
    }
  }
  computeFrequencies(true);
  updateMatrices_();
}

CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel(
    shared_ptr<const GeneticCode> gCode,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    vector<unique_ptr<CoreCodonSubstitutionModelInterface>>& vpmodel,
    const std::string& name) :
  AbstractParameterAliasable(name + "."),
  AbstractCodonSubstitutionModel(gCode, std::move(pmod1), std::move(pmod2), std::move(pmod3), name + "."),
  vModel_(),
  name_(name),
  freqSet_()
{
  for (auto& model : vpmodel)
  {
    if (model != NULL)
    {
      model->setNamespace(name + ".");
      if (model->hasCodonFrequencySet())
      {
        if (freqSet_)
          throw Exception("CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel : two sub models with FrequencySet");

        freqSet_ = unique_ptr<CodonFrequencySetInterface>(model->codonFrequencySet().clone());
      }

      addParameters_(model->getParameters());
      vModel_.push_back(std::move(model));
    }
  }

  computeFrequencies(true);
  updateMatrices_();
}

CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel(const CodonAdHocSubstitutionModel& model) :
  AbstractParameterAliasable(model),
  AbstractCodonSubstitutionModel(model),
  vModel_(),
  name_(model.name_),
  freqSet_()
{
  for (auto& mod : model.vModel_)
  {
    vModel_.emplace_back(mod->clone());
  }

  for (auto& mod : vModel_)
  {
    if (mod->hasCodonFrequencySet())
    {
      if (freqSet_)
        throw Exception("CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel : two sub models with FrequencySet");

      freqSet_ = unique_ptr<CodonFrequencySetInterface>(mod->codonFrequencySet().clone());
    }
  }
}

CodonAdHocSubstitutionModel& CodonAdHocSubstitutionModel::operator=(const CodonAdHocSubstitutionModel& model)
{
  AbstractParameterAliasable::operator=(model);
  AbstractCodonSubstitutionModel::operator=(model);
  name_ = model.name_;

  vModel_.clear();
  freqSet_ = 0;

  for (auto& mod : model.vModel_)
  {
    vModel_.emplace_back(mod->clone());
  }

  for (auto& mod : vModel_)
  {
    if (mod->hasCodonFrequencySet())
    {
      if (freqSet_)
        throw Exception("CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel : two sub models with FrequencySet");

      freqSet_ = unique_ptr<CodonFrequencySetInterface>(mod->codonFrequencySet().clone());
    }
  }

  return *this;
}

void CodonAdHocSubstitutionModel::setFreq(std::map<int, double>& frequencies)
{
  for (auto& model : vModel_)
  {
    model->setFreq(frequencies);
    matchParametersValues(model->getParameters());
  }
}

void CodonAdHocSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  for (auto& model : vModel_)
  {
    model->matchParametersValues(parameters);
  }

  // Beware: must be called at the end
  AbstractCodonSubstitutionModel::fireParameterChanged(parameters);
}

double CodonAdHocSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  double x(1);

  for (auto& model : vModel_)
  {
    x *= model->getCodonsMulRate(i, j);
  }

  return x;
}
