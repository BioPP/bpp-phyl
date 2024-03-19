// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractCodonSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

AbstractCodonSubstitutionModel::AbstractCodonSubstitutionModel(
    shared_ptr<const GeneticCode> gCode,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    const string& prefix,
    bool paramRates) :
  AbstractParameterAliasable(prefix),
  AbstractWordSubstitutionModel(
    gCode->getSourceAlphabet(),
    std::shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
    prefix),
  hasParametrizedRates_(paramRates),
  gCode_(gCode)
{
  enableEigenDecomposition(true);

  shared_ptr<NucleotideSubstitutionModelInterface> pmod2 = std::move(pmod);
  for (size_t i = 0; i < 3; ++i)
  {
    VSubMod_.push_back(pmod2);
    VnestedPrefix_.push_back(pmod2->getNamespace());
  }

  pmod2->setNamespace(prefix + "123_" + VnestedPrefix_[0]);
  pmod2->enableEigenDecomposition(0);

  addParameters_(pmod2->getParameters());

  Vrate_.resize(3);
  for (size_t i = 0; i < 3; ++i)
  {
    Vrate_[i] = 1.0 / 3;
  }


  if (hasParametrizedRates_)
  {
    // relative rates
    for (size_t i = 0; i < 2; ++i)
    {
      addParameter_(new Parameter(prefix + "relrate" + TextTools::toString(i + 1), 1.0 / double(3 - i), Parameter::PROP_CONSTRAINT_EX));
    }
  }
}

AbstractCodonSubstitutionModel::AbstractCodonSubstitutionModel(
    shared_ptr<const GeneticCode> gCode,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    const std::string& prefix,
    bool paramRates) :
  AbstractParameterAliasable(prefix),
  AbstractWordSubstitutionModel(
    gCode->getSourceAlphabet(),
    shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
    prefix),
  hasParametrizedRates_(paramRates),
  gCode_(gCode)
{
  enableEigenDecomposition(1);

  VSubMod_.push_back(std::move(pmod1));
  VnestedPrefix_.push_back(VSubMod_[0]->getNamespace());
  VSubMod_[0]->setNamespace(prefix + "1_" + VnestedPrefix_[0]);
  VSubMod_[0]->enableEigenDecomposition(0);
  addParameters_(VSubMod_[0]->getParameters());

  VSubMod_.push_back(std::move(pmod2));
  VnestedPrefix_.push_back(VSubMod_[1]->getNamespace());
  VSubMod_[1]->setNamespace(prefix + "2_" + VnestedPrefix_[1]);
  VSubMod_[1]->enableEigenDecomposition(0);
  addParameters_(VSubMod_[1]->getParameters());

  VSubMod_.push_back(std::move(pmod3));
  VnestedPrefix_.push_back(VSubMod_[2]->getNamespace());
  VSubMod_[2]->setNamespace(prefix + "3_" + VnestedPrefix_[2]);
  VSubMod_[2]->enableEigenDecomposition(0);
  addParameters_(VSubMod_[2]->getParameters());

  Vrate_.resize(3);
  for (size_t i = 0; i < 3; ++i)
  {
    Vrate_[i] = 1.0 / 3;
  }

  if (hasParametrizedRates_)
  {
    // relative rates
    for (int i = 0; i < 2; ++i)
    {
      addParameter_(new Parameter(prefix + "relrate" + TextTools::toString(i + 1), 1.0 / double(3 - i), Parameter::PROP_CONSTRAINT_EX));
    }
  }
}

void AbstractCodonSubstitutionModel::updateMatrices_()
{
  if (hasParametrizedRates_)
  {
    size_t i, nbmod = VSubMod_.size();
    double x;
    size_t k;
    for (k = nbmod - 1; k > 0; k--)
    {
      x = 1.0;
      for (i = 0; i < k; i++)
      {
        x *= 1 - getParameterValue("relrate" + TextTools::toString(i + 1));
      }
      if (k != nbmod - 1)
        x *= getParameterValue("relrate" + TextTools::toString(k + 1));
      Vrate_[k] = x;
    }
  }

  AbstractWordSubstitutionModel::updateMatrices_();
}

void AbstractCodonSubstitutionModel::completeMatrices_()
{
  size_t i, j;
  size_t salph = getNumberOfStates();

  for (i = 0; i < salph; i++)
  {
    for (j = 0; j < salph; j++)
    {
      if (gCode_->isStop(getAlphabetStateAsInt(i)) || gCode_->isStop(getAlphabetStateAsInt(j)))
      {
        generator_(i, j) = 0;
      }
      else
        generator_(i, j) *= getCodonsMulRate(i, j);
    }
  }
}
