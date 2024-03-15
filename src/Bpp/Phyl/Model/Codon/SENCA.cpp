// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SENCA.h"

using namespace bpp;
using namespace std;

SENCA::SENCA(
  shared_ptr<const GeneticCode> gCode,
  unique_ptr<NucleotideSubstitutionModelInterface> pmod,
  unique_ptr<FrequencySetInterface> pfit,
  shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("SENCA."),
  AbstractCodonSubstitutionModel(gCode, std::move(pmod), "SENCA."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "SENCA."),
  AbstractCodonFitnessSubstitutionModel(std::move(pfit), gCode, "SENCA.")
{
  computeFrequencies(true);
  updateMatrices_();
}

SENCA::SENCA(
  shared_ptr<const GeneticCode> gCode,
  unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
  unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
  unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
  unique_ptr<FrequencySetInterface> pfit,
  shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("SENCA."),
  AbstractCodonSubstitutionModel(gCode, std::move(pmod1), std::move(pmod2), std::move(pmod3), "SENCA."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "SENCA."),
  AbstractCodonFitnessSubstitutionModel(std::move(pfit), gCode, "SENCA.")
{
  computeFrequencies(true);
  updateMatrices_();
}

void SENCA::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);
  AbstractCodonFitnessSubstitutionModel::fireParameterChanged(parameters);

  AbstractCodonSubstitutionModel::fireParameterChanged(parameters);
}

double SENCA::getCodonsMulRate(size_t i, size_t j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i, j)
         * AbstractCodonFitnessSubstitutionModel::getCodonsMulRate(i, j);
}

void SENCA::setNamespace(const std::string& st)
{
  AbstractParameterAliasable::setNamespace(st);
  AbstractCodonSubstitutionModel::setNamespace(st);
  AbstractCodonDistanceSubstitutionModel::setNamespace(st);
  AbstractCodonFitnessSubstitutionModel::setNamespace(st);
}

void SENCA::setFreq(map<int, double>& frequencies)
{
  AbstractSubstitutionModel::setFreq(frequencies);
  const Vdouble& freq1 = AbstractCodonSubstitutionModel::getFrequencies();
  auto alphabet = getAlphabet();

  map<int, double> freq2;
  double s = 0;

  for (auto& it : frequencies)
  {
    freq2[it.first] = (freq1[alphabet->getStateIndex(it.first) - 1] != 0 ? it.second / freq1[alphabet->getStateIndex(it.first) - 1] : 0);
    s += freq2[it.first];
  }

  for (auto& it : freq2)
  {
    freq2[it.first] /= s;
  }

  AbstractCodonFitnessSubstitutionModel::setFreq(freq2);

  updateMatrices_();
}
