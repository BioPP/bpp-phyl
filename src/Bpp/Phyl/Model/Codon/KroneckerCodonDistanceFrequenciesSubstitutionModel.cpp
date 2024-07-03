// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "KroneckerCodonDistanceFrequenciesSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

KroneckerCodonDistanceFrequenciesSubstitutionModel::KroneckerCodonDistanceFrequenciesSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    std::unique_ptr<CodonFrequencySetInterface> pfreq,
    std::shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("KronCodonDistFreq."),
  AbstractKroneckerWordSubstitutionModel(
      gCode->getSourceAlphabet(),
      shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
      "KronCodonDistFreq."),
  AbstractKroneckerCodonSubstitutionModel(gCode, std::move(pmod), "KronCodonDistFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "KronCodonDistFreq."),
  AbstractCodonFrequenciesSubstitutionModel(std::move(pfreq), "KronCodonDistFreq.")
{
  computeFrequencies(false);
  updateMatrices_();
}

KroneckerCodonDistanceFrequenciesSubstitutionModel::KroneckerCodonDistanceFrequenciesSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    const std::vector<std::set<size_t>>& vPos,
    std::unique_ptr<CodonFrequencySetInterface> pfreq,
    std::shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("KronCodonDistFreq."),
  AbstractKroneckerWordSubstitutionModel(
      gCode->getSourceAlphabet(),
      shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
      "KronCodonDistFreq."),
  AbstractKroneckerCodonSubstitutionModel(gCode, std::move(pmod), vPos, "KronCodonDistFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "KronCodonDistFreq."),
  AbstractCodonFrequenciesSubstitutionModel(std::move(pfreq), "KronCodonDistFreq.")
{
  computeFrequencies(false);
  updateMatrices_();
}

KroneckerCodonDistanceFrequenciesSubstitutionModel::KroneckerCodonDistanceFrequenciesSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    std::unique_ptr<CodonFrequencySetInterface> pfreq,
    std::shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("KronCodonDistFreq."),
  AbstractKroneckerWordSubstitutionModel(
      gCode->getSourceAlphabet(),
      shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
      "KronCodonDistFreq."),
  AbstractKroneckerCodonSubstitutionModel(gCode, std::move(pmod1), std::move(pmod2), std::move(pmod3), "KronCodonDistFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "KronCodonDistFreq."),
  AbstractCodonFrequenciesSubstitutionModel(std::move(pfreq), "KronCodonDistFreq.")
{
  computeFrequencies(false);
  updateMatrices_();
}

KroneckerCodonDistanceFrequenciesSubstitutionModel::KroneckerCodonDistanceFrequenciesSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    const std::vector<std::set< size_t>>& vPos,
    std::unique_ptr<CodonFrequencySetInterface> pfreq,
    std::shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("KronCodonDistFreq."),
  AbstractKroneckerWordSubstitutionModel(
      gCode->getSourceAlphabet(),
      shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
      "KronCodonDistFreq."),
  AbstractKroneckerCodonSubstitutionModel(gCode, std::move(pmod1), std::move(pmod2), std::move(pmod3), vPos, "KronCodonDistFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "KronCodonDistFreq."),
  AbstractCodonFrequenciesSubstitutionModel(std::move(pfreq), "KronCodonDistFreq.")
{
  computeFrequencies(false);
  updateMatrices_();
}

std::string KroneckerCodonDistanceFrequenciesSubstitutionModel::getName() const
{
  return "KronCodonDistFreq";
}

void KroneckerCodonDistanceFrequenciesSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);
  AbstractCodonFrequenciesSubstitutionModel::fireParameterChanged(parameters);

  // Beware: must be called last
  AbstractKroneckerCodonSubstitutionModel::fireParameterChanged(parameters);
}

double KroneckerCodonDistanceFrequenciesSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i, j)
         * AbstractKroneckerCodonSubstitutionModel::getCodonsMulRate(i, j)
         * AbstractCodonFrequenciesSubstitutionModel::getCodonsMulRate(i, j);
}

void KroneckerCodonDistanceFrequenciesSubstitutionModel::setNamespace(const std::string& st)
{
  AbstractParameterAliasable::setNamespace(st);
  AbstractKroneckerCodonSubstitutionModel::setNamespace(st);
  AbstractCodonDistanceSubstitutionModel::setNamespace(st);
  AbstractCodonFrequenciesSubstitutionModel::setNamespace(st);
}

void KroneckerCodonDistanceFrequenciesSubstitutionModel::setFreq(map<int, double>& frequencies)
{
  AbstractCodonFrequenciesSubstitutionModel::setFreq(frequencies);
}
