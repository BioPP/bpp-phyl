// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "KroneckerCodonDistanceSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

KroneckerCodonDistanceSubstitutionModel::KroneckerCodonDistanceSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    std::shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("KronCodonDist."),
  AbstractKroneckerWordSubstitutionModel(
      gCode->getSourceAlphabet(),
      shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
      "KronCodonDist."),
  AbstractKroneckerCodonSubstitutionModel(gCode, std::move(pmod), "KronCodonDist."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "KronCodonDist.")
{
  computeFrequencies(true);
  updateMatrices_();
}

KroneckerCodonDistanceSubstitutionModel::KroneckerCodonDistanceSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    std::shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("KronCodonDist."),
  AbstractKroneckerWordSubstitutionModel(
      gCode->getSourceAlphabet(),
      shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
      "KronCodonDist."),
  AbstractKroneckerCodonSubstitutionModel(gCode, std::move(pmod1), std::move(pmod2), std::move(pmod3), "KronCodonDist."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "KronCodonDist.")
{
  computeFrequencies(true);
  updateMatrices_();
}

KroneckerCodonDistanceSubstitutionModel::KroneckerCodonDistanceSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    const vector<std::set< size_t>>& vPos,
    std::shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("KronCodonDist."),
  AbstractKroneckerWordSubstitutionModel(
      gCode->getSourceAlphabet(),
      shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
      "KronCodonDist."),
  AbstractKroneckerCodonSubstitutionModel(gCode, std::move(pmod), vPos, "KronCodonDist."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "KronCodonDist.")
{
  computeFrequencies(true);
  updateMatrices_();
}

KroneckerCodonDistanceSubstitutionModel::KroneckerCodonDistanceSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    const vector<std::set<size_t>>& vPos,
    std::shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("KronCodonDist."),
  AbstractKroneckerWordSubstitutionModel(
      gCode->getSourceAlphabet(),
      shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
      "KronCodonDist."),
  AbstractKroneckerCodonSubstitutionModel(gCode, std::move(pmod1), std::move(pmod2), std::move(pmod3), vPos, "KronCodonDist."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "KronCodonDist.")
{
  computeFrequencies(true);
  updateMatrices_();
}

string KroneckerCodonDistanceSubstitutionModel::getName() const
{
  return "KronCodonDist";
}

void KroneckerCodonDistanceSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);

  // Beware: must be call at the end
  AbstractKroneckerCodonSubstitutionModel::fireParameterChanged(parameters);
}

double KroneckerCodonDistanceSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i, j)
         * AbstractKroneckerCodonSubstitutionModel::getCodonsMulRate(i, j);
}

void KroneckerCodonDistanceSubstitutionModel::setNamespace(const std::string& st)
{
  AbstractParameterAliasable::setNamespace(st);
  AbstractKroneckerCodonSubstitutionModel::setNamespace(st);
  AbstractCodonDistanceSubstitutionModel::setNamespace(st);
}
