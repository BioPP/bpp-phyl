// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "CodonDistanceFrequenciesSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

CodonDistanceFrequenciesSubstitutionModel::CodonDistanceFrequenciesSubstitutionModel(
    shared_ptr<const GeneticCode> gCode,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    unique_ptr<CodonFrequencySetInterface> pfreq,
    shared_ptr<const AlphabetIndex2> pdist,
    bool paramSynRate) :
  AbstractParameterAliasable("CodonDistFreq."),
  AbstractCodonSubstitutionModel(gCode, std::move(pmod), "CodonDistFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "CodonDistFreq.", paramSynRate),
  AbstractCodonFrequenciesSubstitutionModel(std::move(pfreq), "CodonDistFreq.")
{
  computeFrequencies(true); // for initialization
  updateMatrices_();
  computeFrequencies(false);
}

CodonDistanceFrequenciesSubstitutionModel::CodonDistanceFrequenciesSubstitutionModel(
    shared_ptr<const GeneticCode> gCode,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    unique_ptr<CodonFrequencySetInterface> pfreq,
    shared_ptr<const AlphabetIndex2> pdist,
    bool paramSynRate) :
  AbstractParameterAliasable("CodonDistFreq."),
  AbstractCodonSubstitutionModel(gCode, std::move(pmod1), std::move(pmod2), std::move(pmod3), "CodonDistFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "CodonDistFreq.", paramSynRate),
  AbstractCodonFrequenciesSubstitutionModel(std::move(pfreq), "CodonDistFreq.")
{
  computeFrequencies(true); // for initialization
  updateMatrices_();
  computeFrequencies(false);
}

std::string CodonDistanceFrequenciesSubstitutionModel::getName() const
{
  return "CodonDistFreq";
}

void CodonDistanceFrequenciesSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);
  AbstractCodonFrequenciesSubstitutionModel::fireParameterChanged(parameters);
  getFrequencies_() = AbstractCodonFrequenciesSubstitutionModel::codonFrequencySet().getFrequencies();

  // Beware: must be call at the end
  AbstractCodonSubstitutionModel::fireParameterChanged(parameters);
}

double CodonDistanceFrequenciesSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i, j)
         * AbstractCodonFrequenciesSubstitutionModel::getCodonsMulRate(i, j);
}

void CodonDistanceFrequenciesSubstitutionModel::setNamespace(const std::string& st)
{
  AbstractParameterAliasable::setNamespace(st);
  AbstractCodonSubstitutionModel::setNamespace(st);
  AbstractCodonDistanceSubstitutionModel::setNamespace(st);
  AbstractCodonFrequenciesSubstitutionModel::setNamespace(st);
}

void CodonDistanceFrequenciesSubstitutionModel::setFreq(map<int, double>& frequencies)
{
  AbstractCodonFrequenciesSubstitutionModel::setFreq(frequencies);
  getFrequencies_() = AbstractCodonFrequenciesSubstitutionModel::codonFrequencySet().getFrequencies();
}
