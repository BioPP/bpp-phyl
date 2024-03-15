// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "CodonDistancePhaseFrequenciesSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

CodonDistancePhaseFrequenciesSubstitutionModel::CodonDistancePhaseFrequenciesSubstitutionModel(
    shared_ptr<const GeneticCode> gCode,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    unique_ptr<CodonFrequencySetInterface> pfreq,
    shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("CodonDistPhasFreq."),
  AbstractCodonSubstitutionModel(gCode, move(pmod), "CodonDistPhasFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "CodonDistPhasFreq."),
  AbstractCodonPhaseFrequenciesSubstitutionModel(move(pfreq), "CodonDistPhasFreq.")
{
  computeFrequencies(true); // for init
  updateMatrices_();
  computeFrequencies(false);
}

CodonDistancePhaseFrequenciesSubstitutionModel::CodonDistancePhaseFrequenciesSubstitutionModel(
    shared_ptr<const GeneticCode> gCode,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    unique_ptr<CodonFrequencySetInterface> pfreq,
    shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("CodonDistPhasFreq."),
  AbstractCodonSubstitutionModel(gCode, move(pmod1), move(pmod2), move(pmod3), "CodonDistPhasFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "CodonDistPhasFreq."),
  AbstractCodonPhaseFrequenciesSubstitutionModel(move(pfreq), "CodonDistPhasFreq.")
{
  computeFrequencies(true);
  updateMatrices_();
  computeFrequencies(false);
}

std::string CodonDistancePhaseFrequenciesSubstitutionModel::getName() const
{
  return "CodonDistPhasFreq";
}

void CodonDistancePhaseFrequenciesSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);
  AbstractCodonPhaseFrequenciesSubstitutionModel::fireParameterChanged(parameters);
  getFrequencies_() = AbstractCodonPhaseFrequenciesSubstitutionModel::codonFrequencySet().getFrequencies();

  // Beware: must be call at the end
  AbstractCodonSubstitutionModel::fireParameterChanged(parameters);
}

double CodonDistancePhaseFrequenciesSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i, j)
         * AbstractCodonPhaseFrequenciesSubstitutionModel::getCodonsMulRate(i, j);
}

void CodonDistancePhaseFrequenciesSubstitutionModel::setNamespace(const std::string& st)
{
  AbstractParameterAliasable::setNamespace(st);
  AbstractCodonSubstitutionModel::setNamespace(st);
  AbstractCodonDistanceSubstitutionModel::setNamespace(st);
  AbstractCodonPhaseFrequenciesSubstitutionModel::setNamespace(st);
}

void CodonDistancePhaseFrequenciesSubstitutionModel::setFreq(map<int, double>& frequencies)
{
  AbstractCodonPhaseFrequenciesSubstitutionModel::setFreq(frequencies);
  getFrequencies_() = AbstractCodonPhaseFrequenciesSubstitutionModel::codonFrequencySet().getFrequencies();
}
