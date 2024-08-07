// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "CodonDistancePhaseFrequenciesSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

CodonDistancePhaseFrequenciesSubstitutionModel::CodonDistancePhaseFrequenciesSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    std::unique_ptr<CodonFrequencySetInterface> pfreq,
    std::shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("CodonDistPhasFreq."),
  AbstractCodonSubstitutionModel(gCode, std::move(pmod), "CodonDistPhasFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "CodonDistPhasFreq."),
  AbstractCodonPhaseFrequenciesSubstitutionModel(std::move(pfreq), "CodonDistPhasFreq.")
{
  computeFrequencies(true); // for init
  updateMatrices_();
  computeFrequencies(false);
}

CodonDistancePhaseFrequenciesSubstitutionModel::CodonDistancePhaseFrequenciesSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    std::unique_ptr<CodonFrequencySetInterface> pfreq,
    std::shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("CodonDistPhasFreq."),
  AbstractCodonSubstitutionModel(gCode, std::move(pmod1), std::move(pmod2), std::move(pmod3), "CodonDistPhasFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "CodonDistPhasFreq."),
  AbstractCodonPhaseFrequenciesSubstitutionModel(std::move(pfreq), "CodonDistPhasFreq.")
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
