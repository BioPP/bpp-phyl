// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "DFPDistanceFrequenciesSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

DFPDistanceFrequenciesSubstitutionModel::DFPDistanceFrequenciesSubstitutionModel(
    shared_ptr<const GeneticCode> gCode,
    unique_ptr<CodonFrequencySetInterface> pfreq,
    shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("DFPDistFreq."),
  AbstractDFPSubstitutionModel(gCode, "DFPDistFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "DFPDistFreq."),
  AbstractCodonFrequenciesSubstitutionModel(std::move(pfreq), "DFPDistFreq.")
{
  updateMatrices_();
}

std::string DFPDistanceFrequenciesSubstitutionModel::getName() const
{
  return "DFPDistFreq";
}

void DFPDistanceFrequenciesSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);
  AbstractCodonFrequenciesSubstitutionModel::fireParameterChanged(parameters);

  // Beware: must be called last, since it calls updateMatrices
  AbstractDFPSubstitutionModel::fireParameterChanged(parameters);
}

double DFPDistanceFrequenciesSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i, j)
         * AbstractDFPSubstitutionModel::getCodonsMulRate(i, j)
         * AbstractCodonFrequenciesSubstitutionModel::getCodonsMulRate(i, j);
}

void DFPDistanceFrequenciesSubstitutionModel::setNamespace(const std::string& st)
{
  AbstractParameterAliasable::setNamespace(st);
  AbstractDFPSubstitutionModel::setNamespace(st);
  AbstractCodonDistanceSubstitutionModel::setNamespace(st);
  AbstractCodonFrequenciesSubstitutionModel::setNamespace(st);
}

void DFPDistanceFrequenciesSubstitutionModel::setFreq(map<int, double>& frequencies)
{
  AbstractCodonFrequenciesSubstitutionModel::setFreq(frequencies);
  getFrequencies_() = AbstractCodonFrequenciesSubstitutionModel::codonFrequencySet().getFrequencies();
}
