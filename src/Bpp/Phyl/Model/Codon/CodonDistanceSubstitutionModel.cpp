// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "CodonDistanceSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

CodonDistanceSubstitutionModel::CodonDistanceSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    std::shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("CodonDist."),
  AbstractCodonSubstitutionModel(gCode, std::move(pmod), "CodonDist."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "CodonDist.")
{
  computeFrequencies(true);
  updateMatrices_();
}

CodonDistanceSubstitutionModel::CodonDistanceSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    std::shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("CodonDist."),
  AbstractCodonSubstitutionModel(gCode, std::move(pmod1), std::move(pmod2), std::move(pmod3), "CodonDist."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "CodonDist.")
{
  computeFrequencies(true);
  updateMatrices_();
}

std::string CodonDistanceSubstitutionModel::getName() const
{
  return "CodonDist";
}

void CodonDistanceSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);

  // Beware: must be call at the end
  AbstractCodonSubstitutionModel::fireParameterChanged(parameters);
}

double CodonDistanceSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i, j);
}
