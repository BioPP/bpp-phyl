// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_CODONDISTANCEPHASEFREQUENCIESSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_CODONDISTANCEPHASEFREQUENCIESSUBSTITUTIONMODEL_H


#include "AbstractCodonDistanceSubstitutionModel.h"
#include "AbstractCodonPhaseFrequenciesSubstitutionModel.h"
#include "AbstractCodonSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Class for asynonymous substitution models on codons with
 * parameterized equilibrium frequencies and nucleotidic basic models.
 *
 * @author Laurent GuÃÂ©guen
 *
 * This class should be used with models which equilibrium
 * distribution is fixed, and does not depend on the parameters.
 * Otherwise there may be problems of identifiability of the
 * parameters.
 *
 * See description in AbstractCodonDistanceSubstitutionModel and
 * AbstractCodonPhaseFrequenciesSubstitutionModel class.
 *
 * Only substitutions with one letter changed are accepted. </p>
 *
 * The additional parameter to CodonPhaseFrequenciesSubstitutionModel
 * is the ratio of nonsynonymous over synonymous substitutions.
 *
 * If a distance @f$d@f$ between amino-acids is defined, the ratio between
 * non-synonymous and synonymous substitutions rates is, if the codied
 * amino-acids are @f$x@f$ and @f$y@f$, @f$\beta*\exp(-\alpha.d(x,y))@f$ with
 * non-negative parameter \c "alpha" and positive parameter \c "beta".
 *
 * If such a distance is not defined, the ratio between non-synonymous
 * and synonymous substitutions rates is @f$\beta@f$ with positive
 * parameter \c "beta".
 */
class CodonDistancePhaseFrequenciesSubstitutionModel :
  public AbstractCodonSubstitutionModel,
  public AbstractCodonDistanceSubstitutionModel,
  public AbstractCodonPhaseFrequenciesSubstitutionModel
{
public:
  /**
   * @brief Build a new CodonDistancePhaseFrequenciesSubstitutionModel object
   * from three pointers to AbstractSubstitutionModels. NEW
   * AbstractSubstitutionModels are copied from the given ones.
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod pointer to the NucleotideSubstitutionModel to use in
   *        the three positions. It is owned by the instance.
   * @param pfreq pointer to the std::shared_ptr<FrequencySet> equilibrium frequencies
   * @param pdist optional pointer to the AlphabetIndex2 amino-acids distance object.
   */
  CodonDistancePhaseFrequenciesSubstitutionModel(
      std::shared_ptr<const GeneticCode> gCode,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
      std::unique_ptr<CodonFrequencySetInterface> pfreq,
      std::shared_ptr<const AlphabetIndex2> pdist = nullptr);

  /**
   * @brief Build a new CodonDistancePhaseFrequenciesSubstitutionModel object
   * from three pointers to AbstractSubstitutionModels. NEW
   * AbstractSubstitutionModels are copied from the given ones.
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod1, pmod2, pmod3 are pointers to the
   *   NucleotideSubstitutionModel to use in the three positions.
   *   All the models must be different objects to avoid redundant
   *   parameters.  They are owned by the instance.
   * @param pfreq pointer to the std::shared_ptr<FrequencySet> equilibrium frequencies
   * @param pdist optional pointer to the AlphabetIndex2 amino-acids distance object.
   */
  CodonDistancePhaseFrequenciesSubstitutionModel(
      std::shared_ptr<const GeneticCode> gCode,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
      std::unique_ptr<CodonFrequencySetInterface> pfreq,
      std::shared_ptr<const AlphabetIndex2> pdist = nullptr);

  virtual ~CodonDistancePhaseFrequenciesSubstitutionModel() {}

  CodonDistancePhaseFrequenciesSubstitutionModel* clone() const override
  {
    return new CodonDistancePhaseFrequenciesSubstitutionModel(*this);
  }

public:
  void fireParameterChanged(const ParameterList& parameterlist) override;

  std::string getName() const override;

  double getCodonsMulRate(size_t i, size_t j) const override;

  void setNamespace(const std::string&) override;

  void setFreq(std::map<int, double>& frequencies) override;

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    return AbstractCodonPhaseFrequenciesSubstitutionModel::codonFrequencySet();
  }

  const FrequencySetInterface& frequencySet() const override
  {
    return AbstractCodonPhaseFrequenciesSubstitutionModel::frequencySet();
  }

  bool hasCodonFrequencySet() const override
  {
    return AbstractCodonPhaseFrequenciesSubstitutionModel::hasCodonFrequencySet();
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_CODONDISTANCEPHASEFREQUENCIESSUBSTITUTIONMODEL_H
