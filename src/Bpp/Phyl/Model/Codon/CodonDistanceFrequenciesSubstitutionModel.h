// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_CODONDISTANCEFREQUENCIESSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_CODONDISTANCEFREQUENCIESSUBSTITUTIONMODEL_H


#include "AbstractCodonDistanceSubstitutionModel.h"
#include "AbstractCodonFrequenciesSubstitutionModel.h"
#include "AbstractCodonSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Class for asynonymous substitution models on codons with
 * parameterized equilibrium frequencies and nucleotidic models.
 * @author Laurent GuÃÂ©guen
 *
 * This class should be used with models which equilibrium
 * distribution is fixed, and does not depend on the parameters.
 * Otherwise there may be problems of identifiability of the
 * parameters.
 *
 * See description in AbstractCodonDistanceSubstitutionModel and
 * AbstractCodonFrequenciesSubstitutionModel class.
 *
 * Only substitutions with one letter changed are accepted. </p>
 *
 * The additional parameters to CodonFrequenciesSubstitutionModel are
 * the rates of nonsynonymous over synonymous substitutions.
 *
 * If a distance @f$d@f$ between amino-acids is defined, the
 *  non-synonymous rate is multiplied with, if the coded amino-acids
 *  are @f$x@f$ and @f$y@f$, @f$\beta*\exp(-\alpha.d(x,y))@f$ with
 *  non-negative parameter \c "alpha" and positive parameter \c
 *  "beta".
 *
 * If such a distance is not defined, the non-synonymous substitution
 *  rate is multiplied with @f$\beta@f$ with positive parameter \c
 *  "beta" (ie @f$d=0@f$).
 *
 * If such a distance is not defined, the ratio between non-synonymous
 * and synonymous substitutions rates is @f$\beta@f$ with positive
 * parameter \c "beta".
 *
 * If paramSynRate is true, the synonymous substitution rate is
 *  multiplied with @f$\gamma@f$ (with optional positive parameter \c
 *  "gamma"), else it is multiplied with 1.
 *
 */

class CodonDistanceFrequenciesSubstitutionModel :
  public AbstractCodonSubstitutionModel,
  public AbstractCodonDistanceSubstitutionModel,
  public AbstractCodonFrequenciesSubstitutionModel
{
public:
  /**
   * @brief Build a new CodonDistanceFrequenciesSubstitutionModel object
   * from three pointers to AbstractSubstitutionModels. NEW
   * AbstractSubstitutionModels are copied from the given ones.
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod pointer to the NucleotideSubstitutionModel to use in
   *        the three positions. It is owned by the instance.
   * @param pfreq pointer to the std::shared_ptr<FrequencySet> equilibrium frequencies
   * @param pdist optional pointer to the AlphabetIndex2 amino-acids
   *        distance object.
   * @param paramSynRate is true iff synonymous rate is parametrised
   *        (default=false).
   */
  CodonDistanceFrequenciesSubstitutionModel(
      std::shared_ptr<const GeneticCode> gCode,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
      std::unique_ptr<CodonFrequencySetInterface> pfreq,
      std::shared_ptr<const AlphabetIndex2> pdist = nullptr,
      bool paramSynRate = false);

  /**
   * @brief Build a new CodonDistanceFrequenciesSubstitutionModel object
   * from three pointers to AbstractSubstitutionModels. NEW
   * AbstractSubstitutionModels are copied from the given ones.
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod1, pmod2, pmod3 are pointers to the
   *   NucleotideSubstitutionModel to use in the three positions.
   *   All the models must be different objects to avoid redundant
   *   parameters.  They are owned by the instance.
   * @param pfreq pointer to the std::shared_ptr<FrequencySet> equilibrium frequencies
   * @param pdist optional pointer to the AlphabetIndex2 amino-acids
   *   distance object.
   * @param paramSynRate is true iff synonymous rate is parametrised
   *   (default=false).
   */
  CodonDistanceFrequenciesSubstitutionModel(
      std::shared_ptr<const GeneticCode> gCode,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
      std::unique_ptr<CodonFrequencySetInterface> pfreq,
      std::shared_ptr<const AlphabetIndex2> pdist = nullptr,
      bool paramSynRate = false);

  virtual ~CodonDistanceFrequenciesSubstitutionModel() {}

  CodonDistanceFrequenciesSubstitutionModel* clone() const override
  {
    return new CodonDistanceFrequenciesSubstitutionModel(*this);
  }

public:
  void fireParameterChanged(const ParameterList& parameterlist) override;

  std::string getName() const override;

  double getCodonsMulRate(size_t i, size_t j) const override;

  void setNamespace(const std::string&) override;

  void setFreq(std::map<int, double>& frequencies) override;

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    return AbstractCodonFrequenciesSubstitutionModel::codonFrequencySet();
  }

  const FrequencySetInterface& frequencySet() const override
  {
    return AbstractCodonFrequenciesSubstitutionModel::frequencySet();
  }

  bool hasCodonFrequencySet() const override
  {
    return AbstractCodonFrequenciesSubstitutionModel::hasCodonFrequencySet();
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_CODONDISTANCEFREQUENCIESSUBSTITUTIONMODEL_H
