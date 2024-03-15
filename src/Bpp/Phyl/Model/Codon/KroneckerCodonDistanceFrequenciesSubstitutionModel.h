// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_KRONECKERCODONDISTANCEFREQUENCIESSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_KRONECKERCODONDISTANCEFREQUENCIESSUBSTITUTIONMODEL_H


#include "AbstractCodonDistanceSubstitutionModel.h"
#include "AbstractCodonFrequenciesSubstitutionModel.h"
#include "AbstractKroneckerCodonSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Class for non-synonymous substitution models on codons with
 * parameterized equilibrium frequencies and nucleotidic models,
 * with allowed multiple substitutions.
 *
 * Rates of multiple substitutions equal the product of single
 * substitutions involved, before mgmt of selection and removing stop
 * codons.
 *
 * @author Laurent GuÃÂ©guen
 *
 * This class should be used with models which equilibrium
 * distribution is fixed, ans does not depend on the parameters.
 * Otherwise there may be problems of identifiability of the
 * parameters.
 *
 * See description in AbstractKroneckerCodonDistanceSubstitutionModel
 * and AbstractCodonFrequenciesSubstitutionModel class.
 *
 * The additional parameters to
 * AbstractCodonFrequenciesSubstitutionModel are the rates of
 * nonsynonymous over synonymous substitutions.
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
 */

class KroneckerCodonDistanceFrequenciesSubstitutionModel :
  public AbstractKroneckerCodonSubstitutionModel,
  public AbstractCodonDistanceSubstitutionModel,
  public AbstractCodonFrequenciesSubstitutionModel
{
public:
  /**
   * @brief Build a new KroneckerCodonDistanceFrequenciesSubstitutionModel object
   * from three pointers to AbstractSubstitutionModels. NEW
   * AbstractSubstitutionModels are copied from the given ones.
   *
   * Any number of any positions can change simultaneously
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod pointer to the NucleotideSubstitutionModel to use in
   *        the three positions. It is owned by the instance.
   * @param pfreq pointer to the std::shared_ptr<FrequencySet> equilibrium frequencies
   * @param pdist optional pointer to the AlphabetIndex2 amino-acids
   *        distance object.
   */
  KroneckerCodonDistanceFrequenciesSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    std::unique_ptr<CodonFrequencySetInterface> pfreq,
    std::shared_ptr<const AlphabetIndex2> pdist = nullptr);

  /**
   * @brief Build a new KroneckerCodonDistanceFrequenciesSubstitutionModel object
   * from three pointers to AbstractSubstitutionModels. NEW
   * AbstractSubstitutionModels are copied from the given ones.
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod pointer to the NucleotideSubstitutionModel to use in
   *        the three positions. It is owned by the instance.
   * @param vPos a vector of sets of simultaneously changing
   *   positions.
   * @param pfreq pointer to the std::shared_ptr<FrequencySet> equilibrium frequencies
   * @param pdist optional pointer to the AlphabetIndex2 amino-acids
   *        distance object.
   */
  KroneckerCodonDistanceFrequenciesSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    const std::vector<std::set< size_t> >& vPos,
    std::unique_ptr<CodonFrequencySetInterface> pfreq,
    std::shared_ptr<const AlphabetIndex2> pdist = nullptr);

  /**
   * @brief Build a new KroneckerCodonDistanceFrequenciesSubstitutionModel object
   * from three pointers to AbstractSubstitutionModels. NEW
   * AbstractSubstitutionModels are copied from the given ones.
   *
   * Any number of any positions can change simultaneously
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod1, pmod2, pmod3 are pointers to the
   *   NucleotideSubstitutionModel to use in the three positions.
   *   All the models must be different objects to avoid redundant
   *   parameters.  They are owned by the instance.
   * @param pfreq pointer to the std::shared_ptr<FrequencySet> equilibrium frequencies
   * @param pdist optional pointer to the AlphabetIndex2 amino-acids
   *   distance object.
   */
  KroneckerCodonDistanceFrequenciesSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    std::unique_ptr<CodonFrequencySetInterface> pfreq,
    std::shared_ptr<const AlphabetIndex2> pdist = nullptr);

  /**
   * @brief Build a new KroneckerCodonDistanceFrequenciesSubstitutionModel object
   * from three pointers to AbstractSubstitutionModels. NEW
   * AbstractSubstitutionModels are copied from the given ones.
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod1, pmod2, pmod3 are pointers to the
   *   NucleotideSubstitutionModel to use in the three positions.
   *   All the models must be different objects to avoid redundant
   *   parameters.  They are owned by the instance.
   * @param vPos a vector of sets of simultaneously changing
   *   positions.
   * @param pfreq pointer to the std::shared_ptr<FrequencySet> equilibrium frequencies
   * @param pdist optional pointer to the AlphabetIndex2 amino-acids
   *   distance object.
   */
  KroneckerCodonDistanceFrequenciesSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    const std::vector<std::set< size_t> >& vPos,
    std::unique_ptr<CodonFrequencySetInterface> pfreq,
    std::shared_ptr<const AlphabetIndex2> pdist = nullptr);

  virtual ~KroneckerCodonDistanceFrequenciesSubstitutionModel() {}

  KroneckerCodonDistanceFrequenciesSubstitutionModel* clone() const override
  {
    return new KroneckerCodonDistanceFrequenciesSubstitutionModel(*this);
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

  bool hasCodonFrequencySet() const override
  {
    return AbstractCodonFrequenciesSubstitutionModel::hasCodonFrequencySet();
  }

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_KRONECKERCODONDISTANCEFREQUENCIESSUBSTITUTIONMODEL_H
