// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_SENCA_H
#define BPP_PHYL_MODEL_CODON_SENCA_H


#include "AbstractCodonDistanceSubstitutionModel.h"
#include "AbstractCodonFitnessSubstitutionModel.h"
#include "AbstractCodonSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Class for non-synonymous and synonymous substitution models
 * on codons with parameterized equilibrium frequencies and
 * nucleotidic basic models.
 *
 * @author Fanny Pouyet, Laurent Guéguen
 *
 * See description in AbstractCodonDistanceSubstitutionModel
 * class, AbstractCodonFitnessSubstitutionModel class.
 *
 * Only substitutions with one letter changed are accepted. </p>
 *
 * The ratio between non-synonymous and synonymous substitutions
 * rates is @f$\beta@f$ with positive parameter \c "beta".
 *
 * The fitness of a codon is a value between 0 and 1 defining the
 * relative advantage of a codon, compared to others. If a codon
 * @f$i@f$ has a fitness @f$\phi_i@f$ and another one (@f$j@f$)
 * has a fitness @f$\phi_j@f$, the substitution rate from codon
 * @f$i@f$ to codon @f$j@f$ is multiplied by
 * \f[-\frac{ \log(\frac{\phi_i}{\phi_j})}{1-\frac{\phi_i}{\phi_j}}\f]
 *
 * The set of fitnesses is implemented through a Codon
 * FrequencySet object. The parameters are named \c
 * "fit_NameOfTheParameterInTheFrequencySet".
 *
 * Note: equilibrium frequencies are computed from the generator. To
 * be done : implement analytic computation.
 *
 * Reference:
 * - Pouyet & al, Genome Biology and Evolution, 2016
 */
class SENCA :
  public virtual SubstitutionModelInterface,
  public AbstractCodonSubstitutionModel,
  public AbstractCodonDistanceSubstitutionModel,
  public AbstractCodonFitnessSubstitutionModel
{
public:
  SENCA(
      std::shared_ptr<const GeneticCode> gCode,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
      std::unique_ptr<FrequencySetInterface> pfit,
      std::shared_ptr<const AlphabetIndex2> pdist = nullptr);

  SENCA(
      std::shared_ptr<const GeneticCode> gCode,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
      std::unique_ptr<FrequencySetInterface> pfit,
      std::shared_ptr<const AlphabetIndex2> pdist = nullptr);

  virtual ~SENCA() {}

  SENCA* clone() const override
  {
    return new SENCA(*this);
  }

public:
  void fireParameterChanged(const ParameterList& parameterlist) override;

  std::string getName() const override
  {
    return "SENCA";
  }

  double getCodonsMulRate(size_t i, size_t j) const override;

  void setNamespace(const std::string&) override;

  /**
   * @brief set the fitness of the model from
   * given frequencies, such that the equilibrium frequencies of the
   * model matches at best the given ones.
   *
   * Matching is done in two steps : first, frequencies from nucleotide substitution model are
   * matched at best, then the resulting discrepancy (in terms of
   * ratios between the given one and the one computed by the pahse
   * frequencies) is given for matching to the fitness.
   *
   * @param frequencies the frequencies to match on.
   */
  void setFreq(std::map<int, double>& frequencies) override;

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    return AbstractCodonFitnessSubstitutionModel::codonFrequencySet();
  }

  bool hasCodonFrequencySet() const override
  {
    return AbstractCodonFitnessSubstitutionModel::hasCodonFrequencySet();
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_SENCA_H
