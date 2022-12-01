//
// File: SENCA.h
// Authors:
//   Fanny Pouyet
// Created: mars 2012
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

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
    std::shared_ptr<NucleotideSubstitutionModelInterface> pmod,
    std::shared_ptr<FrequencySetInterface> pfit,
    const AlphabetIndex2* pdist = nullptr);

  SENCA(
    std::shared_ptr<const GeneticCode> gCode,
    std::shared_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::shared_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::shared_ptr<NucleotideSubstitutionModelInterface> pmod3,
    std::shared_ptr<FrequencySetInterface> pfit,
    std::shared_ptr<const AlphabetIndex2> pdist = nullptr);

  virtual ~SENCA() {}

  SENCA* clone() const override
  {
    return new SENCA(*this);
  }

public:
  void fireParameterChanged(const ParameterList& parameterlist) override;

  std::string getName() const override;

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

  const FrequencySetInterface& frequencySet() const override
  {
    return AbstractCodonFitnessSubstitutionModel::frequencySet();
  }

  std::shared_ptr<const FrequencySetInterface> getFrequencySet() const override
  {
    return AbstractCodonFitnessSubstitutionModel::getFrequencySet();
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_SENCA_H
