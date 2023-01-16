//
// File: AbstractCodonPhaseFrequenciesSubstitutionModel.h
// Authors:
//   vendredi 23 septembre 2011, Ã 16h 29
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

#ifndef BPP_PHYL_MODEL_CODON_ABSTRACTCODONPHASEFREQUENCIESSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_ABSTRACTCODONPHASEFREQUENCIESSUBSTITUTIONMODEL_H


#include "../FrequencySet/CodonFrequencySet.h"
#include "CodonSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Abstract Class for substitution models on codons
 *  parametrized by a frequency.
 *
 * @author Laurent GuÃÂ©guen
 *
 * This class should be used with models which equilibrium
 * distribution is fixed, ans does not depend on the parameters.
 * Otherwise there may be problems of identifiability of the
 * parameters.
 *
 * If we denote @f$F@f$ the given frequencies for codons,
 * @f$F_{j_k}@f$ is the frequency of letter @f$j@f$ in phase
 * @f$k@f$.
 *
 * For codons @f$i=i_1i_2i_3@f$ and @f$j=j_1j_2j_3@f$, the generator
 * term defined from inherited and inheriting classes,
 * @f$Q_{ij})@f$, is multiplied by the product of the @f$F_{j_k}@f$
 * for each @f$k \in 1, 2, 3@f$ such that @f$i_k \neq j_k@f$.
 */
class AbstractCodonPhaseFrequenciesSubstitutionModel :
  public virtual CoreCodonSubstitutionModelInterface,
  public virtual AbstractParameterAliasable
{
private:
  /**
   * @brief Position dependent version of Codon Frequencies Set
   */
  std::unique_ptr<CodonFrequencySetInterface> posFreqSet_;
  std::string freqName_;

public:
  /**
   * @brief Build a AbstractCodonPhaseFrequenciesSubstitutionModel instance
   *
   * @param pfreq pointer to the AbstractFrequencySet equilibrium frequencies.
   *        It is owned by the instance.
   * @param prefix the Namespace
   */
  AbstractCodonPhaseFrequenciesSubstitutionModel(
    std::unique_ptr<CodonFrequencySetInterface> pfreq,
    const std::string& prefix);

  AbstractCodonPhaseFrequenciesSubstitutionModel(const AbstractCodonPhaseFrequenciesSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    posFreqSet_(model.posFreqSet_->clone()),
    freqName_(model.freqName_)
  {}

  AbstractCodonPhaseFrequenciesSubstitutionModel& operator=(const AbstractCodonPhaseFrequenciesSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    posFreqSet_.reset(model.posFreqSet_->clone());
    freqName_ = model.freqName_;

    return *this;
  }

  AbstractCodonPhaseFrequenciesSubstitutionModel* clone() const override
  {
    return new AbstractCodonPhaseFrequenciesSubstitutionModel(*this);
  }

  virtual ~AbstractCodonPhaseFrequenciesSubstitutionModel();

  void fireParameterChanged(const ParameterList& parameters) override;

  void setFreq(std::map<int, double>& frequencies) override;

  void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix);
    posFreqSet_->setNamespace(prefix + freqName_);
  }

  double getCodonsMulRate(size_t, size_t) const override;

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    return *posFreqSet_;
  }

  bool hasCodonFrequencySet() const override
  {
    return (posFreqSet_ != nullptr);
  }

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTCODONPHASEFREQUENCIESSUBSTITUTIONMODEL_H
