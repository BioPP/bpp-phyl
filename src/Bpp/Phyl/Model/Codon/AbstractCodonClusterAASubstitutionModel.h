// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_ABSTRACTCODONCLUSTERAASUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_ABSTRACTCODONCLUSTERAASUBSTITUTIONMODEL_H

#include <Bpp/Numeric/AbstractParameterAliasable.h>

#include "CodonSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Abstract class for modelling of non-synonymous and
 *  synonymous substitution rates in codon models, with AA clustered.
 *
 * @author Laurent GuÃÂ©guen
 *
 * Non-synonymous rates between amino-acids in the same cluster are
 *  multiplied with @f$\omega_C@f$, and non-synonymous rates between
 *  amino-acids in different clusters are multiplied with @f$\omega_R@f$.
 *
 *  Parameters names are \c "omegaR" and "omegaC".
 *
 *  Clusters can be defined as a vector of assignations, amino-acids
 *  ordered as 3-letters abbreviates ("Ala", "Arg", ...). Default
 *  cluster splits according to polarity and size: "AGPV", "RQEHKWY",
 *  "NDCST", "ILMF", which produces vector:
 *
 * (1,2,3,3,3,2,2,1,2,4,4,2,4,4,1,3,3,2,2,1)
 *
 * References:
 *
 * Sainudiin et al., 2005, Detecting Site-Specific Physicochemical
 *    Selective Pressures: Applications to the Class I HLA of the
 *    Human Major Histocompatibility Complex and the SRK of the Plant
 *    Sporophytic Self-Incompatibility System, Journal of Molecular
 *    Evolution 60(3):315-26
 *
 *
 * Claudia C Weber, Simon Whelan, 2019, Physicochemical Amino Acid
 *    Properties Better Describe Substitution Rates in Large
 *    Populations, Molecular Biology and Evolution, Volume 36, Issue
 *    4, April 2019, Pages 679Ã¢ÂÂ690,
 *    https://doi.org/10.1093/molbev/msz003
 */
class AbstractCodonClusterAASubstitutionModel :
  public virtual CoreCodonSubstitutionModelInterface,
  public virtual AbstractParameterAliasable
{
private:
  std::shared_ptr<const GeneticCode> pgencode_;

  double omegaR_, omegaC_;

  std::vector<unsigned int> assign_;

  std::shared_ptr<const StateMapInterface> stateMap_;

public:
  /**
   * @brief Build a new AbstractCodonClusterAASubstitutionModel object.
   *
   * @param pgencode the genetic code
   * @param prefix the Namespace
   * @param assign an paramSynRate is true iff synonymous rate is parametrised
   *       (default categories:   "AGPV", "RQEHKWY", "NDCST", "ILMF")
   */
  AbstractCodonClusterAASubstitutionModel(
      std::shared_ptr<const GeneticCode> pgencode,
      const std::string& prefix,
      const std::vector<unsigned int>& assign = {1, 2, 3, 3, 3, 2, 2, 1, 2, 4, 4, 2, 4, 4, 1, 3, 3, 2, 2, 1});

  AbstractCodonClusterAASubstitutionModel(const AbstractCodonClusterAASubstitutionModel& model) :
    AbstractParameterAliasable(model),
    pgencode_(model.pgencode_),
    omegaR_(model.omegaR_),
    omegaC_(model.omegaC_),
    assign_(model.assign_),
    stateMap_(model.stateMap_)
  {}

  AbstractCodonClusterAASubstitutionModel& operator=(
      const AbstractCodonClusterAASubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    pgencode_ = model.pgencode_;
    omegaR_ = model.omegaR_;
    omegaC_ = model.omegaC_;
    assign_ = model.assign_;
    stateMap_ = model.stateMap_;

    return *this;
  }

  AbstractCodonClusterAASubstitutionModel* clone() const override
  {
    return new AbstractCodonClusterAASubstitutionModel(*this);
  }

  virtual ~AbstractCodonClusterAASubstitutionModel() {}

public:
  void fireParameterChanged(const ParameterList& parameters) override;

  double getCodonsMulRate(size_t i, size_t j) const override;

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    throw NullPointerException("AbstractCodonClusterAASubstitutionModel::frequencySet. No associated FrequencySet.");
  }

  bool hasCodonFrequencySet() const override
  {
    return false;
  }

  const std::vector<unsigned int>& getAssign() const
  {
    return assign_;
  }

  void setFreq(std::map<int, double>& frequencies) override {}
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTCODONCLUSTERAASUBSTITUTIONMODEL_H
