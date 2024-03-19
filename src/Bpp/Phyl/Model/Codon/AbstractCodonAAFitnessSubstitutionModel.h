// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_ABSTRACTCODONAAFITNESSSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_ABSTRACTCODONAAFITNESSSUBSTITUTIONMODEL_H


# ifndef _ABSTRACTCODON_AA_FITNESSSUBSTITUTIONMODEL_H_
# define _ABSTRACTCODON_AA_FITNESSSUBSTITUTIONMODEL_H_

# include "CodonSubstitutionModel.h"

#include "../FrequencySet/ProteinFrequencySet.h"

namespace bpp
{
/**
 * @brief Abstract class for modelling of ratios of substitution
 * rates between codons, whatever they are synonymous or not.
 *
 * @author Laurent GuÃÂ©guen
 *
 * The fitness of an amino acid is a value between 0 and 1 defining
 * the relative advantage of an amino acid, compared to others. If
 * an amino acid @f$i@f$ has a fitness @f$\phi_i@f$ and another one
 * (@f$j@f$) has a fitness @f$\phi_j@f$, the substitution rate from
 * codon @f$i@f$ to codon @f$j@f$ is multiplied by
 * \f[-\frac{\log\left((\frac{\phi_i}{\phi_j})^{s}\right)}{1-\left(\frac{\phi_i}{\phi_j}\right)^{s}}\f]
 *
 * where @f$s@f$ is a positive Ns coefficient (default value:
 * 1).
 *
 * The set of fitnesses is implemented through a Protein
 * FrequencySet object. The parameters are named \c
 * "fit_NameOfTheParameterInTheFrequencySet".
 */
class AbstractCodonAAFitnessSubstitutionModel :
  public virtual CoreCodonSubstitutionModelInterface,
  public virtual AbstractParameterAliasable
{
private:
  std::shared_ptr<FrequencySetInterface> pfitset_;

  std::shared_ptr<const GeneticCode> pgencode_;

  std::string fitName_;

  std::shared_ptr<const StateMapInterface> stateMap_;

  std::shared_ptr<const StateMapInterface> protStateMap_;

  /**
   * @brief The Ns of the model (default: 1),  The generator (and all
   * its vectorial components) is independent of the rate, since it
   * should be normalized.
   */
  double Ns_;

public:
  AbstractCodonAAFitnessSubstitutionModel(
      std::shared_ptr<FrequencySetInterface> pfitset,
      std::shared_ptr<const GeneticCode> pgencode,
      const std::string& prefix);

  AbstractCodonAAFitnessSubstitutionModel(const AbstractCodonAAFitnessSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    pfitset_(model.pfitset_->clone()),
    pgencode_(model.pgencode_),
    fitName_(model.fitName_),
    stateMap_(model.stateMap_),
    protStateMap_(pfitset_->getStateMap()),
    Ns_(1)
  {}

  AbstractCodonAAFitnessSubstitutionModel& operator=(const AbstractCodonAAFitnessSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    pfitset_.reset(model.pfitset_->clone());
    pgencode_ = model.pgencode_;
    fitName_ = model.fitName_;
    stateMap_ = model.stateMap_;
    protStateMap_ = pfitset_->getStateMap();
    Ns_ = 1;

    return *this;
  }

  AbstractCodonAAFitnessSubstitutionModel* clone() const override
  {
    return new AbstractCodonAAFitnessSubstitutionModel(*this);
  }

  virtual ~AbstractCodonAAFitnessSubstitutionModel();

public:
  void fireParameterChanged(const ParameterList& parameters) override;

  void setFreq(std::map<int, double>& frequencies) override;

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    throw NullPointerException("AbstractCodonAAFitnessSubstitutionModel::codonFrequencySet. This model does not take codon frequencies. See aaFitness.");
  }

  bool hasCodonFrequencySet() const override
  {
    return false;
  }

  void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix);
    pfitset_->setNamespace(prefix + fitName_);
  }

  double getCodonsMulRate(size_t i, size_t j) const override;

  const FrequencySetInterface& aaFitness() const { return *pfitset_; }

  std::shared_ptr<const FrequencySetInterface> getAAFitness() const { return pfitset_; }

  void addNsParameter()
  {
    addParameter_(new Parameter(getNamespace() + "Ns", 1, std::make_shared<IntervalConstraint>(NumConstants::MILLI(), 100, true, true)));
  }
};
} // end of namespace bpp

# endif // _ABSTRACTCODON_AA_FITNESSSUBSTITUTIONMODEL_H_
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTCODONAAFITNESSSUBSTITUTIONMODEL_H
