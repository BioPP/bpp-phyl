// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_ABSTRACTCODONFITNESSSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_ABSTRACTCODONFITNESSSUBSTITUTIONMODEL_H



# ifndef _ABSTRACTCODONFITNESSSUBSTITUTIONMODEL_H_
# define _ABSTRACTCODONFITNESSSUBSTITUTIONMODEL_H_

# include "CodonSubstitutionModel.h"
#include "../FrequencySet/CodonFrequencySet.h"
namespace bpp
{
/**
 * @brief Abstract class for modelling of ratios of substitution
 * rates between codons, whatever they are synonymous or not.
 *
 * @author Fanny Pouyet, Laurent Guéguen
 *
 * The fitness of a codon is a value between 0 and 1 defining the
 * relative advantage of a codon, compared to others. If a codon
 * @f$i@f$ has a fitness @f$\phi_i@f$ and another one (@f$j@f$) has
 * a fitness @f$\phi_j@f$, the substitution rate from codon @f$i@f$
 * to codon @f$j@f$ is multiplied by
 * \f[-\frac{\log(\frac{\phi_i}{\phi_j})}{1-\frac{\phi_i}{\phi_j}}\f]
 *
 * The set of fitnesses is implemented through a Codon
 * FrequencySet object. The parameters are named \c
 * "fit_NameOfTheParameterInTheFrequencySet".
 */
class AbstractCodonFitnessSubstitutionModel :
  public virtual CoreCodonSubstitutionModelInterface,
  public virtual AbstractParameterAliasable
{
private:
  std::unique_ptr<FrequencySetInterface> pfitset_;

  std::shared_ptr<const GeneticCode> pgencode_;

  std::string fitName_;

public:
  AbstractCodonFitnessSubstitutionModel(
    std::unique_ptr<FrequencySetInterface> pfitset,
    std::shared_ptr<const GeneticCode> pgencode,
    const std::string& prefix);

  AbstractCodonFitnessSubstitutionModel(const AbstractCodonFitnessSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    pfitset_(model.pfitset_->clone()),
    pgencode_(model.pgencode_),
    fitName_(model.fitName_)
  {}

  AbstractCodonFitnessSubstitutionModel& operator=(const AbstractCodonFitnessSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    pfitset_.reset(model.pfitset_->clone());
    pgencode_ = model.pgencode_;
    fitName_ = model.fitName_;
    return *this;
  }

  AbstractCodonFitnessSubstitutionModel* clone() const override
  {
    return new AbstractCodonFitnessSubstitutionModel(*this);
  }

  virtual ~AbstractCodonFitnessSubstitutionModel();

public:
  void fireParameterChanged (const ParameterList& parameters) override;

  void setFreq(std::map<int, double>& frequencies) override;

  //const FrequencySet& getFreq() const { return *pfitset_; }

  void setNamespace (const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix);
    pfitset_->setNamespace(prefix + fitName_);
  }

  double getCodonsMulRate(size_t i, size_t j) const override;

  const FrequencySetInterface& fitness() const { return *pfitset_; }
  
  //TODO (jdutheil 30/12/22) not allowed if fully encapsulated.
  //std::shared_ptr<const FrequencySetInterface> getFitness() const { return pfitset_; }

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    throw NullPointerException("AbstractCodonFitnessSubstitutionModel::frequencySet. No associated FrequencySet.");
  }

  bool hasCodonFrequencySet() const override
  {
    return false;
  }

};
} // end of namespace bpp
# endif
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTCODONFITNESSSUBSTITUTIONMODEL_H
