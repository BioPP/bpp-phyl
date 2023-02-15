//
// File: AbstractCodonAAFitnessSubstitutionModel.h
// Authors:
//   Laurent Guéguen
// Created: mercredi 8 novembre 2017, ÃÂ  21h 11
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

  const CodonFrequencySetInterface& codonFrequencySet() const override {
    throw NullPointerException("AbstractCodonAAFitnessSubstitutionModel::codonFrequencySet. This model does not take codon frequencies. See aaFitness.");
  }
  
  bool hasCodonFrequencySet() const override {
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

# endif//_ABSTRACTCODON_AA_FITNESSSUBSTITUTIONMODEL_H_
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTCODONAAFITNESSSUBSTITUTIONMODEL_H
