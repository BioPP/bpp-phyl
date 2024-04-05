// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_G2001_H
#define BPP_PHYL_MODEL_G2001_H

#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

#include "MarkovModulatedSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Galtier's 2001 covarion model.
 *
 * This model is a subclass of the so-called Markov-modulated substitution models,
 * with a Jukes-Cantor rate matrix, of parameter @f$\nu@f$.
 * the original version uses a discrete @f$\Gamma@f$ distribution for rates, but you can
 * use it with virtually any rate distribution.
 *
 * @see MarkovModulatedSubstitutionModel
 *
 * Galtier N., Maximum-likelihood phylogenetic analysis under a covarion-like model (2001).
 * _Molecular Biology and Evolution_, 18:866-73.
 */
class G2001 :
  public MarkovModulatedSubstitutionModel
{
private:
  std::unique_ptr<DiscreteDistributionInterface> rDist_;

  std::string nestedRatePrefix_;

public:
  /**
   * @brief Build a new G2001 substitution model.
   *
   * @param model The substitution model to use. May be of any alphabet type.
   * @param rDist The discrete distribution for rates. The class will own the DiscreteDistribution object,
   * which will be deleted together with this instance.
   * @param nu    The rate matrix parameter.
   * @param normalizeRateChanges Tell if the rate transition matrix should be normalized.
   */
  G2001(
      std::unique_ptr<ReversibleSubstitutionModelInterface> model,
      std::unique_ptr<DiscreteDistributionInterface> rDist,
      double nu = 1.,
      bool normalizeRateChanges = false) :
    MarkovModulatedSubstitutionModel(std::move(model), static_cast<unsigned int>(rDist->getNumberOfCategories()), normalizeRateChanges, "G01."),
    rDist_(std::move(rDist)),
    nestedRatePrefix_("rdist_" + rDist->getNamespace())
  {
    ratesFreq_ = std::vector<double>(nbRates_, 1. / static_cast<double>(nbRates_));
    rDist_->setNamespace(getNamespace() + nestedRatePrefix_);
    addParameters_(rDist_->getIndependentParameters());
    addParameter_(new Parameter("G01.nu", nu, Parameter::R_PLUS));
    updateRatesModel_();
    updateMatrices_();
  }

  G2001(const G2001& model) :
    MarkovModulatedSubstitutionModel(model),
    rDist_(model.rDist_->clone()),
    nestedRatePrefix_(model.nestedRatePrefix_)
  {}

  G2001& operator=(const G2001& model)
  {
    MarkovModulatedSubstitutionModel::operator=(model);
    rDist_.reset(model.rDist_->clone());
    nestedRatePrefix_ = model.nestedRatePrefix_;
    return *this;
  }

  virtual ~G2001() { }

  G2001* clone() const override { return new G2001(*this); }

public:
  std::string getName() const override { return "G01"; }

  /**
   * @brief Re-definition of the super-class method to update the rate distribution too.
   *
   * @param parameters The parameters that have been modified.
   */
  void fireParameterChanged(const ParameterList& parameters) override
  {
    rDist_->matchParametersValues(parameters);
    MarkovModulatedSubstitutionModel::fireParameterChanged(parameters);
  }

  /**
   * @return The rate distribution associated to this instance.
   */
  const DiscreteDistributionInterface& rateDistribution() const { return *rDist_; }

  void setNamespace(const std::string& prefix) override
  {
    MarkovModulatedSubstitutionModel::setNamespace(prefix);
    // We also need to update the namespace of the nested distribution:
    rDist_->setNamespace(prefix + nestedRatePrefix_);
  }


  double getRate() const override {  return 1.; }

  void setRate(double rate) override {}

  void addRateParameter() override {}

protected:
  void updateRatesModel_() override
  {
    double nu = getParameterValue("nu");
    for (size_t i = 0; i < nbRates_; ++i)
    {
      rates_(i, i) = rDist_->getCategory(i);
      for (size_t j = 0; j < nbRates_; ++j)
      {
        if (i == j)
        {
          ratesExchangeability_(i, j) = -static_cast<double>(nbRates_) * nu;
        }
        else
        {
          ratesExchangeability_(i, j) = static_cast<double>(nbRates_) * nu / static_cast<double>(nbRates_ - 1);
        }
      }
    }
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_G2001_H
