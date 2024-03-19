// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_ABSTRACTDISCRETERATESACROSSSITESTREELIKELIHOOD_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_ABSTRACTDISCRETERATESACROSSSITESTREELIKELIHOOD_H


#include "../../Model/SubstitutionModel.h"
#include "AbstractTreeLikelihood.h"
#include "DiscreteRatesAcrossSitesTreeLikelihood.h"

namespace bpp
{
/**
 * @brief Partial implementation of the DiscreteRatesAcrossSitesTreeLikelihood interface.
 *
 * It contains a pointer toward a DiscreteDistribution object.
 * This object may be shared by several instances of the class.
 */
class AbstractDiscreteRatesAcrossSitesTreeLikelihood :
  public AbstractTreeLikelihood,
  public virtual DiscreteRatesAcrossSitesTreeLikelihoodInterface
{
protected:
  std::shared_ptr<DiscreteDistributionInterface> rateDistribution_;

public:
  AbstractDiscreteRatesAcrossSitesTreeLikelihood(
      std::shared_ptr<DiscreteDistributionInterface> rDist,
      bool verbose = true
      );

  AbstractDiscreteRatesAcrossSitesTreeLikelihood(
      const AbstractDiscreteRatesAcrossSitesTreeLikelihood& tl) :
    AbstractTreeLikelihood(tl),
    rateDistribution_(tl.rateDistribution_)
  {}

  AbstractDiscreteRatesAcrossSitesTreeLikelihood& operator=(
      const AbstractDiscreteRatesAcrossSitesTreeLikelihood& tl)
  {
    AbstractTreeLikelihood::operator=(tl);
    rateDistribution_ = tl.rateDistribution_;
    return *this;
  }

  virtual ~AbstractDiscreteRatesAcrossSitesTreeLikelihood() {}

public:
  /**
   * @name The TreeLikelihood interface.
   *
   * Other methods are implemented in the AbstractTreeLikelihood class.
   *
   * @{
   */
  double getLikelihoodForASiteForAState (size_t site, int state) const;
  double getLogLikelihoodForASiteForAState(size_t site, int state) const;
  ParameterList getDerivableParameters() const;
  ParameterList getNonDerivableParameters() const;
  VVdouble getTransitionProbabilities(int nodeId, size_t siteIndex) const;
  /** @} */

  /**
   * @name The DiscreteRatesAcrossSites interface implementation:
   *
   * @{
   */
  std::shared_ptr<const DiscreteDistributionInterface> getRateDistribution() const
  {
    return rateDistribution_;
  }

  std::shared_ptr<DiscreteDistributionInterface> getRateDistribution()
  {
    return rateDistribution_;
  }

  const DiscreteDistributionInterface& rateDistribution() const
  {
    return *rateDistribution_;
  }

  DiscreteDistributionInterface& rateDistribution()
  {
    return *rateDistribution_;
  }

  size_t getNumberOfClasses() const
  {
    return rateDistribution_->getNumberOfCategories();
  }

  ParameterList getRateDistributionParameters() const;
  VVdouble getLikelihoodPerSitePerRateClass() const;
  VVdouble getLogLikelihoodPerSitePerRateClass() const;
  VVVdouble getLikelihoodPerSitePerRateClassPerState() const;
  VVVdouble getLogLikelihoodPerSitePerRateClassPerState() const;
  VVdouble getPosteriorProbabilitiesPerRate() const;
  Vdouble getRateWithMaxPostProbPerSite() const;
  std::vector<size_t> getRateClassWithMaxPostProbPerSite() const;
  Vdouble getPosteriorRatePerSite() const;
  /** @} */

  /**
   * @name Generic tools to deal with likelihood arrays
   *
   * @{
   */

  /**
   * @brief Set all conditional likelihoods to 1.
   *
   * @param likelihoodArray the likelihood array.
   */
  static void resetLikelihoodArray(VVVdouble& likelihoodArray);

  /**
   * @brief Print the likelihood array to terminal (debugging tool).
   *
   * @param likelihoodArray the likelihood array.
   */
  static void displayLikelihoodArray(const VVVdouble& likelihoodArray);

  /** @} */
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_ABSTRACTDISCRETERATESACROSSSITESTREELIKELIHOOD_H
