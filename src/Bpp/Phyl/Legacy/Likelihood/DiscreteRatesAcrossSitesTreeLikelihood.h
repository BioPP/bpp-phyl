// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_DISCRETERATESACROSSSITESTREELIKELIHOOD_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_DISCRETERATESACROSSSITESTREELIKELIHOOD_H

#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

#include "TreeLikelihood.h"

namespace bpp
{
/**
 * @brief Interface for rate across sites (RAS) implementation.
 *
 * This interface provides methods for dealing with RAS models.
 */
class DiscreteRatesAcrossSitesTreeLikelihoodInterface :
  public virtual TreeLikelihoodInterface
{
public:
  DiscreteRatesAcrossSitesTreeLikelihoodInterface() {}
  virtual ~DiscreteRatesAcrossSitesTreeLikelihoodInterface() {}

public:
  /**
   * @brief Get the rate distribution used for the computation.
   *
   * @return A const pointer toward the rate distribution of this instance.
   */
  virtual std::shared_ptr<const DiscreteDistributionInterface> getRateDistribution() const = 0;

  /**
   * @brief Get the rate distribution used for the computation.
   *
   * @return A pointer toward the rate distribution of this instance.
   */
  virtual std::shared_ptr<DiscreteDistributionInterface> getRateDistribution() = 0;

  /**
   * @brief Get the rate distribution used for the computation.
   *
   * @return A const reference toward the rate distribution of this instance.
   */
  virtual const DiscreteDistributionInterface& rateDistribution() const = 0;

  /**
   * @brief Get the rate distribution used for the computation.
   *
   * @return A reference toward the rate distribution of this instance.
   */
  virtual DiscreteDistributionInterface& rateDistribution() = 0;

  /**
   * @brief Get the likelihood for a site knowing its rate class.
   *
   * @param site      The site index.
   * @param rateClass The rate class index.
   * @return The likelihood for the specified site and rate class.
   */
  virtual double getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const = 0;

  /**
   * @brief Get the logarithm of the likelihood for a site knowing its rate class.
   *
   * @param site      The site index.
   * @param rateClass The rate class index.
   * @return The logarithm of the likelihood for the specified site and rate class.
   */
  virtual double getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const = 0;

  /**
   * @brief Get the likelihood for a site knowing its rate class and its ancestral state.
   *
   * @param site      The site index.
   * @param rateClass The rate class index.
   * @param state     The ancestral state.
   * @return The likelihood for the specified site and rate class and ancestral state.
   */
  virtual double getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const = 0;

  /**
   * @brief Get the logarithm of the likelihood for a site knowing its rate class and its ancestral state.
   *
   * @param site      The site index.
   * @param rateClass The rate class index.
   * @param state     The ancestral state.
   * @return The logarithm of the likelihood for the specified site and rate class and ancestral state..
   */
  virtual double getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const = 0;

  /**
   * @brief Get the likelihood for each site and each rate class.
   *
   * @return A two-dimension vector with all likelihoods.
   */
  virtual VVdouble getLikelihoodPerSitePerRateClass() const = 0;

  /**
   * @brief Get the logarithm of the likelihood for each site and each rate class.
   *
   * @return A two-dimension vector with all log likelihoods:
   * <code>V[i][j] =</code> likelihood of site i and rate class j.
   */
  virtual VVdouble getLogLikelihoodPerSitePerRateClass() const = 0;

  /**
   * @brief Get the likelihood for each site and each rate class and each state.
   *
   * @return A three-dimension vector with all likelihoods.
   */
  virtual VVVdouble getLikelihoodPerSitePerRateClassPerState() const = 0;

  /**
   * @brief Get the logarithm of the likelihood for each site and each rate class and each state.
   *
   * @return A three-dimension vector with all log likelihoods:
   * <code>V[i][j][k} =</code> likelihood of site i and rate class j and state k.
   */
  virtual VVVdouble getLogLikelihoodPerSitePerRateClassPerState() const = 0;

  /**
   * @brief Get the posterior probability for each site of belonging to a
   * particular rate class.
   *
   * @return A two-dimension vector with all posterior probabilities:
   * <code>V[i][j] =</code> probablity for site i of belonging to rate class j.
   */
  virtual VVdouble getPosteriorProbabilitiesPerRate() const = 0;

  /**
   * @brief Get the posterior rate class (the one with maximum posterior
   * probability) for each site.
   *
   * @return A vector with all rate classes indexes.
   */
  virtual std::vector<size_t> getRateClassWithMaxPostProbPerSite() const = 0;

  /**
   * @brief Get the posterior rate (the one with maximum posterior
   * probability) for each site.
   *
   * @return A vector with all rate classes indexes.
   */
  virtual Vdouble getRateWithMaxPostProbPerSite() const = 0;

  /**
   * @brief Get the posterior rate, i.e. averaged over all classes
   * and weighted with posterior probabilities, for each site.
   *
   * @return A vector with all rates.
   */
  virtual Vdouble getPosteriorRatePerSite() const = 0;

  /**
   * @brief Get the parameters associated to the rate distirbution.
   *
   * @return A ParameterList object with all rate distribution parameters.
   */
  virtual ParameterList getRateDistributionParameters() const = 0;

  /**
   * @brief Get the number of classes.
   *
   * @return The number of classes.
   */
  virtual size_t getNumberOfClasses() const = 0;

  /**
   * @brief Retrieves all Pij(t) for a particular branch, defined by the upper node.
   *
   * These intermediate results may be used by other methods.
   *
   * @param nodeId The node defining the branch of interest.
   * @param siteIndex The position in the alignment.
   * @return An array of dimension 3, where a[c][x][y] is the probability of substituting from x to y while being in rate class c.
   */
  virtual VVVdouble getTransitionProbabilitiesPerRateClass(int nodeId, size_t siteIndex) const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_DISCRETERATESACROSSSITESTREELIKELIHOOD_H
