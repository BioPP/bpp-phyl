// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PHYLOLIKELIHOOD_H

#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/ParameterList.h>

#include "../DataFlow/LikelihoodCalculation.h"

namespace bpp
{
/**
 * @brief The PhyloLikelihood interface, for phylogenetic likelihood.
 *
 * This interface defines the common methods needed to compute a likelihood
 * from a sequence alignement, usually involving one or more phylogenetic trees.
 */
class PhyloLikelihoodInterface :
  public virtual SecondOrderDerivable
{
public:
  PhyloLikelihoodInterface() {}
  virtual ~PhyloLikelihoodInterface() {}

  PhyloLikelihoodInterface* clone() const = 0;

public:
  /**
   * @name The data functions
   *
   * @{
   */

  /**
   * @return 'true' is the likelihood function has been initialized.
   */
  virtual bool isInitialized() const = 0;

  virtual const Context& context() const = 0;

  virtual Context& context() = 0;

  /**
   * @}
   */

  /**
   * @name The likelihood functions.
   *
   * @{
   */

  /**
   * @return the LikDF node where the Likelihood is computed.
   */
  virtual ValueRef<DataLik> getLikelihoodNode() const = 0;

  /**
   * @brief Get the logarithm of the likelihood for the whole dataset.
   *
   * @return The logarithm of the likelihood of the dataset.
   */
  double getLogLikelihood() const
  {
    return -getValue();
  }

  /**
   * @brief Get the derivates of the LogLikelihood.
   *
   */

  // virtual double getDLogLikelihood(const std::string& variable) const = 0;

  // virtual double getD2LogLikelihood(const std::string& variable) const = 0;

  /** @} */

  /**
   * @name Retrieve some particular independent parameters subsets.
   *
   * @{
   */

  virtual ParameterList getNonDerivableParameters() const = 0;

  virtual ParameterList getDerivableParameters() const = 0;

  /**
   * @brief Get the independent branch lengths parameters.
   *
   * @return A ParameterList with all branch lengths.
   */
  virtual ParameterList getBranchLengthParameters() const = 0;

  /**
   * @brief Get the independent parameters associated to substitution model(s).
   *
   * @return A ParameterList.
   */
  virtual ParameterList getSubstitutionModelParameters() const = 0;

  /**
   * @brief Get the independent parameters associated to the rate distribution(s).
   *
   * @return A ParameterList.
   */
  virtual ParameterList getRateDistributionParameters() const = 0;

  /**
   * @brief Get the independent parameters associated to the root frequencies(s).
   *
   * @return A ParameterList.
   */
  virtual ParameterList getRootFrequenciesParameters() const = 0;

  /** @} */

  /**
   * @return The LikelihoodCalculation.
   */
  virtual LikelihoodCalculation& likelihoodCalculation() const = 0;

  /**
   * @return A shared pointer toward the LikelihoodCalculation.
   */
  virtual std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation() const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_PHYLOLIKELIHOOD_H
