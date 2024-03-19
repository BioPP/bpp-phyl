// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_PAIREDSITELIKELIHOODS_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_PAIREDSITELIKELIHOODS_H


// From the STL:
#include <vector>
#include <string>

// From Bio++
#include "TreeLikelihood.h"
#include <Bpp/Exceptions.h>

namespace bpp
{
/**
 * @brief A container for paired-site likelihoods (likelihoods over
 * the same sites for different models, especially topologies).
 * An instance of this class is, roughly, a list of models, each of
 * them having a name (stored in the <i>modelNames</i> attribute) and
 * a set of site likelihoods (stored in the <i>logLikelihoods</i> attribute).
 */
class PairedSiteLikelihoods
{
private:
  std::vector<std::vector<double>> logLikelihoods_;
  std::vector<std::string> modelNames_;

public:
  PairedSiteLikelihoods();

  /**
   * @brief Build a new object from a site likelihoods array.
   *
   * @param siteLogLikelihoods An nmodels*nsites array of loglikelihoods.
   * @param modelNames <i>(Optional)</i> The names of the models.
   *
   * @throw Exception If the number of sites differ between the models,
   * or if the number of names and site loglikelihood records differ.
   */
  PairedSiteLikelihoods(
      const std::vector<std::vector<double>>& siteLogLikelihoods,
      const std::vector<std::string>& modelNames = std::vector<std::string>()
      );

  ~PairedSiteLikelihoods() {}

  /**
   * @brief Append a model.
   *
   * @param siteLogLikelihoods The loglikelihoods of the sites under this model.
   * @param modelName The name of the model.
   *
   * @throw Exception If the number of sites is not the same as in the container.
   */
  void appendModel(
    const std::vector<double>& siteLogLikelihoods,
    const std::string& modelName = "");

  /**
   * @brief Append a model.
   *
   * @param treelikelihood A TreeLikelihood record.
   *
   * @throw Exception If the number of sites is not the same as in the container.
   */
  void appendModel(const bpp::TreeLikelihoodInterface& treelikelihood);

  /**
   * @brief Append models by concatenation.
   *
   * @param psl the PairedSiteLikelihoods object to append to the caller.
   *
   * @throw Exception If the number of sites in the two object is not equal.
   */
  void appendModels(const PairedSiteLikelihoods& psl);

  /**
   * @return The site-likelihoods of all models.
   */
  const std::vector<std::vector<double>>& getLikelihoods() const
  {
    return logLikelihoods_;
  }

  /**
   * @return The model names.
   */
  const std::vector<std::string>& getModelNames() const
  {
    return modelNames_;
  }

  /** @brief Get the number of models in the container. */
  size_t getNumberOfModels() const
  {
    return logLikelihoods_.size();
  }

  /**
   * @return The number of sites for each model.
   * @throw Exception If the container is empty.
   */
  std::size_t getNumberOfSites() const
  {
    try
    {
      return logLikelihoods_.at(0).size();
    }
    catch (std::out_of_range&)
    {
      throw Exception("PairedSiteLikelihoods::nsites: The container is empty, there isn't a number of sites.");
    }
  }

  /**
   * @brief Set the name of a model.
   *
   * @param pos The position of the target model.
   * @param name The new name.
   */
  void setName(std::size_t pos, std::string& name)
  {
    modelNames_.at(pos) = name;
  }

  /**
   * @brief Compute the Expected Likelihood Weights of the models.
   *
   * The weight \f$W_m\f$ of a model is :
   * \f[
   * W_m
   * = \frac{1}{B} \sum_{b \in B} \frac{L_m^{(b)}}{\sum L_k^{(b)}}
   * = \frac{1}{B} \sum_{b \in B} \frac{exp(Y^{(b)}_m - Ymax^{(b)})}{\sum exp(Y_k^{(b)} - Y_{max}^{(b)})}
   * \f]
   * where \f$Y_k^{(b)}\f$ is the loglikelihood of model k for replicate b.
   * @return A pair of vectors containing the names and weights of the models.
   *
   * @param replicates The number of pseudoreplicates over which the weights are to be averaged.
   */
  std::pair< std::vector<std::string>, std::vector<double>> computeExpectedLikelihoodWeights(int replicates = 10000) const;

  /**
   * @brief Draw a nonparametric pseudoreplicate
   *
   * @param length The length of the data.
   * @param scaling The length of the pseudoreplicate, in fraction of
   *  the length of the data.
   *
   * @return A vector of the same length as the data, containing the count
   * of each element in the pseudoreplicate.
   */
  static std::vector<int> bootstrap(std::size_t length, double scaling = 1);
};
} // namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_PAIREDSITELIKELIHOODS_H
