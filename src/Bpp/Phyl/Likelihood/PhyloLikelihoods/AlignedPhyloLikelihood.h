// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_ALIGNEDPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_ALIGNEDPHYLOLIKELIHOOD_H


// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>

// from bpp-core

#include <Bpp/Numeric/AbstractParametrizable.h>

#include "PhyloLikelihood.h"
#include "AbstractPhyloLikelihood.h"

namespace bpp
{
/**
 * @brief The AlignedPhyloLikelihood interface, for phylogenetic likelihood.
 *
 * This interface defines the common methods needed to compute a
 * likelihood from an alignment of data (involving several process
 * & alphabets, if needed), where all sites are independent (eq it
 * does not match for HMM phylolikelihoods).
 *
 */
class AlignedPhyloLikelihoodInterface :
  public virtual PhyloLikelihoodInterface
{
public:
  AlignedPhyloLikelihoodInterface() {}
  virtual ~AlignedPhyloLikelihoodInterface() {}

  virtual AlignedPhyloLikelihoodInterface* clone() const = 0;

public:
  /**
   *
   * @name The data functions
   *
   * @{
   */

  /**
   * @brief Get the number of sites in the dataset.
   *
   * @return the number of sites in the dataset.
   */
  virtual size_t getNumberOfSites() const = 0;

  /**
   * @}
   */

  /**
   * @name The likelihood functions.
   *
   * @{
   */

  /**
   * @return The LikelihoodCalculation.
   */
  virtual AlignedLikelihoodCalculation& alignedLikelihoodCalculation() const = 0;

  /**
   * @return A shared pointer toward the LikelihoodCalculation.
   */
  virtual std::shared_ptr<AlignedLikelihoodCalculation> getAlignedLikelihoodCalculation() const = 0;

  /**
   * @brief Get the likelihood for a site.
   *
   * @param site The site index to analyse.
   * @return The likelihood for site <i>site</i>.
   */

  virtual DataLik getLikelihoodForASite(size_t site) const = 0;

  /**
   * @brief Get the log likelihood for a site, and its derivatives.
   *
   * @param site The site index to analyse.
   * @return The (D)log likelihood for site <i>site</i>.
   */

  virtual double getLogLikelihoodForASite(size_t site) const = 0;

  /**
   * @brief Get the likelihood for each site
   *
   * @return A vector with all site likelihoods.
   */

  virtual VDataLik getLikelihoodPerSite() const = 0;

/** @} */

  friend class AlignedPhyloLikelihoodSet;
};


class AbstractAlignedPhyloLikelihood :
  public virtual AlignedPhyloLikelihoodInterface,
  public virtual AbstractPhyloLikelihood
{
protected:
  size_t nbSites_;

public:
  AbstractAlignedPhyloLikelihood(Context& context, size_t nbSites) :
    AbstractPhyloLikelihood(context),
    nbSites_(nbSites)
  {}

protected:
  AbstractAlignedPhyloLikelihood(const AbstractAlignedPhyloLikelihood& aasd) :
    AbstractPhyloLikelihood(aasd),
    nbSites_(aasd.nbSites_)
  {}

  AbstractAlignedPhyloLikelihood& operator=(const AbstractAlignedPhyloLikelihood& aasd)
  {
    AbstractPhyloLikelihood::operator=(aasd);
    nbSites_ = aasd.nbSites_;
    return *this;
  }

public:
  virtual ~AbstractAlignedPhyloLikelihood() {}

  size_t getNumberOfSites() const { return nbSites_; }

  /**
   * @brief Get the likelihood for a site (on uncompressed data)
   *
   * @param site The site index to analyse.
   * @return The likelihood for site <i>site</i>.
   */
  DataLik getLikelihoodForASite(size_t site) const
  {
    return alignedLikelihoodCalculation().getLikelihoodForASite(site);
  }

  /**
   * @brief Get the log likelihood for a site, and its derivatives.
   *
   * @param site The site index to analyse.
   * @return The (D)log likelihood for site <i>site</i>.
   */
  double getLogLikelihoodForASite(size_t site) const
  {
    using namespace numeric;

    return convert(getAlignedLikelihoodCalculation()->getLogLikelihoodForASite(site));
  }

  /**
   * @brief Get the likelihood for each site.
   *
   *@return A vector with all site likelihoods.
   *
   */
  VDataLik getLikelihoodPerSite() const
  {
    return getAlignedLikelihoodCalculation()->getLikelihoodPerSite();
  }

protected:
  void setNumberOfSites(size_t nbSites)
  {
    nbSites_ = nbSites;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_ALIGNEDPHYLOLIKELIHOOD_H
