// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SETOFALIGNEDPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SETOFALIGNEDPHYLOLIKELIHOOD_H



// From bpp-seq:
#include <Bpp/Seq/Container/AlignmentData.h>

#include "AlignedPhyloLikelihood.h"
#include "PhyloLikelihoodSet.h"

namespace bpp
{
/**
 * @brief Joint interface SetOf+Aligned PhylLikelihood
 */
class AlignedPhyloLikelihoodSetInterface :
  public virtual PhyloLikelihoodSetInterface,
  public virtual AlignedPhyloLikelihoodInterface
{
public:
  AlignedPhyloLikelihoodSetInterface* clone() const override = 0;

  virtual std::shared_ptr<const AlignedPhyloLikelihoodInterface> getAlignedPhyloLikelihood(size_t nPhyl) const = 0;

  virtual std::shared_ptr<AlignedPhyloLikelihoodInterface> getAlignedPhyloLikelihood(size_t nPhyl) = 0;

  virtual const AlignedPhyloLikelihoodInterface& alignedPhyloLikelihood(size_t nPhyl) const = 0;

  virtual AlignedPhyloLikelihoodInterface& alignedPhyloLikelihood(size_t nPhyl) = 0;

  /**
   * @brief Get the likelihood for a site for an aligned
   * phyloLikelihood
   *
   * @param site The site index to analyse.
   * @param nPhyl the phyloLikelihood index.
   * @return The likelihood for site <i>site</i>.
   */
  virtual DataLik getLikelihoodForASiteForAPhyloLikelihood(size_t site, size_t nPhyl) const = 0;

  /**
   * @brief Get the log likelihood for a site for an aligned
   * phyloLikelihood
   *
   * @param site The site index to analyse.
   * @param nPhyl the phyloLikelihood index.
   * @return The log likelihood for site <i>site</i>.
   */
  virtual double getLogLikelihoodForASiteForAPhyloLikelihood(size_t site, size_t nPhyl) const = 0;

};

/**
 * @brief The AlignedPhyloLikelihoodSet abstract class.
 *
 * This class defines the common methods needed to compute a
 * likelihood from aligned phylogenies.
 */
class AbstractAlignedPhyloLikelihoodSet :
  public virtual AlignedPhyloLikelihoodSetInterface,
  public virtual AbstractPhyloLikelihoodSet,
  public virtual AbstractAlignedPhyloLikelihood
{
public:
  AbstractAlignedPhyloLikelihoodSet(
      Context& context,
      std::shared_ptr<PhyloLikelihoodContainer> pC,
      bool inCollection = true,
      const std::string& prefix = "");

  AbstractAlignedPhyloLikelihoodSet(
      Context& context,
      std::shared_ptr<PhyloLikelihoodContainer> pC,
      const std::vector<size_t>& nPhylo, 
      bool inCollection = true,
      const std::string& prefix = "");

protected:

  AbstractAlignedPhyloLikelihoodSet(const AbstractAlignedPhyloLikelihoodSet& soap) :
    AbstractPhyloLikelihood(soap),
    AbstractPhyloLikelihoodSet(soap),
    AbstractAlignedPhyloLikelihood(soap)
  {}

  AbstractAlignedPhyloLikelihoodSet& operator=(const AbstractAlignedPhyloLikelihoodSet& soap)
  {
    AbstractPhyloLikelihoodSet::operator=(soap);
    AbstractAlignedPhyloLikelihood::operator=(soap);
    return *this;
  }

public:

  virtual ~AbstractAlignedPhyloLikelihoodSet() {}

  /**
   *
   * @brief adds a PhyloLikelihood already stored in the m ap, iff
   * it is an AlignedPhyloLikelihood of the same size
   *
   * @param nPhyl  number of the phylolikelihood
   * @param suff for parameters names if use specific parameters names
   *
   * @return if the PhyloLikelihood has been added.
   */
  bool addPhyloLikelihood(size_t nPhyl, const std::string& suff) override;

  std::shared_ptr<const AlignedPhyloLikelihoodInterface> getAlignedPhyloLikelihood(size_t nPhyl) const override
  {
    return std::dynamic_pointer_cast<const AlignedPhyloLikelihoodInterface>((*pPhyloCont_)[nPhyl]);
  }

  std::shared_ptr<AlignedPhyloLikelihoodInterface> getAlignedPhyloLikelihood(size_t nPhyl) override
  {
    return std::dynamic_pointer_cast<AlignedPhyloLikelihoodInterface>((*pPhyloCont_)[nPhyl]);
  }

  const AlignedPhyloLikelihoodInterface& alignedPhyloLikelihood(size_t nPhyl) const override
  {
    return dynamic_cast<const AlignedPhyloLikelihoodInterface&>(*(*pPhyloCont_)[nPhyl]);
  }

  AlignedPhyloLikelihoodInterface& alignedPhyloLikelihood(size_t nPhyl) override
  {
    return dynamic_cast<AlignedPhyloLikelihoodInterface&>(*(*pPhyloCont_)[nPhyl]);
  }

  DataLik getLikelihoodForASiteForAPhyloLikelihood(size_t site, size_t nPhyl) const override
  {
    return alignedPhyloLikelihood(nPhyl).getLikelihoodForASite(site);
  }

  double getLogLikelihoodForASiteForAPhyloLikelihood(size_t site, size_t nPhyl) const override
  {
    return alignedPhyloLikelihood(nPhyl).getLogLikelihoodForASite(site);
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SETOFALIGNEDPHYLOLIKELIHOOD_H
