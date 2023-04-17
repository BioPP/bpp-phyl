//
// File: SetOfAlignedPhyloLikelihood.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: mercredi 7 octobre 2015, ÃÂ  13h 30
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

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SETOFALIGNEDPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SETOFALIGNEDPHYLOLIKELIHOOD_H



// From bpp-seq:
#include <Bpp/Seq/Container/AlignmentData.h>

#include "AlignedPhyloLikelihood.h"
#include "SetOfPhyloLikelihood.h"

namespace bpp
{
/**
 * @brief Joint interface SetOf+Aligned PhylLikelihood
 */
class SetOfAlignedPhyloLikelihoodInterface :
  public virtual SetOfPhyloLikelihoodInterface,
  public virtual AlignedPhyloLikelihoodInterface
{
public:
  SetOfAlignedPhyloLikelihoodInterface* clone() const override = 0;

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
 * @brief The SetOfAlignedPhyloLikelihood abstract class.
 *
 * This class defines the common methods needed to compute a
 * likelihood from aligned phylogenies.
 */
class AbstractSetOfAlignedPhyloLikelihood :
  public virtual SetOfAlignedPhyloLikelihoodInterface,
  public virtual AbstractSetOfPhyloLikelihood,
  public virtual AbstractAlignedPhyloLikelihood
{
public:
  AbstractSetOfAlignedPhyloLikelihood(
      Context& context,
      std::shared_ptr<PhyloLikelihoodContainer> pC,
      bool inCollection = true,
      const std::string& prefix = "");

  AbstractSetOfAlignedPhyloLikelihood(
      Context& context,
      std::shared_ptr<PhyloLikelihoodContainer> pC,
      const std::vector<size_t>& nPhylo, 
      bool inCollection = true,
      const std::string& prefix = "");

protected:

  AbstractSetOfAlignedPhyloLikelihood(const AbstractSetOfAlignedPhyloLikelihood& soap) :
    AbstractPhyloLikelihood(soap),
    AbstractSetOfPhyloLikelihood(soap),
    AbstractAlignedPhyloLikelihood(soap)
  {}

  AbstractSetOfAlignedPhyloLikelihood& operator=(const AbstractSetOfAlignedPhyloLikelihood& soap)
  {
    AbstractSetOfPhyloLikelihood::operator=(soap);
    AbstractAlignedPhyloLikelihood::operator=(soap);
    return *this;
  }

public:

  virtual ~AbstractSetOfAlignedPhyloLikelihood() {}

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
