// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SINGLEDATAPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SINGLEDATAPHYLOLIKELIHOOD_H


// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignmentData.h>

// from bpp-core

#include <Bpp/Numeric/AbstractParametrizable.h>

#include "AlignedPhyloLikelihood.h"

namespace bpp
{
/**
 * @brief The SingleDataPhyloLikelihood interface, for phylogenetic likelihood.
 *
 * This interface defines the common methods needed to compute a
 * likelihood from aligned sequences.
 *
 */
class SingleDataPhyloLikelihoodInterface :
  public virtual AlignedPhyloLikelihoodInterface
{
public:
  SingleDataPhyloLikelihoodInterface() {}
  virtual ~SingleDataPhyloLikelihoodInterface() {}

  virtual SingleDataPhyloLikelihoodInterface* clone() const = 0;

public:
  /**
   * @name The data functions
   *
   * @{
   */

  /**
   * @brief Set the dataset for which the likelihood must be evaluated.
   *
   * @param sites The data set to use.
   * @param nData the number of the data
   */
  virtual void setData(std::shared_ptr<const AlignmentDataInterface> sites, size_t nData = 0) = 0;

  /**
   * @brief Get the dataset for which the likelihood must be evaluated.
   *
   * @return A pointer toward the site container where the sequences are stored.
   */
  virtual std::shared_ptr<const AlignmentDataInterface> getData() const = 0;

  /**
   * @brief Get the number of dataset concerned.
   */
  virtual size_t getNData() const = 0;

  /**
   * @brief Get the number the states.
   */
  virtual size_t getNumberOfStates() const = 0;

  /**
   * @brief Get the alphabet associated to the dataset.
   *
   * @return the alphabet associated to the dataset.
   */
  virtual std::shared_ptr<const Alphabet> getAlphabet() const = 0;

  /**
   * @}
   */
};


class AbstractSingleDataPhyloLikelihood :
  public virtual SingleDataPhyloLikelihoodInterface,
  public virtual AbstractAlignedPhyloLikelihood
{
protected:
  size_t nbStates_;

  /**
   * @brief Number of the concerned data.
   */
  size_t nData_;

public:
  AbstractSingleDataPhyloLikelihood(Context& context, size_t nbSites, size_t nbStates, size_t nData = 0) :
    AbstractPhyloLikelihood(context),
    AbstractAlignedPhyloLikelihood(context, nbSites),
    nbStates_(nbStates),
    nData_(nData)
  {}

  AbstractSingleDataPhyloLikelihood(const AbstractSingleDataPhyloLikelihood& asdpl) :
    AbstractPhyloLikelihood(asdpl),
    AbstractAlignedPhyloLikelihood(asdpl),
    nbStates_(asdpl.nbStates_),
    nData_(asdpl.nData_)
  {}

  AbstractSingleDataPhyloLikelihood& operator=(const AbstractSingleDataPhyloLikelihood& asdpl)
  {
    AbstractAlignedPhyloLikelihood::operator=(asdpl);
    nbStates_ = asdpl.nbStates_;
    nData_ = asdpl.nData_;
    return *this;
  }

  virtual ~AbstractSingleDataPhyloLikelihood() {}

  virtual void setData(std::shared_ptr<const AlignmentDataInterface> sites, size_t nData = 0)
  {
    setNumberOfSites(sites->getNumberOfSites());
    nbStates_ = sites->alphabet().getSize();
    nData_ = nData;
  }

  size_t getNData() const
  {
    return nData_;
  }

  void setNData(size_t nData)
  {
    nData_ = nData;
  }

  size_t getNumberOfStates() const { return nbStates_; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SINGLEDATAPHYLOLIKELIHOOD_H
