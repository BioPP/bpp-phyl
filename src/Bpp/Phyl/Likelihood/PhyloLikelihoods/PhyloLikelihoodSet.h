// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SETOFPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SETOFPHYLOLIKELIHOOD_H


#include "AbstractPhyloLikelihood.h"
#include "PhyloLikelihood.h"
#include "PhyloLikelihoodContainer.h"

namespace bpp
{
/**
 * @brief The PhyloLikelihoodSet interface, to manage a
 * subset of PhyloLikelihoods from a given
 * PhyloLikelihoodContainer
 */
class PhyloLikelihoodSetInterface :
  public virtual PhyloLikelihoodInterface
{
public:
  virtual PhyloLikelihoodSetInterface* clone() const = 0;

public:
  virtual std::shared_ptr<PhyloLikelihoodContainer> getPhyloContainer() = 0;

  virtual std::shared_ptr<const PhyloLikelihoodContainer> getPhyloContainer() const = 0;

  virtual const std::vector<size_t>& getNumbersOfPhyloLikelihoods() const = 0;

  /**
   *
   * @brief adds a PhyloLikelihood already stored in the
   * PhyloLikelihoodContainer, iff it is an
   * AbstractPhyloLikelihood.
   *
   * @param nPhyl  number of the phylolikelihood
   * @param suff for parameters names if use specific parameters names
   *
   * @return if the PhyloLikelihood has been added.
   */
  virtual bool addPhyloLikelihood(size_t nPhyl, const std::string& suff = "") = 0;

  /**
   *
   * @name The PhyloLikelihood storage.
   *
   * @{
   */
  virtual bool hasPhyloLikelihood(size_t nPhyl) = 0;

  virtual std::shared_ptr<const PhyloLikelihoodInterface> getPhyloLikelihood(size_t nPhyl) const = 0;

  virtual std::shared_ptr<PhyloLikelihoodInterface> getPhyloLikelihood(size_t nPhyl) = 0;
};


/**
 * @brief The PhyloLikelihoodSet class, to manage a
 * subset of PhyloLikelihoods from a given
 * PhyloLikelihoodContainer
 */
class AbstractPhyloLikelihoodSet :
  public virtual PhyloLikelihoodSetInterface,
  public virtual AbstractPhyloLikelihood,
  public virtual AbstractParametrizable
{
protected:
  /**
   * @brief pointer to a  PhyloLikelihoodContainer
   */
  std::shared_ptr<PhyloLikelihoodContainer> pPhyloCont_;

  /**
   * @brief vector of AbstractPhyloLikelihood numbers
   */
  std::vector<size_t> nPhylo_;

  /**
   * vector of pointers towards LikelihoodCalculation, used
   * for the global likelihood.
   */
  mutable std::vector<std::shared_ptr<LikelihoodCalculation>> vLikCal_;

public:
  /**
   * @param inCollection : avoid suffix addition to parameter names
   */
  AbstractPhyloLikelihoodSet(
      Context& context,
      std::shared_ptr<PhyloLikelihoodContainer> pC,
      bool inCollection = true,
      const std::string& prefix = "");

  AbstractPhyloLikelihoodSet(
      Context& context,
      std::shared_ptr<PhyloLikelihoodContainer> pC,
      const std::vector<size_t>& nPhylo,
      bool inCollection = true,
      const std::string& prefix = "");

  virtual ~AbstractPhyloLikelihoodSet() {}

protected:
  AbstractPhyloLikelihoodSet(const AbstractPhyloLikelihoodSet& sd) :
    AbstractPhyloLikelihood(sd),
    AbstractParametrizable(sd),
    pPhyloCont_(sd.pPhyloCont_),
    nPhylo_(sd.nPhylo_),
    vLikCal_(sd.vLikCal_)
  {}

  AbstractPhyloLikelihoodSet& operator=(const AbstractPhyloLikelihoodSet& sd)
  {
    AbstractPhyloLikelihood::operator=(sd);
    AbstractParametrizable::operator=(sd);
    pPhyloCont_ = sd.pPhyloCont_;
    nPhylo_ = sd.nPhylo_;
    vLikCal_ = sd.vLikCal_;
    return *this;
  }

public:
  std::shared_ptr<PhyloLikelihoodContainer> getPhyloContainer() override
  {
    return pPhyloCont_;
  }

  std::shared_ptr<const PhyloLikelihoodContainer> getPhyloContainer() const override
  {
    return pPhyloCont_;
  }

  const std::vector<size_t>& getNumbersOfPhyloLikelihoods() const override
  {
    return nPhylo_;
  }

  /**
   *
   * @brief adds a PhyloLikelihood already stored in the
   * PhyloLikelihoodContainer, iff it is an
   * AbstractPhyloLikelihood.
   *
   * @param nPhyl  number of the phylolikelihood
   * @param suff for parameters names if use specific parameters names
   *
   * @return if the PhyloLikelihood has been added.
   */

  virtual bool addPhyloLikelihood(size_t nPhyl, const std::string& suff = "") override;

  /**
   *
   * @name The PhyloLikelihood storage.
   *
   * @{
   */
  bool hasPhyloLikelihood(size_t nPhyl) override
  {
    return std::find(nPhylo_.begin(), nPhylo_.end(), nPhyl) != nPhylo_.end();
  }

  std::shared_ptr<const PhyloLikelihoodInterface> getPhyloLikelihood(size_t nPhyl) const override
  {
    return (*pPhyloCont_)[nPhyl];
  }


  std::shared_ptr<PhyloLikelihoodInterface> getPhyloLikelihood(size_t nPhyl) override
  {
    return (*pPhyloCont_)[nPhyl];
  }

  /**
   *
   * @}
   *
   */

  /**
   *
   * @name Inherited from PhyloLikelihood
   *
   * @{
   */


  /**
   * @return 'true' is the likelihood function has been initialized.
   */
  bool isInitialized() const override
  {
    for (auto nPhylo: nPhylo_)
    {
      if (!getPhyloLikelihood(nPhylo)->isInitialized())
        return false;
    }
    return true;
  }

public:
  virtual void fireParameterChanged(const ParameterList& params) override
  {
    for (auto nPhylo: nPhylo_)
    {
      getPhyloLikelihood(nPhylo)->matchParametersValues(params);

      // to ensure each phylolikelihood is recomputed, such as in
      // case of total aliasing
      // getAbstractPhyloLikelihood(nPhylo_[i])->update();
    }
  }

  /**
   * @name Retrieve some particular parameters subsets.
   *
   * @{
   */


  ParameterList getNonDerivableParameters() const override;

  ParameterList getDerivableParameters() const override;

  /**
   * @brief Get the branch lengths parameters.
   *
   * @return A ParameterList with all branch lengths.
   */
  ParameterList getBranchLengthParameters() const override;

  /**
   * @brief Get the parameters associated to substitution model(s).
   *
   * @return A ParameterList.
   */
  ParameterList getSubstitutionModelParameters() const override;

  /**
   * @brief Get the parameters associated to the rate distribution(s).
   *
   * @return A ParameterList.
   */

  ParameterList getRateDistributionParameters() const override;

  /**
   * @brief Get the parameters associated to the root frequencies(s).
   *
   * @return A ParameterList.
   */

  ParameterList getRootFrequenciesParameters() const override;

  /** @} */
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SETOFPHYLOLIKELIHOOD_H
