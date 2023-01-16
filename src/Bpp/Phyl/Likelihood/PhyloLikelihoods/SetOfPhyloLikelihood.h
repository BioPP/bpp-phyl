//
// File: SetOfAbstractPhyloLikelihood.h
// Authors:
//   Laurent Guéguen
// Created: jeudi 8 octobre 2015, ÃÂ  14h 33
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

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SETOFPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SETOFPHYLOLIKELIHOOD_H


#include "AbstractPhyloLikelihood.h"
#include "PhyloLikelihood.h"
#include "PhyloLikelihoodContainer.h"

namespace bpp
{
/**
 * @brief The SetOfPhyloLikelihood interface, to manage a
 * subset of PhyloLikelihoods from a given
 * PhyloLikelihoodContainer
 */
class SetOfPhyloLikelihoodInterface :
  public virtual PhyloLikelihoodInterface
{
public:
  virtual SetOfPhyloLikelihoodInterface* clone() const = 0;

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
 * @brief The SetOfPhyloLikelihood class, to manage a
 * subset of PhyloLikelihoods from a given
 * PhyloLikelihoodContainer
 */
class AbstractSetOfPhyloLikelihood :
  public virtual SetOfPhyloLikelihoodInterface,
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
  mutable std::vector<std::shared_ptr<LikelihoodCalculation> > vLikCal_;

public:
  /**
   * @param inCollection : avoid suffix addition to parameter names
   */
  AbstractSetOfPhyloLikelihood(
      Context& context,
      std::shared_ptr<PhyloLikelihoodContainer> pC,
      bool inCollection = true,
      const std::string& prefix = "");

  AbstractSetOfPhyloLikelihood(
      Context& context,
      std::shared_ptr<PhyloLikelihoodContainer> pC,
      const std::vector<size_t>& nPhylo,
      bool inCollection = true,
      const std::string& prefix = "");

  virtual ~AbstractSetOfPhyloLikelihood() {}

protected:

  AbstractSetOfPhyloLikelihood(const AbstractSetOfPhyloLikelihood& sd) :
    AbstractPhyloLikelihood(sd),
    AbstractParametrizable(sd),
    pPhyloCont_(sd.pPhyloCont_),
    nPhylo_(sd.nPhylo_),
    vLikCal_(sd.vLikCal_)
  {}

  AbstractSetOfPhyloLikelihood& operator=(const AbstractSetOfPhyloLikelihood& sd)
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
