//
// File: MultiPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: jeudi 11 juillet 2013, à 21h 51
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef _MULTIPHYLOLIKELIHOOD_H_
#define _MULTIPHYLOLIKELIHOOD_H_

#include "PhyloLikelihood.h"
#include "SinglePhyloLikelihood.h"
#include "SubstitutionProcessCollection.h"
#include "../Tree.h"
#include "../TreeTemplate.h"

#include <Bpp/Numeric/AbstractParametrizable.h>

// From SeqLib:
#include <Bpp/Seq/Container/SiteContainer.h>

using namespace std;

namespace bpp
{
namespace newlik
{
/**
 * @brief Partial implementation of the Likelihood interface for multiple processes.
 *
 * This class uses several TreeLikelihoodCalculation instances to compute a the global
 * likelihood of the data set, as well as a collection of SubstitutionProcess.
 * It implements the Function interface and manages the parameters of all substitution processes.
 */
class MultiPhyloLikelihood :
  public virtual PhyloLikelihood,
  public AbstractParametrizable
{
protected:
  std::auto_ptr<const SiteContainer> data_;
  std::auto_ptr<SubstitutionProcessCollection> processColl_;

  /**
   * @brief recursivity of the tree likelihood computations
   *
   * S for simple (default) or D for double (not implemented yet).
   */
  char recursivity_;
  bool computeFirstOrderDerivatives_;
  bool computeSecondOrderDerivatives_;
  bool initialized_;
  bool verbose_;
  size_t nbSites_;
  size_t nbStates_;

  double minusLogLik_;

  /**
   * vector of pointers towards Treelikelihoods, used for the
   * global likelihood.
   */
  
  std::vector<TreeLikelihoodCalculation*> vpTreelik_;

public:
  MultiPhyloLikelihood(SubstitutionProcessCollection* processColl,
                       bool verbose = true,
                       bool patterns = true);

  MultiPhyloLikelihood(
    const SiteContainer& data,
    SubstitutionProcessCollection* processColl,
    bool verbose = true,
    bool patterns = true);

  MultiPhyloLikelihood(const MultiPhyloLikelihood& lik) :
    AbstractParametrizable(lik),
    data_(0),
    processColl_(0),
    recursivity_('S'),
    computeFirstOrderDerivatives_(lik.computeFirstOrderDerivatives_),
    computeSecondOrderDerivatives_(lik.computeSecondOrderDerivatives_),
    initialized_(lik.initialized_),
    verbose_(lik.verbose_),
    nbSites_(lik.nbSites_),
    nbStates_(lik.nbStates_),
    minusLogLik_(lik.minusLogLik_),
    vpTreelik_()
  {
    if (lik.data_.get())
      data_.reset(lik.data_->clone());
    if (lik.processColl_.get())
      processColl_.reset(lik.processColl_->clone());
    for (size_t i = 0; i < lik.vpTreelik_.size(); i++)
    {
      vpTreelik_.push_back(lik.vpTreelik_[i]->clone());
    }
  }

  MultiPhyloLikelihood& operator=(const MultiPhyloLikelihood& lik)
  {
    AbstractParametrizable::operator=(lik);
    if (lik.data_.get())
      data_.reset(lik.data_->clone());
    else
      data_.reset();
    if (lik.processColl_.get())
      processColl_.reset(lik.processColl_->clone());
    else
      processColl_.reset();

    recursivity_ = 'S';
    computeFirstOrderDerivatives_  = lik.computeFirstOrderDerivatives_;
    computeSecondOrderDerivatives_ = lik.computeSecondOrderDerivatives_;
    initialized_                   = lik.initialized_;
    verbose_                       = lik.verbose_;
    nbSites_                       = lik.nbSites_;
    nbStates_                      = lik.nbStates_;
    minusLogLik_                   = lik.minusLogLik_;


    for (size_t i = 0; i < vpTreelik_.size(); i++)
    {
      if (vpTreelik_[i])
        delete vpTreelik_[i];
    }

    vpTreelik_.empty();

    for (size_t i = 0; i < lik.vpTreelik_.size(); i++)
    {
      vpTreelik_.push_back(lik.vpTreelik_[i]->clone());
    }

    return *this;
  }

  /**
   * @brief Abstract class destructor
   *
   */
  virtual ~MultiPhyloLikelihood()
  {
    for (size_t i = 0; i < vpTreelik_.size(); i++)
    {
      if (vpTreelik_[i])
        delete vpTreelik_[i];
    }
    vpTreelik_.empty();
  }

public:
  /**
   * @name The Likelihood interface.
   *
   * @{
   */
  const SiteContainer* getData() const { return data_.get(); }

  const Alphabet* getAlphabet() const { return data_->getAlphabet(); }

  virtual double getLogLikelihood() const = 0;

  virtual Vdouble getLikelihoodForEachSite() const;

  virtual double getDLogLikelihood() const = 0;

  virtual double getDLogLikelihoodForASite(size_t site) const;

  virtual double getD2LogLikelihood() const = 0;
  
  virtual double getD2LogLikelihoodForASite(size_t site) const;

  /*
   * @}
   */

  
  /**
   * @brief The collection
   *
   */

  const SubstitutionProcessCollection* getCollection() const { return processColl_.get(); }

protected:
  virtual void computeDLikelihood_(const std::string& variable) const = 0;

  virtual void computeD2Likelihood_(const std::string& variable) const = 0;

public:
  /**
   * @brief To be defined in inheriting classes.
   *
   */

  /**
   * @brief Get the likelihood for a site.
   *
   * @param site The site index to analyse.
   * @return The likelihood for site <i>site</i>.
   */

  virtual double getLikelihoodForASite(size_t site) const = 0;

  /**
   * @brief Get the first order derivate of the likelihood for a site. 
   *
   * This derivate should have been first computed through
   * getFirstOrderDerivative function.
   *
   * @param site The site index to analyse.
   * @return The first order derivate likelihood for site <i>site</i>.
   */

  virtual double getDLikelihoodForASite(size_t site) const = 0;

  /**
   * @brief Get the second order derivate of the likelihood for a site. 
   *
   * This derivate should have been first computed through
   * getSecondOrderDerivative function.
   *
   * @param site The site index to analyse.
   * @return The second order derivate likelihood for site <i>site</i>.
   */

  virtual double getD2LikelihoodForASite(size_t site) const = 0;

  
  double getLikelihoodForASiteForAProcess(size_t i, size_t p) const
  {
    return vpTreelik_[p]->getLikelihoodForASite(i);
  }

  void computeDLikelihoodForAProcess(std::string& variable, size_t p) const
  {
    vpTreelik_[p]->computeTreeDLikelihood(variable);
  }

  void computeD2LikelihoodForAProcess(std::string& variable, size_t p) const
  {
    vpTreelik_[p]->computeTreeD2Likelihood(variable);
  }

  double getDLikelihoodForASiteForAProcess(size_t i, size_t p) const
  {
    return vpTreelik_[p]->getDLikelihoodForASite(i);
  }

  double getD2LikelihoodForASiteForAProcess(size_t i, size_t p) const
  {
    return vpTreelik_[p]->getD2LikelihoodForASite(i);
  }

  VVdouble getLikelihoodForEachSiteForEachProcess() const;

  virtual VVdouble getPosteriorProbabilitiesForEachSiteForEachProcess() const = 0;

  
  /**
   * @brief Set the dataset for which the likelihood must be evaluated.
   *
   * @param sites The data set to use.
   */
  void setData(const SiteContainer& sites);

  size_t getNumberOfSites() const { return data_->getNumberOfSites(); }
  size_t getNumberOfStates() const { return data_->getAlphabet()->getSize(); }

  size_t getNumberOfSubstitutionProcess() const { return getCollection()->getNumberOfSubstitutionProcess();}

  void enableDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
  void enableFirstOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = yn; }
  void enableSecondOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
  bool enableFirstOrderDerivatives() const { return computeFirstOrderDerivatives_; }
  bool enableSecondOrderDerivatives() const { return computeSecondOrderDerivatives_; }

  bool isInitialized() const { return initialized_; }

  ParameterList getSubstitutionProcessParameters() const;

  ParameterList getSubstitutionModelParameters() const;

  ParameterList getRateDistributionParameters() const;

  ParameterList getRootFrequenciesParameters() const;

  ParameterList getBranchLengthsParameters() const;

  std::vector<const TreeTemplate<Node>* > getTrees() const;
    
  //    ParameterList getTransitionProbabilitiesParameters() const { return process_->getTransitionProbabilitiesParameters(); }
  // TODO: this has to be modified to deal with special cases...

  void fireParameterChanged(const ParameterList& parameters);

  void setParameters(const ParameterList& parameters)   throw (ParameterNotFoundException, ConstraintException);

  double getValue() const throw (Exception);

  /**
   * @name DerivableFirstOrder interface.
   *
   * @{
   */
  double getFirstOrderDerivative(const std::string& variable) const throw (Exception);

  /** @} */

  /**
   * @name DerivableSecondOrder interface.
   *
   * @{
   */
  double getSecondOrderDerivative(const std::string& variable) const throw (Exception);
  double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const throw (Exception) { return 0; } // Not implemented for now.

  /** @} */

  /** @} */
};
}
} // end of namespace bpp.

#endif  // _MULTIPHYLOLIKELIHOOD_H_

