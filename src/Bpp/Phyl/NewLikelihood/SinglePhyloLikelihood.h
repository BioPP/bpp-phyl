//
// File: TSinglePhyloLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Oct 17 17:36:44 2003
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

#ifndef _SINGLEPHYLOLIKELIHOOD_H_
#define _SINGLEPHYLOLIKELIHOOD_H_

#include "../Node.h"
#include "../Tree.h"
#include "../Model/SubstitutionModel.h"
#include "TreeLikelihoodData.h"
#include "ModelIterator.h"
#include "SitePartition.h"

#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/VectorTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/SiteContainer.h>

#include "PhyloLikelihood.h"
#include "SingleRecursiveTreeLikelihoodCalculation.h"

namespace bpp
{
namespace newlik
{
/**
 * @brief The SinglePhyloLikelihood class: phylogenetic likelihood computation with a single process.
 *
 * This class implements likelihood calculation with a single process/tree.
 * It uses a unique TreeLikelihoodCalculation instance, and implements the
 * Function interface, dealing with parameters from the associated SubstitutionProcess.
 */
class SinglePhyloLikelihood :
  public virtual PhyloLikelihood,
  public AbstractParametrizable
{
protected:
  std::auto_ptr<TreeLikelihoodCalculation> tlComp_;
  std::auto_ptr<SubstitutionProcess> process_;
  bool computeFirstOrderDerivatives_;
  bool computeSecondOrderDerivatives_;
  double minusLogLik_;
  bool verbose_;

public:
  SinglePhyloLikelihood(SubstitutionProcess* process, TreeLikelihoodCalculation* tlComp, bool verbose = true);

  SinglePhyloLikelihood(const SinglePhyloLikelihood& lik) :
    AbstractParametrizable(lik),
    tlComp_(0),
    process_(0),
    computeFirstOrderDerivatives_(lik.computeFirstOrderDerivatives_),
    computeSecondOrderDerivatives_(lik.computeSecondOrderDerivatives_),
    minusLogLik_(lik.minusLogLik_),
    verbose_(lik.verbose_)
  {
    if (lik.tlComp_.get()) tlComp_.reset(lik.tlComp_->clone());
    if (lik.process_.get()) process_.reset(lik.process_->clone());
  }

  SinglePhyloLikelihood& operator=(const SinglePhyloLikelihood& lik)
  {
    AbstractParametrizable::operator=(lik);
    if (lik.tlComp_.get()) tlComp_.reset(lik.tlComp_->clone());
    else tlComp_.reset();
    if (lik.process_.get()) process_.reset(lik.process_->clone());
    else process_.reset();
    computeFirstOrderDerivatives_  = lik.computeFirstOrderDerivatives_;
    computeSecondOrderDerivatives_ = lik.computeSecondOrderDerivatives_;
    minusLogLik_                   = lik.minusLogLik_;
    verbose_                       = lik.verbose_;
    return *this;
  }


  virtual ~SinglePhyloLikelihood() {}

  SinglePhyloLikelihood* clone() const { return new SinglePhyloLikelihood(*this); }

public:

  /**
   * @name Handling of data
   *
   * @{
   */
  void setData(const SiteContainer& sites) throw (Exception) {
    tlComp_->setData(sites); //This automatically calls computeTreeLikelihood().
    minusLogLik_ = - tlComp_->getLogLikelihood();
  }
  
  const SiteContainer* getData() const {
    return tlComp_->getData();
  }

  const Alphabet* getAlphabet() const { return tlComp_->getAlphabet(); }

  size_t getNumberOfSites() const { return tlComp_->getNumberOfSites(); }

  /** @} */

  // TODO: need to acount for model classes
  // ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const
  // {
  //  return new ConstNoPartitionBranchModelIterator(modelSet_->getModelForNode(nodeId), nbDistinctSites_);
  // }


  /**
   * @name Handling of substitution process
   *
   * @{
   */

  /**
   * @brief Get the number of states in the alphabet associated to the dataset.
   *
   * @return the number of states in the alphabet associated to the dataset.
   */
  size_t getNumberOfStates() const { return getAlphabet()->getSize(); }

  /**
   * @brief Get the number of model classes.
   *
   * @return The Number of model classes.
   */
  size_t getNumberOfClasses() const { return process_->getNumberOfClasses(); }

  /**
   * @brief Get the tree (topology and branch lengths).
   *
   * @return The tree of this SinglePhyloLikelihood object.
   */
  const Tree& getTree() const { return process_->getTree(); }

  ParameterList getSubstitutionProcessParameters() const { return process_->getParameters(); }

  const SubstitutionProcess& getSubstitutionProcess() const { return *process_; }

  ParameterList getBranchLengthsParameters() const { return process_->getParametrizableTree().getParameters(); }

  ParameterList getRootFrequenciesParameters() const
  {
    return process_->getRootFrequenciesSet() ? process_->getRootFrequenciesSet()->getParameters() : ParameterList();
  }

  ParameterList getRateDistributionParameters() const
  {
    return process_->getRateDistributionParameters();
  }

  ParameterList getSubstitutionModelParameters() const
  {
    return process_->getSubstitutionModelParameters();
  }

  /** @} */


  /**
   * @brief Implements the Function interface.
   *
   * Update the parameter list and call the applyParameters() method.
   * Then compute the likelihoods at each node (computeLikelihood() method)
   * and call the getLogLikelihood() method.
   *
   * If a subset of the whole parameter list is passed to the function,
   * only these parameters are updated and the other remain constant (i.e.
   * equal to their last value).
   *
   * @param parameters The parameter list to pass to the function.
   */
  void setParameters(const ParameterList& parameters) throw (ParameterNotFoundException, ConstraintException);
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


  void computeTreeLikelihood();

protected:
  void fireParameterChanged(const ParameterList& params);
  
  void computeDLikelihood_(const std::string& variable) const;

  void computeD2Likelihood_(const std::string& variable) const;

public:
  /**
   * @return The underlying likelihood data structure.
   */
  virtual TreeLikelihoodData* getLikelihoodData() { return tlComp_->getLikelihoodData(); }

  /**
   * @return The underlying likelihood data structure.
   */
  virtual const TreeLikelihoodData* getLikelihoodData() const { return tlComp_->getLikelihoodData(); }

  void enableDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
  void enableFirstOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = yn; }
  void enableSecondOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
  bool enableFirstOrderDerivatives() const { return computeFirstOrderDerivatives_; }
  bool enableSecondOrderDerivatives() const { return computeSecondOrderDerivatives_; }
  bool isInitialized() const { return tlComp_->isInitialized(); }

  //    ParameterList getTransitionProbabilitiesParameters() const { return process_->getTransitionProbabilitiesParameters(); }
  // TODO: this has to be modified to deal with special cases...
  ParameterList getDerivableParameters() const { return getBranchLengthsParameters(); }
  ParameterList getNonDerivableParameters() const { return getSubstitutionProcessParameters(); }


  double getLogLikelihood() const {
    return tlComp_->getLogLikelihood();
  }

  double getLikelihoodForASite(size_t siteIndex) const {
    return tlComp_->getLikelihoodForASite(siteIndex);
  }

  /**
   * @brief Get the likelihood for each site.
   *
   * @return A 1d vector with all likelihoods for each site.
   */
  Vdouble getLikelihoodForEachSite() const;

  /**
   * @brief Get the likelihood for each site and for each state.
   *
   * @return A 2d vector with all likelihoods for each site and for each state.
   */
  VVdouble getLikelihoodForEachSiteForEachState() const;

  /**
   * @brief Get the likelihood for each site and each model class.
   *
   * @return A two-dimension vector with all likelihoods:
   * <code>V[i][j] =</code> likelihood of site i and model class j.
   */
  VVdouble getLikelihoodForEachSiteForEachClass() const;

  /**
   * @brief Get the likelihood for each site and each model class and each state.
   *
   * @return A three-dimension vector with all likelihoods:
   * <code>V[i][j][k} =</code> likelihood of site i and model class j and state k.
   */
  VVVdouble getLikelihoodForEachSiteForEachClassForEachState() const;

  /** @} */

  /**
   * @brief Get the posterior model class (the one with maximum posterior
   * probability) for each site.
   *
   * @return A vector with all model classes indexes.
   */
  std::vector<size_t> getClassWithMaxPostProbOfEachSite() const;

  VVdouble getPosteriorProbabilitiesOfEachClass() const;

  /**
   * @brief Get the index (used for inner computations) of a given site (original alignment column).
   *
   * @param site An alignment position.
   * @return The site index corresponding to the given input alignment position.
   */
  size_t getSiteIndex(size_t site) const throw (IndexOutOfBoundsException) {
    return tlComp_->getSiteIndex(site);
  }

  /**
   * @name Iterators
   * @{
   */
  // TODO jdutheil on 21/04/13: need to account for model classes!
  // virtual ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const = 0;

  // jdutheil on 21/04/13: I think we will drop this type of iterator, which were never used before and are difficult to implement in the new framework...
  // virtual ConstSiteModelIterator* getNewSiteModelIterator(size_t siteIndex) const = 0;

  // 19/07/13 jdutheil: copied for AbstractTreeLikelihood:

  // TODO jdutheil on 08.04.13 we drop model iterators for now
  // ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const {
  //  return new ConstNoPartitionBranchModelIterator(model_.get(), sitePartition_->getNumberOfPatternsForPartition(0));
  // }

  // ConstSiteModelIterator* getNewSiteModelIterator(size_t siteIndex) const {
  //  return new ConstHomogeneousSiteModelIterator(*pTree_, model_.get());
  // }


  /* @} */

};
} // end of namespace newlik.
} // end of namespace bpp.

#endif  // _SINGLEPHYLOLIKELIHOOD_H_

