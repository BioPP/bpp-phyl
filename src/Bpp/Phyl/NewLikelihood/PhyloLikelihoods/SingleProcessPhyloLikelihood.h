//
// File: SingleProcessPhyloLikelihood.h
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

#ifndef _SINGLEPROCESSPHYLOLIKELIHOOD_H_
#define _SINGLEPROCESSPHYLOLIKELIHOOD_H_

#include "../../Tree/PhyloNode.h"
#include "../../Tree/PhyloTree.h"
#include "../../Model/SubstitutionModel.h"
#include "../ModelIterator.h"
#include "../SitePartition.h"

#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/VectorTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>

#include "SingleDataPhyloLikelihood.h"
#include "../RecursiveLikelihoodTreeCalculation.h"
#include "AbstractPhyloLikelihood.h"

namespace bpp
{
/**
 * @brief The SingleProcessPhyloLikelihood class: phylogenetic
 * likelihood computation with a single process.
 *
 * This class implements likelihood calculation with a single
 * process/tree. It uses a unique LikelihoodTreeCalculation instance,
 * and implements the Function interface, dealing with parameters from
 * the associated SubstitutionProcess.
 */

  class SingleProcessPhyloLikelihood :
    public AbstractSingleDataPhyloLikelihood,
    public AbstractParametrizable
  { 
  protected:
    mutable std::unique_ptr<LikelihoodTreeCalculation> tlComp_;
    SubstitutionProcess* process_;

    /**
     * @brief the Substitution Process number
     *
     **/

    size_t nProc_;

  public:
    SingleProcessPhyloLikelihood(SubstitutionProcess* process, LikelihoodTreeCalculation* tlComp, size_t nProc = 0, size_t nData = 0);

    SingleProcessPhyloLikelihood(const SingleProcessPhyloLikelihood& lik) :
      AbstractPhyloLikelihood(lik),
      AbstractAlignedPhyloLikelihood(lik),
      AbstractSingleDataPhyloLikelihood(lik),
      AbstractParametrizable(lik),
      tlComp_(),
      process_(lik.process_),
      nProc_(lik.nProc_)
    {
      if (lik.tlComp_.get()) tlComp_.reset(lik.tlComp_->clone());
    }

    SingleProcessPhyloLikelihood& operator=(const SingleProcessPhyloLikelihood& lik)
    {
      AbstractSingleDataPhyloLikelihood::operator=(lik);

      AbstractParametrizable::operator=(lik);

      if (lik.tlComp_.get()) tlComp_.reset(lik.tlComp_->clone());
      else tlComp_.reset();

      process_=lik.process_;
      nProc_=lik.nProc_;
    
      return *this;
    }

    virtual ~SingleProcessPhyloLikelihood() {}

    SingleProcessPhyloLikelihood* clone() const { return new SingleProcessPhyloLikelihood(*this); }

  public:

    /**
     * @name Handling of data
     *
     * @{
     */
    void setData(const AlignedValuesContainer& sites, size_t nData = 0)
    {
      AbstractSingleDataPhyloLikelihood::setData(sites, nData);

      update();
                
      tlComp_->setData(sites);
    }

    /**
     * @brief return a pointer to the compressed data. 
     *
     */
      
    const AlignedValuesContainer* getData() const {
      return tlComp_->getData();
    }

    const Alphabet* getAlphabet() const {
      return tlComp_->getAlphabet();
    }

    /** @} */

    /**
     * @name Handling of substitution process
     *
     * @{
     */

    /**
     * @brief Get the number of model classes.
     *
     * @return The Number of model classes.
     */
    size_t getNumberOfClasses() const { return process_->getNumberOfClasses(); }

    /**
     * @brief Get the tree (topology and branch lengths).
     *
     * @return The tree of this SingleProcessPhyloLikelihood object.
     */
    const ParametrizablePhyloTree& getTree() const { return process_->getParametrizablePhyloTree(); }

    const SubstitutionProcess& getSubstitutionProcess() const { return *process_; }

    size_t getSubstitutionProcessNumber() const { return nProc_; }

    /** @} */

    /**
     * @name Handling of parameters
     *
     * @{
     */

    ParameterList getBranchLengthParameters() const
    {
      return process_->getBranchLengthParameters(true);
    }

    ParameterList getRootFrequenciesParameters() const
    {
      return process_->getRootFrequenciesParameters(true);
    }

    ParameterList getRateDistributionParameters() const
    {
      return process_->getRateDistributionParameters(true);
    }

    ParameterList getSubstitutionModelParameters() const
    {
      return process_->getSubstitutionModelParameters(true);
    }

    ParameterList getNonDerivableParameters() const
    {
      ParameterList pl=getSubstitutionModelParameters();
      
      pl.includeParameters(getRateDistributionParameters());
      pl.includeParameters(getRootFrequenciesParameters());

      return pl;
    }

        
    /** @} */


    void updateLikelihood() const
    {
      if (computeLikelihoods_)
        tlComp_->updateLikelihood();
    }
      
    void computeLikelihood() const
    {
      if (computeLikelihoods_)
      {
        tlComp_->computeTreeLikelihood();
        computeLikelihoods_=false;
      }
    }

    bool isInitialized() const {
      return tlComp_->isInitialized();
    }

  protected:
    void fireParameterChanged(const ParameterList& params);
  
    void computeDLogLikelihood_(const std::string& variable) const;

    void computeD2LogLikelihood_(const std::string& variable) const;

  public:

    double getFirstOrderDerivative(const std::string& variable) const 
    {
      if (dValues_.find(variable)==dValues_.end())
        computeDLogLikelihood_(variable);
      
      if (dValues_.find(variable)==dValues_.end() || std::isnan(dValues_[variable]))
        dValues_[variable]=-getDLogLikelihood(variable);
      
      return dValues_[variable];
    }

    double getSecondOrderDerivative(const std::string& variable) const 
    {
      if (d2Values_.find(variable)==d2Values_.end())
        computeD2LogLikelihood_(variable);
      
      if (d2Values_.find(variable)==d2Values_.end() || std::isnan(d2Values_[variable]))
        d2Values_[variable]=-getD2LogLikelihood(variable);
      
      return d2Values_[variable];
    }
    
    double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const { return 0; } // Not implemented for now.
    
    /**
     * @return The underlying likelihood computation structure.
     */
      
    LikelihoodTreeCalculation* getLikelihoodCalculation() { return tlComp_.get(); }

    /**
     * @return The underlying likelihood data structure.
     */
    virtual LikelihoodTree* getLikelihoodData() { return &tlComp_->getLikelihoodData(); }

    /**
     * @return The underlying likelihood data structure.
     */
    virtual const LikelihoodTree* getLikelihoodData() const { return &tlComp_->getLikelihoodData(); }

    /**
     * @brief set it arrays should be computed in log.
     *
     */

    void setUseLog(bool useLog)
    {
      tlComp_->setAllUseLog(useLog);
    }
      
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
     */
      

    double getLogLikelihood() const
    {
      updateLikelihood();
      computeLikelihood();
      return tlComp_->getLogLikelihood();
    }

    double getDLogLikelihood(const std::string& variable) const
    {
      return tlComp_->getDLogLikelihood();
    }

    double getD2LogLikelihood(const std::string& variable) const
    {
      return tlComp_->getD2LogLikelihood();
    }

    double getLikelihoodForASite(size_t site) const
    {
      updateLikelihood();
      computeLikelihood();
      return tlComp_->getLikelihoodForASite(site);
    }

    double getLogLikelihoodForASite(size_t site) const
    {
      updateLikelihood();
      computeLikelihood();
      return tlComp_->getLogLikelihoodForASite(site);
    }

    double getDLogLikelihoodForASite(const std::string& variable, size_t site) const {
      if (dValues_.find(variable)!=dValues_.end())
        return tlComp_->getDLogLikelihoodForASite(site);
      else
        return 0;
    }

    double getD2LogLikelihoodForASite(const std::string& variable, size_t site) const {
      if (dValues_.find(variable)!=dValues_.end())
        return tlComp_->getD2LogLikelihoodForASite(site);
      else
        return 0;
    }

    /**
     * @brief Get the likelihood for each site and for each state.
     *
     * @return A 2d vector with all likelihoods for each site and for each state.
     */
    VVdouble getLikelihoodPerSitePerState() const;

    /**
     * @brief Get the likelihood for each site and each model class.
     *
     * @return A two-dimension vector with all likelihoods:
     * <code>V[i][j] =</code> likelihood of site i and model class j.
     */
    VVdouble getLikelihoodPerSitePerClass() const;

    /**
     * @brief Get the likelihood for a site and each model class.
     * @param i the index of the site
     *
     * @return A  vector with all likelihoods:
     * <code>V[j] =</code> likelihood of site i and model class j.
     */

    Vdouble getLikelihoodForSitePerClass(size_t i) const;

    /**
     * @brief Get the likelihood for each site and each model class and each state.
     *
     * @return A three-dimension vector with all likelihoods:
     * <code>V[i][j][k} =</code> likelihood of site i and model class j and state k.
     */
    VVVdouble getLikelihoodPerSitePerClassPerState() const;
      
    /** @} */
      
    /**
     * Utilities
     *
     */

    /*
     * @brief Compute and return the Posterior Probabilities Of Rate
     *        Classes on all sites (array site X classes=)
     */
    
    VVdouble getPosteriorProbabilitiesPerClass() const;

    /*
     * @brief Compute and return the Posterior Probabilities Of Rate
     *        Classes on a site
     * @param i the index of the site
     */
    
    Vdouble getPosteriorProbabilitiesForSitePerClass(size_t i) const;

    /**
     * @brief Get the posterior model class (the one with maximum posterior
     * probability) for each site.
     *
     * @return A vector with all model classes indexes.
     */

    std::vector<size_t> getClassWithMaxPostProbPerSite() const;
      
    Vdouble getPosteriorRatePerSite() const;

    /* @} */

  };
} // end of namespace bpp.

#endif  // _SINGLEPROCESSPHYLOLIKELIHOOD_H_
