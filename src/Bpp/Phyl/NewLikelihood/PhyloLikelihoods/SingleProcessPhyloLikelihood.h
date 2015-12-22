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

#include "../../Tree/Node.h"
#include "../../Tree/Tree.h"
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
#include "PhyloLikelihood.h"

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
      mutable std::auto_ptr<LikelihoodTreeCalculation> tlComp_;
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
      tlComp_(0),
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
      void setData(const SiteContainer& sites, size_t nData = 0)
      {
        AbstractSingleDataPhyloLikelihood::setData(sites, nData);

        update();
                
        tlComp_->setData(sites);
      }

      /**
       * @brief return a pointer to the compressed data. 
       *
       */
      
      const SiteContainer* getData() const {
        return tlComp_->getData();
      }

      const Alphabet* getAlphabet() const {
        return tlComp_->getAlphabet();
      }

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
      const Tree& getTree() const { return process_->getTree(); }

      const SubstitutionProcess& getSubstitutionProcess() const { return *process_; }

      size_t getSubstitutionProcessNumber() const { return nProc_; }

      /** @} */

      /**
       * @name Handling of parameters
       *
       * @{
       */

      bool hasDerivableParameter(const std::string& name) const
      {
        return process_->hasDerivableParameter(name);
      }
      
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

      //    ParameterList getTransitionProbabilitiesParameters() const { return process_->getTransitionProbabilitiesParameters(); }
      // TODO: this has to be modified to deal with special cases...
      ParameterList getDerivableParameters() const {
        // patch, to be fixed properly later
        return ParameterList();

        return getBranchLengthParameters();
      }

      ParameterList getNonDerivableParameters() const;

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

      double getFirstOrderDerivative(const std::string& variable) const throw (Exception)
      {
        if (!hasParameter(variable))
          throw ParameterNotFoundException("SingleProcessPhyloLikelihood::getFirstOrderDerivative().", variable);
//        if (!hasDerivableParameter(variable))
        {
          throw Exception("SingleProcessPhyloLikelihood::getFirstOrderDerivative : Derivative is not implemented for " + variable + " parameter.");
        }
        
        computeDLogLikelihood_(variable);
        return -getDLogLikelihood();
      }

      double getSecondOrderDerivative(const std::string& variable) const throw (Exception)
      {
        if (!hasParameter(variable))
          throw ParameterNotFoundException("SingleProcessPhyloLikelihood::getSecondOrderDerivative().", variable);
//        if (!hasDerivableParameter(variable))
        {
          throw Exception("SingleProcessPhyloLikelihood::getSecondOrderDerivative : Derivative is not implemented for " + variable + " parameter.");
        }
        computeD2LogLikelihood_(variable);
        return -getD2LogLikelihood();
      }

      double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const throw (Exception) { return 0; } // Not implemented for now.

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

      double getDLogLikelihood() const
      {
        return tlComp_->getDLogLikelihood();
      }

      double getD2LogLikelihood() const
      {
        return tlComp_->getD2LogLikelihood();
      }

      double getLikelihoodForASite(size_t siteIndex) const
      {
        updateLikelihood();
        computeLikelihood();
        return tlComp_->getLikelihoodForASite(siteIndex);
      }

      double getLogLikelihoodForASite(size_t siteIndex) const
      {
        updateLikelihood();
        computeLikelihood();
        return tlComp_->getLogLikelihoodForASite(siteIndex);
      }

      double getDLogLikelihoodForASite(size_t siteIndex) const {
        return tlComp_->getDLogLikelihoodForASite(siteIndex);
      }

      double getD2LogLikelihoodForASite(size_t siteIndex) const {
        return tlComp_->getD2LogLikelihoodForASite(siteIndex);
      }

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
       * @brief Get the index (used for inner computations) of a given site (original alignment column).
       *
       * @param site An alignment position.
       * @return The site index corresponding to the given input alignment position.
       */
      size_t getSiteIndex(size_t site) const throw (IndexOutOfBoundsException) {
        return tlComp_->getSiteIndex(site);
      }
      
      /**
       * Utilities
       *
       */
      
      VVdouble getPosteriorProbabilitiesOfEachClass() const;
      
      /**
       * @brief Get the posterior model class (the one with maximum posterior
       * probability) for each site.
       *
       * @return A vector with all model classes indexes.
       */
      std::vector<size_t> getClassWithMaxPostProbOfEachSite() const;
      
      Vdouble getPosteriorRateOfEachSite() const;

      /* @} */

    };
} // end of namespace bpp.

#endif  // _SINGLEPROCESSPHYLOLIKELIHOOD_H_

