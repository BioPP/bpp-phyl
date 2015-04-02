//
// File: SumOfDataPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: jeudi 11 juillet 2013, à 14h 05
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

#ifndef _SUMOFDATAPHYLOLIKELIHOOD_H_
#define _SUMOFDATAPHYLOLIKELIHOOD_H_

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/SiteContainer.h>

#include "MultiDataPhyloLikelihood.h"

namespace bpp
{
  namespace newlik
  {

    /**
     * @brief The SumOfDataPhyloLikelihood class, for phylogenetic
     * likelihood on several independent data.
     *
     */
    
    class SumOfDataPhyloLikelihood:
      public virtual MultiDataPhyloLikelihood,
      public AbstractParametrizable    
    {
    protected:
      std::map<size_t, AbstractSingleDataPhyloLikelihood*>  mSDP_;
      bool computeFirstOrderDerivatives_;
      bool computeSecondOrderDerivatives_;
      mutable double minusLogLik_;
      
    public:
      SumOfDataPhyloLikelihood();


      /*
       * @brief Build from a map of AbstractSingleDataPhyloLikelihood.
       *
       * BEWARE : Will own the AbstractSingleDataPhyloLikelihood objects.
       *
       */
      
      SumOfDataPhyloLikelihood(std::map<size_t, SingleDataPhyloLikelihood*>& mSDP);

      ~SumOfDataPhyloLikelihood();

      SumOfDataPhyloLikelihood* clone() const
      {
        return new SumOfDataPhyloLikelihood(*this);
      }

      SumOfDataPhyloLikelihood(const SumOfDataPhyloLikelihood& sd);
        
      SumOfDataPhyloLikelihood& operator=(const SumOfDataPhyloLikelihood& sd);

    public:

      /**
       *
       * @name The data functions
       *
       * @{
       */
      
      /**
       * @brief Set the dataset for which the likelihood must be evaluated.
       *
       * @param nPhyl The number of the Likelihood.
       * @param sites The data set to use.
       */

      void setData(size_t nPhyl, const SiteContainer* sites)
      {
        if (mSDP_.find(nPhyl)!=mSDP_.end())
          mSDP_[nPhyl]->setData(*sites);
      }
      
    
      /**
       * @brief Get the dataset for which the likelihood must be evaluated.
       *
       * @return A pointer toward the site container where the sequences are stored.
       */

      const SiteContainer* getData(size_t nPhyl) const
      {
        if (mSDP_.find(nPhyl)!=mSDP_.end())
          return mSDP_.find(nPhyl)->second->getData();
        else
          return 0;
      }
      
      std::vector<size_t> getNumbersOfSingleDataPhyloLikelihoods() const;
      
      /**
       * @}
       */

      /**
       *
       * @name The single data Phylolikelihood storage.
       *
       * @{
       */

      /*
       * @brief Adds a SingleDataPhyloLikelihood if it is an
       * AbstractSingleDataPhyloLikelihood object.
       *
       * Gets ownership of the SingleDataPhyloLikelihood object.
       *
       * @param nPhyl the number of the phylolikelihood object.
       * @param SDP a pointer to the phylolikelihood object.
       */
      
      void addSingleDataPhylolikelihood(size_t nPhyl, SingleDataPhyloLikelihood* SDP);

      const SingleDataPhyloLikelihood* getSingleDataPhylolikelihood(size_t nPhyl) const
      {
        if (mSDP_.find(nPhyl)!=mSDP_.end())
          return mSDP_.find(nPhyl)->second;
        else
          return 0;
      }
      
      SingleDataPhyloLikelihood* getSingleDataPhylolikelihood(size_t nPhyl)
      {
        if (mSDP_.find(nPhyl)!=mSDP_.end())
          return mSDP_.find(nPhyl)->second;
        else
          return 0;
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
      
      bool isInitialized() const 
      {
        std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it;
        
        for (it=mSDP_.begin(); it != mSDP_.end(); it++)
          if (! it->second->isInitialized())
            return false;
        return true;
      }

      /**
       * @name The likelihood functions.
       *
       */

      /**
       * @brief Get the logarithm of the likelihood for the whole dataset.
       *
       * @return The logarithm of the likelihood of the dataset.
       */

      double getLogLikelihood() const
      {
        double x=0;
        for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
          x += it->second->getLogLikelihood();
        
        minusLogLik_=-x;
        
        return x;
      }

      double getValue() const throw (Exception)
      {
        if (!isInitialized())
          throw Exception("SingleProcessPhyloLikelihood::getValue(). Instance is not initialized.");
        minusLogLik_=0;
        for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
          minusLogLik_ += it->second->getValue();
        
        return minusLogLik_;
      }
        
      /**
       * @brief Get the derivates of the LogLikelihood.
       *
       */

      double getDLogLikelihood() const
      {
        double x=0;
        for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
          x += it->second->getDLogLikelihood();
        return x;
      }


      double getD2LogLikelihood() const
      {
        double x=0;
        for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
          x += it->second->getD2LogLikelihood();
        return x;
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

    protected:
      
      void fireParameterChanged(const ParameterList& params)
      {
        for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
          it->second->matchParametersValues(params);
        
        getValue();
      }
      
      /**
       * @name Retrieve some particular parameters subsets.
       *
       * @{
       */
    
      /**
       * @brief Get the branch lengths parameters.
       *
       * @return A ParameterList with all branch lengths.
       */

      ParameterList getBranchLengthParameters() const;
      
      /**
      * @brief Get the parameters associated to substitution model(s).
      * 
      * @return A ParameterList.
      */

      ParameterList getSubstitutionModelParameters() const;

      /**
       * @brief Get the parameters associated to the rate distribution(s).
       *
       * @return A ParameterList.
       */

      ParameterList getRateDistributionParameters() const;

      /**
       * @brief Get the parameters associated to the root frequencies(s).
       *
       * @return A ParameterList.
       */

      ParameterList getRootFrequenciesParameters() const;

      /**
       * @brief All derivable parameters.
       *
       * Usually, this contains all branch lengths parameters.
       *
       * @return A ParameterList.
       */

      ParameterList getDerivableParameters() const;

      /**
       * @brief All non derivable parameters.
       *
       * Usually, this contains all substitution model parameters and rate distribution.
       *
       * @return A ParameterList.
       */

      ParameterList getNonDerivableParameters() const;

      /** @} */

      /**
       * @brief Tell if derivatives must be computed.
       *
       * This methods calls the enableFirstOrderDerivatives and enableSecondOrderDerivatives.
       *
       * @param yn Yes or no.
       */
      
      void enableDerivatives(bool yn);

      void enableFirstOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = yn; }
      void enableSecondOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }

      bool enableFirstOrderDerivatives() const { return computeFirstOrderDerivatives_; }
      
      bool enableSecondOrderDerivatives() const { return computeSecondOrderDerivatives_;
      }
      
    };

  }; //end of namespace newlik.
} //end of namespace bpp.

#endif  //_MULTIDATAPHYLOLIKELIHOOD_H_

