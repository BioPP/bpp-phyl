//
// File: PhyloLikelihood.h
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

#ifndef _PHYLOLIKELIHOOD_H_
#define _PHYLOLIKELIHOOD_H_

#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Function/Functions.h>

namespace bpp
{

    /**
     * @brief The PhyloLikelihood interface, for phylogenetic likelihood.
     *
     * This interface defines the common methods needed to compute a likelihood
     * from a sequence alignement, usually involving one or more phylogenetic trees.
     */
    class PhyloLikelihood:
      public virtual DerivableSecondOrder
    {
    public:
      PhyloLikelihood() {}
      virtual ~PhyloLikelihood() {}

      PhyloLikelihood* clone() const = 0;

    public:

      /**
       *
       * @name The data functions
       *
       * @{
       */
      
      /**
       * @return initialize the likelihood function.
       */
      
      virtual void initialize() = 0;
      
      /**
       * @return 'true' is the likelihood function has been initialized.
       */
      virtual bool isInitialized() const = 0;
      
      /**
       * @}
       */

      /**
       * @brief set it arrays should be computed in log.
       *
       */

      virtual void setUseLog(bool useLog) = 0;

      /**
       * @name The likelihood functions.
       *
       * @{
       */
      
      /**
       * @brief update the likelihood to get ready for computation
       *
       */

      virtual void updateLikelihood() const = 0;

      /**
       * @brief compute the likelihood
       *
       */

      virtual void computeLikelihood() const = 0;

      /**
       * @brief Get the logarithm of the likelihood for the whole dataset.
       *
       * @return The logarithm of the likelihood of the dataset.
       */
      virtual double getLogLikelihood() const = 0;
      
      /**
       * @brief Get the derivates of the LogLikelihood.
       *
       */

      virtual double getDLogLikelihood() const = 0;

      virtual double getD2LogLikelihood() const = 0;

      virtual void computeDLogLikelihood_(const std::string& variable) const = 0;
      
      virtual void computeD2LogLikelihood_(const std::string& variable) const = 0;


      /** @} */

      /**
       * @name Retrieve some particular independent parameters subsets.
       *
       * @{
       */
    
      /**
       * @brief Get the independent branch lengths parameters.
       *
       * @return A ParameterList with all branch lengths.
       */

       virtual ParameterList getBranchLengthParameters() const = 0;
    
       /**
        * @brief Get the independent parameters associated to substitution model(s).
        *
        * @return A ParameterList.
        */

       virtual ParameterList getSubstitutionModelParameters() const = 0;

       /**
        * @brief Get the independent parameters associated to the rate distribution(s).
        *
        * @return A ParameterList.
        */

       virtual ParameterList getRateDistributionParameters() const = 0;

       /**
        * @brief Get the independent parameters associated to the root frequencies(s).
        *
        * @return A ParameterList.
        */

       virtual ParameterList getRootFrequenciesParameters() const = 0;

      /**
       * @brief All independent derivable parameters.
       *
       * Usually, this contains all branch lengths parameters.
       *
       * @return A ParameterList.
       */

      virtual ParameterList getDerivableParameters() const = 0;

      /**
       * @brief All independent non derivable parameters.
       *
       * Usually, this contains all substitution model parameters and rate distribution.
       *
       * @return A ParameterList.
       */

      virtual ParameterList getNonDerivableParameters() const = 0;

      /** @} */

      /**
       * @brief Tell if derivatives must be computed.
       *
       * This methods calls the enableFirstOrderDerivatives and enableSecondOrderDerivatives.
       *
       * @param yn Yes or no.
       */
      virtual void enableDerivatives(bool yn) = 0;

    };


    class AbstractPhyloLikelihood :
      public virtual PhyloLikelihood
    {
    protected:
      
      /**
       * @brief the value
       *
       **/
      
      mutable double minusLogLik_;

      /**
       * @brief sey if derivatives should be computed
       *
       */
      
      bool computeFirstOrderDerivatives_;
      bool computeSecondOrderDerivatives_;

      // say if the Likelihoods should be recomputed
      
      mutable bool computeLikelihoods_;

      // say if initialized
      
      mutable bool initialized_;

    public:
      AbstractPhyloLikelihood() :
        minusLogLik_(0),
        computeFirstOrderDerivatives_(true),
        computeSecondOrderDerivatives_(true),
        computeLikelihoods_(true),
        initialized_(false)
      {
      }

      AbstractPhyloLikelihood(const AbstractPhyloLikelihood& asd) :
        minusLogLik_(asd.minusLogLik_),
        computeFirstOrderDerivatives_(asd.computeFirstOrderDerivatives_),
        computeSecondOrderDerivatives_(asd.computeSecondOrderDerivatives_),
        computeLikelihoods_(asd.computeLikelihoods_),
        initialized_(asd.initialized_)
      {
      }
      
      AbstractPhyloLikelihood& operator=(const AbstractPhyloLikelihood& asd)
      {
        minusLogLik_                   = asd.minusLogLik_;
        computeFirstOrderDerivatives_  = asd.computeFirstOrderDerivatives_;
        computeSecondOrderDerivatives_ = asd.computeSecondOrderDerivatives_;

        computeLikelihoods_ = asd.computeLikelihoods_;
        initialized_ = asd.initialized_;
        
        return *this;
      }

      virtual ~AbstractPhyloLikelihood() {}

      AbstractPhyloLikelihood* clone() const = 0;

    public:
      /**
       * @brief Sets the computeLikelihoods_ to true.
       *
       */
      
      void update()
      {
        computeLikelihoods_ = true;
      }

      virtual void initialize()
      {
        initialized_=true;
      }

    public:

      void setParameters(const ParameterList& parameters) throw (ParameterNotFoundException, ConstraintException)
      {
        setParametersValues(parameters);
      }

      virtual void enableDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
      virtual void enableFirstOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = yn; }
      virtual void enableSecondOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
      bool enableFirstOrderDerivatives() const { return computeFirstOrderDerivatives_; }
      bool enableSecondOrderDerivatives() const { return computeSecondOrderDerivatives_; }

      bool isInitialized() const { return initialized_; }

      /*
       * @brief return the value, ie -loglikelihood
       *
       * !!! check on computeLikelihoods_ is not done here.
       *
       */
      
      double getValue() const throw (Exception)
      {
        if (!isInitialized())
          throw Exception("AbstractPhyloLikelihood::getValue(). Instance is not initialized.");

        minusLogLik_=-getLogLikelihood();
        
        return minusLogLik_;
      }


    };

} //end of namespace bpp.

#endif  //_PHYLOLIKELIHOOD_H_

