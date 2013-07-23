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
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/VectorTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/SiteContainer.h>

namespace bpp
{
  namespace newlik
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
       * @name The data functions
       *
       * @{
       */
      
      /**
       * @brief Set the dataset for which the likelihood must be evaluated.
       *
       * @param sites The data set to use.
       */
      virtual void setData(const SiteContainer& sites) = 0;
    
      /**
       * @brief Get the dataset for which the likelihood must be evaluated.
       *
       * @return A pointer toward the site container where the sequences are stored.
       */
      virtual const SiteContainer* getData() const = 0;

      /**
       * @return 'true' is the likelihood function has been initialized.
       */
      virtual bool isInitialized() const = 0;

      /**
       * @brief Get the number of sites in the dataset.
       *
       * @return the number of sites in the dataset.
       */
      virtual size_t getNumberOfSites() const = 0;

      /**
       * @brief Get the alphabet associated to the dataset.
       *
       * @return the alphabet associated to the dataset.
       */    
      virtual const Alphabet* getAlphabet() const = 0;
 
    
      /**
       * @}
       */
      
      /**
       * @name The likelihood functions.
       *
       * @{
       */

      /**
       * @brief Get the logarithm of the likelihood for the whole dataset.
       *
       * @return The logarithm of the likelihood of the dataset.
       */
      virtual double getLogLikelihood() const = 0;
 
      /**
       * @brief Get the likelihood for a site.
       *
       * @param site The site index to analyse.
       * @return The likelihood for site <i>site</i>.
       */
      virtual double getLikelihoodForASite(size_t site) const = 0;

      /**
       * @brief Get the log likelihood for a site.
       *
       * @param site The site index to analyse.
       * @return The log likelihood for site <i>site</i>.
       */
      
      virtual double getLogLikelihoodForASite(size_t site) const
      {
        return log(getLikelihoodForASite(site));
      }


      /**
       * @brief Get the likelihood for each site.
       *
       * @return A vector with all likelihoods for each site.
       */
      virtual Vdouble getLikelihoodForEachSite() const = 0;

      /** @} */

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

      virtual ParameterList getBranchLengthsParameters() const = 0;
    
      /**
       * @brief Get the parameters associated to substitution process(s).
       *
       * @return A ParameterList.
       */

      virtual ParameterList getSubstitutionProcessParameters() const = 0;

      /**
       * @brief Get the parameters associated to substitution model(s).
       *
       * @return A ParameterList.
       */

      virtual ParameterList getSubstitutionModelParameters() const = 0;

      /**
       * @brief Get the parameters associated to the rate distribution(s).
       *
       * @return A ParameterList.
       */

      virtual ParameterList getRateDistributionParameters() const = 0;

      /**
       * @brief Get the parameters associated to the root frequencies(s).
       *
       * @return A ParameterList.
       */

      virtual ParameterList getRootFrequenciesParameters() const = 0;

      /**
       * @brief All derivable parameters.
       *
       * Usually, this contains all branch lengths parameters.
       *
       * @return A ParameterList.
       */

      virtual ParameterList getDerivableParameters() const = 0;

      /**
       * @brief All non derivable parameters.
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

  } //end of namespace newlik.
} //end of namespace bpp.

#endif  //_PHYLOLIKELIHOOD_H_

