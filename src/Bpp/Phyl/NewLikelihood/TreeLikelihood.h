//
// File: TreeLikelihood.h
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

#ifndef _NEWLIK_TREELIKELIHOOD_H_
#define _NEWLIK_TREELIKELIHOOD_H_

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

namespace bpp
{
namespace newlik
{

/**
 * @brief The TreeLikelihood interface.
 *
 * This interface defines the methods needed for computing the likelihood
 * of a phylogenetic tree, given a dataset.
 */ 
class TreeLikelihood:
  public virtual DerivableSecondOrder
{
  public:
    TreeLikelihood() {}
    virtual ~TreeLikelihood() {}

#ifndef NO_VIRTUAL_COV
    TreeLikelihood* clone() const = 0;
#endif

  public:

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
     * @brief Get the number of states in the alphabet associated to the dataset.
     *
     * @return the number of states in the alphabet associated to the dataset.
     */    
    virtual size_t getNumberOfStates() const = 0;
 
    /**
		 * @brief Get the number of model classes.
		 *
		 * @return The Number of model classes.
		 */
		virtual size_t getNumberOfClasses() const = 0;

    /**
     * @brief Get the tree (topology and branch lengths).
     *
     * @return The tree of this TreeLikelihood object.
      */
    virtual const Tree& getTree() const = 0;
   
    /**
     * @brief Get the alphabet associated to the dataset.
     *
     * @return the alphabet associated to the dataset.
     */    
    virtual const Alphabet* getAlphabet() const = 0;
 
    
    /**
     * @return The underlying likelihood data structure.
     */
    virtual TreeLikelihoodData* getLikelihoodData() = 0;

    /**
     * @return The underlying likelihood data structure.
     */
    virtual const TreeLikelihoodData* getLikelihoodData() const = 0;

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
     * @brief Get the likelihood for a site and for a state.
     *
     * @param site The site index to analyse.
     * @param state The state to consider.
     * @return The logarithm of the likelihood for site <i>site</i> and state <i>state</i>.
     */
    virtual double getLikelihoodForASiteForAState(size_t site, int state) const = 0;

    /**
		 * @brief Get the logarithm of the likelihood for a site knowing its model class.
		 *
		 * @param site      The site index.
		 * @param rateClass The model class index.
		 * @return The logarithm of the likelihood for the specified site and model class.
		 */
		virtual double getLikelihoodForASiteForAClass(size_t site, size_t modelClass) const = 0;
	
		/**
		 * @brief Get the likelihood for a site knowing its model class and its ancestral state.
		 *
		 * @param site      The site index.
		 * @param modelClass The model class index.
		 * @param state     The ancestral state.
		 * @return The logarithm of the likelihood for the specified site and model class and ancestral state..
		 */
		virtual double getLikelihoodForASiteForAClassForAState(size_t site, size_t modelClass, int state) const = 0;
 
    /**
     * @brief Get the likelihood for each site.
     *
     * @return A vector with all likelihoods for each site.
     */
    virtual Vdouble getLikelihoodForEachSite() const = 0;

    /**
     * @brief Get the likelihood for each site and for each state.
     *
     * @return A 2d vector with all log likelihoods for each site and for each state.
     */
    virtual VVdouble getLikelihoodForEachSiteForEachState() const = 0;
    
		/**
		 * @brief Get the likelihood for each site and each model class.
		 *
		 * @return A two-dimension vector with all log likelihoods:
		 * <code>V[i][j] =</code> likelihood of site i and model class j.
		 */
		virtual VVdouble getLikelihoodForEachSiteForEachClass() const = 0;
	
    /**
		 * @brief Get the likelihood for each site and each model class and each state.
		 *
		 * @return A three-dimension vector with all log likelihoods:
		 * <code>V[i][j][k} =</code> likelihood of site i and model class j and state k.
		 */
		virtual VVVdouble getLikelihoodForEachSiteForEachClassForEachState() const = 0;
	  /** @} */

    /**
		 * @brief Get the posterior model class (the one with maximum posterior
		 * probability) for each site.
		 *
		 * @return A vector with all model classes indexes.
		 */
		virtual std::vector<size_t> getClassWithMaxPostProbOfEachSite() const = 0;

   
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
     * @brief Get the parameters associated to substitution model(s).
     *
     * @return A ParameterList.
     */
    virtual ParameterList getSubstitutionProcessParameters() const = 0;

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
     * @brief Get the index (used for inner computations) of a given site (original alignment column).
     *
     * @param site An alignment position.
     * @return The site index corresponding to the given input alignment position.
     */
    virtual size_t getSiteIndex(size_t site) const throw (IndexOutOfBoundsException) = 0;   
    /** @} */

    /**
     * @brief Tell if derivatives must be computed.
     *
     * This methods calls the enableFirstOrderDerivatives and enableSecondOrderDerivatives.
     *
     * @param yn Yes or no.
     */
    virtual void enableDerivatives(bool yn) = 0;

    /**
     * @name Iterators
     * @{
     */
    //TODO jdutheil on 21/04/13: need to account for model classes!
    //virtual ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const = 0;
  
    //jdutheil on 21/04/13: I think we will drop this type of iterator, which were never used before and are difficult to implement in the new framework...
    //virtual ConstSiteModelIterator* getNewSiteModelIterator(size_t siteIndex) const = 0;
    /* @} */


};

} //end of namespace newlik.
} //end of namespace bpp.

#endif  //_TREELIKELIHOOD_H_

