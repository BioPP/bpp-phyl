//
// File: TreeLikelihood.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Fri Oct 17 17:36:44 2003
//

#ifndef _TREELIKELIHOOD_H_
#define _TREELIKELIHOOD_H_

#include "Tree.h"

// From NumCalc:
#include <NumCalc/ParameterList.h>
#include <NumCalc/Parametrizable.h>
#include <NumCalc/Functions.h>
#include <NumCalc/VectorTools.h>

/**
 * @brief The TreeLikelihood interface.
 
 * This interface defines the methods needed for computing the likelihood
 * of a phylogenetic tree, given a dataset.
 *
 * Likelihood computing isvery often memory and CPU expensive,
 * hence many algorithms try to store as information as possible to save computations.
 * This interface separates the computation itself (computeLikelihood() method) and the
 * result (getLikelihood() methods).
 */ 
class TreeLikelihood: public DerivableSecondOrder
{
	public:
		virtual ~TreeLikelihood();
	
	public:
		
		/**
		 * @brief Compute likelihood.
		 */
		virtual void computeTreeLikelihood() = 0;

		/**
		 * @brief Get the likelihood for a site.
		 *
		 * @param site The site index to analyse.
		 * @return The likelihood for site <i>site</i>.
		 */
		virtual double getLikelihoodForASite(unsigned int site) const = 0;

		/**
		 * @brief Get the logarithm of the likelihood for a site.
		 *
		 * @param site The site index to analyse.
		 * @return The logarithm of the likelihood for site <i>site</i>.
		 */
		virtual double getLogLikelihoodForASite(unsigned int site) const = 0;

		/**
		 * @brief Get the likelihood for a site and for a state.
		 *
		 * @param site The site index to analyse.
		 * @param state The state to consider.
		 * @return The likelihood for site <i>site</i> and state <i>state</i>.
		 */
		virtual double getLikelihoodForASiteForAState(unsigned int site, int state) const = 0;

		/**
		 * @brief Get the logarithm of the likelihood for a site and for a state.
		 *
		 * @param site The site index to analyse.
		 * @param state The state to consider.
		 * @return The logarithm of the likelihood for site <i>site</i> and state <i>state</i>.
		 */
		virtual double getLogLikelihoodForASiteForAState(unsigned int site, int state) const = 0;

		/**
		 * @brief Get the likelihood for each site.
		 *
		 * @return A vector with all likelihoods for each site.
		 */
		virtual Vdouble getLikelihoodForEachSite() const = 0;

		/**
		 * @brief Get the logarithm of the likelihood for each site.
		 *
		 * @return A vector with all log likelihoods for each site.
		 */
		virtual Vdouble getLogLikelihoodForEachSite() const = 0;

		/**
		 * @brief Get the likelihood for each site and for each state.
		 *
		 * @return A 2d vector with all likelihoods for each site and for each state.
		 */
		virtual VVdouble getLikelihoodForEachSiteForEachState() const = 0;

		/**
		 * @brief Get the logarithm of the likelihood for each site and for each state.
		 *
		 * @return A 2d vector with all log likelihoods for each site and for each state.
		 */
		virtual VVdouble getLogLikelihoodForEachSiteForEachState() const = 0;
		
		/**
		 * @brief Get the likelihood for the whole dataset.
		 *
		 * @return The likelihood of the dataset.
		 */
		virtual double getLikelihood() const = 0;

		/**
		 * @brief Get the logarithm of the likelihood for the whole dataset.
		 *
		 * @return The logarithm of the likelihood of the dataset.
		 */
		virtual double getLogLikelihood() const = 0;
	
		/**
		 * @brief Get the tree (topology and branch lengths).
		 *
		 * @return The tree of this TreeLikelihood object.
	 	 */
		virtual Tree<Node> * getTree() const = 0;

		/**
		 * @brief Get the number of sites in the dataset.
		 *
		 * @return the number of sites in the dataset.
		 */
		virtual unsigned int getNumberOfSites() const = 0;

		/**
		 * @brief Get the number of states in the alphabet associated to the dataset.
		 *
		 * @return the number of states in the alphabet associated to the dataset.
		 */		
		virtual unsigned int getNumberOfStates() const = 0;
		
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
		 * @brief Get the parameters assoicated to substitution model(s).
		 *
		 * @return A ParameterList.
		 */
		virtual ParameterList getSubstitutionModelParameters() const = 0;
		
		/** @} */
};


#endif	//_TREELIKELIHOOD_H_
