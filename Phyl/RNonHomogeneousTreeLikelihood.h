//
// File: RNonHomogeneousTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Tue Oct 9 16:03 2007
// From file: RHomogeneousTreeLikelihood.h
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _RNONHOMOGENEOUSTREELIKELIHOOD_H_
#define _RNONHOMOGENEOUSTREELIKELIHOOD_H_

#include "AbstractNonHomogeneousTreeLikelihood.h"
#include "SubstitutionModelSet.h"
#include "DRASRTreeLikelihoodData.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>
#include <NumCalc/DiscreteDistribution.h>

namespace bpp
{

/**
 * @brief This class implement the 'traditional' way of computing likelihood for a tree, allowing for non-homogeneous models of substitutions.
 *
 * A non uniform distribution of rates among the sites is allowed (ASRV models).</p>
 *
 * This class uses an instance of the DRASRTreeLikelihoodData for conditionnal likelihood storage.
 *
 * This class can also use a simple or recursive site compression.
 * In the simple case, computations for identical sites are not duplicated.
 * In the recursive case, computations for identical sub-sites (<i>site patterns </i>) are also not duplicated:
 * Following N. Galtier (personal communication), we define a Pattern as a distinct site
 * in a sub-dataset corresponding to the dataset with sequences associated to a particular subtree.
 * The likelihood computation is the same for a given site, hence the idea is to save time from
 * performing many times the same coputation.
 * The network between all patterns is defined by the _patternLinks double map, initialized in the
 * initLikelihoodsWithPatterns() method. This initialisation takes more time than the classic
 * initTreeLikelihood one, where all likelihoods for a given site <i>i</i> are at the <i>i</i> coordinate
 * in the likelihood tensor, but is really faster when computing the likelihood (computeLikelihoods() method).
 * Hence, if you have to compute likelihood many times while holding the tree topology unchanged,
 * you should use patterns.
 * This decreases the likelihood computation time, but at a cost: some time is spent to establish the patterns
 * relationships. Whether to use or not patterns depends on what you actllay need:
 * - The more you compute likelihoods without changing the data or topology, the more patterns are interesting
 *   (this divides the cost of computing patterns by the number of computation performed).
 *   Patterns are hence usefull when you have a high number of computation to perform, while optimizing numerical
 *   parameters for instance).
 * - Patterns are more likely to occur whith small alphabet (nucleotides).
 *
 * Important note: The input tree will be considered as rooted, since the likelihood of non-stationary models
 * depends on the position of the root. If the input tree is not rooted, it will be considered as a rotted tree
 * with a root multifurcation.
 */
class RNonHomogeneousTreeLikelihood :
	public AbstractNonHomogeneousTreeLikelihood
{
	protected:

		mutable DRASRTreeLikelihoodData *_likelihoodData;

	public:
    /**
     * @brief Build a new NonHomogeneousTreeLikelihood object without data.
     *
     * This constructor only initialize the parameters.
     * To compute a likelihood, you will need to call the setData() and the computeTreeLikelihood() methods.
     *
     * @param tree The tree to use.
     * @param modelSet The set of substitution models to use.
     * @param rDist The rate across sites distribution to use.
     * If true, any rooted tree will be unrooted before likelihood computation.
     * @param verbose Should I display some info?
     * @param usePatterns Tell if recursive site compression should be performed.
     * @throw Exception in an error occured.
     */
		RNonHomogeneousTreeLikelihood(
			const Tree & tree,
			SubstitutionModelSet * modelSet,
			DiscreteDistribution * rDist,
			bool verbose = true,
      bool usePatterns = true)
			throw (Exception);
	
    /**
     * @brief Build a new NonHomogeneousTreeLikelihood object and compute the corresponding likelihood.
     *
     * This constructor initializes all parameters, data, and likelihood arrays.
     *
     * @param tree The tree to use.
     * @param data Sequences to use.
     * @param modelSet The set of substitution models to use.
     * @param rDist The rate across sites distribution to use.
     * @param verbose Should I display some info?
     * @param usePatterns Tell if recursive site compression should be performed.
     * @throw Exception in an error occured.
     */
		RNonHomogeneousTreeLikelihood(
			const Tree & tree,
			const SiteContainer & data,
			SubstitutionModelSet * modelSet,
			DiscreteDistribution * rDist,
			bool verbose = true,
      bool usePatterns = true)
			throw (Exception);

    RNonHomogeneousTreeLikelihood(const RNonHomogeneousTreeLikelihood & lik);
    
    RNonHomogeneousTreeLikelihood & operator=(const RNonHomogeneousTreeLikelihood & lik);

    virtual ~RNonHomogeneousTreeLikelihood();

#ifndef NO_VIRTUAL_COV
    RNonHomogeneousTreeLikelihood*
#else
    Clonable*
#endif
    clone() const { return new RNonHomogeneousTreeLikelihood(*this); }
	
	private:

    /**
     * @brief Method called by constructors.
     */
    void _init(bool usePatterns) throw (Exception);
	
  public:

		/**
		 * @name The TreeLikelihood interface.
		 *
		 * Other methods are implemented in the AbstractHomogeneousTreeLikelihood class.
		 *
		 * @{
		 */
    void setData(const SiteContainer & sites) throw (Exception);
  	double getLikelihood () const;
		double getLogLikelihood() const;
		double getLikelihoodForASite (unsigned int site) const;
		double getLogLikelihoodForASite(unsigned int site) const;
 		/** @} */

		
		/**
		 * @name The DiscreteRatesAcrossSites interface implementation:
		 *
		 * @{
		 */
		double getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
		double getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
		double getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const;
		double getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const;
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
		void setParameters(const ParameterList & parameters) throw (ParameterNotFoundException, ConstraintException);
		double getValue() const throw(Exception);
		
		/**
		 * @name DerivableFirstOrder interface.
		 *
		 * @{
		 */
		double getFirstOrderDerivative(const string & variable) const throw (Exception);
		/** @} */

		/**
		 * @name DerivableSecondOrder interface.
		 *
		 * @{
		 */
		double getSecondOrderDerivative(const string & variable) const throw (Exception);
		double getSecondOrderDerivative(const string & variable1, const string & variable2) const throw (Exception) { return 0; } // Not implemented for now.
		/** @} */
	
	public:	// Specific methods:
	
    DRASRTreeLikelihoodData * getLikelihoodData() { return _likelihoodData; }
    const DRASRTreeLikelihoodData * getLikelihoodData() const { return _likelihoodData; }
 
    void computeTreeLikelihood();

		virtual double getDLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;

		virtual double getDLikelihoodForASite(unsigned int site) const;

		virtual double getDLogLikelihoodForASite(unsigned int site) const;
		
		virtual double getDLogLikelihood() const;
		
		virtual void computeTreeDLikelihood(const string & variable);

		virtual double getD2LikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;

		virtual double getD2LikelihoodForASite(unsigned int site) const;

		virtual double getD2LogLikelihoodForASite(unsigned int site) const;
		
		virtual double getD2LogLikelihood() const;
		
		virtual void computeTreeD2Likelihood(const string & variable);

	
	protected:
			
		/**
		 * @brief Compute the likelihood for a subtree defined by the Tree::Node <i>node</i>.
		 *
		 * @param node The root of the subtree.
		 */
		virtual void computeSubtreeLikelihood(const Node * node); //Recursive method.			

		virtual void computeDownSubtreeDLikelihood(const Node *);
		
		virtual void computeDownSubtreeD2Likelihood(const Node *);
	
		void fireParameterChanged(const ParameterList & params);
	
		/**
		 * @brief This method is mainly for debugging purpose.
		 *
		 * @param node The node at which likelihood values must be displayed.
		 */
		virtual void displayLikelihood(const Node * node);

};


} //end of namespace bpp.

#endif	//_RNONHOMOGENEOUSTREELIKELIHOOD_H_

