//
// File: HomogeneousTreeLikelihood.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Fri Oct 17 18:14:51 2003
//

#ifndef _HOMOGENEOUSTREELIKELIHOOD_H_
#define _HOMOGENEOUSTREELIKELIHOOD_H_

#include "AbstractHomogeneousTreeLikelihood.h"
#include "SubstitutionModel.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>
#include <NumCalc/DiscreteDistribution.h>

// From the STL:
#include <map>

using namespace std;

/**
 * @brief This class implement the 'traditional' way of computing likelihood for a tree.
 *
 * The substitution model is constant over the tree (homogeneous model).
 * A non uniform distribution of rates among the sites is allowed (ASRV models).</p>
 *
 * The Felsenstein recursive algorithm is used for conputation.
 * For each node, a likelihood tensor is defined, containing all likelihoods values for
 * for the substree defined by this node. The likelihood tensor dimension are defined as below:
 * <ul>
 * <li>Site</li>
 * <li>Rate class</li>
 * <li>Ancestral state</li>
 * </ul>
 * These tensors are stored into a map with each node as a key (cf. _likelihoods).
 *
 * The computation use <i>site patterns</i> for more efficiency.
 * Following N. Galtier (personal communication ;-), we define a Pattern as a distinct site
 * in a sub-dataset corresponding to the dataset with sequences associated to a particular subtree.
 * The likelihood computation is the same for a given site, hence the idea is to save time from
 * performing many times the same coputation.
 * The network between all patterns is defined by the _patternLinks double map, initialized in the
 * initLikelihoodsWithPatterns() method. This initialisation takes more time than the classic
 * initTreeLikelihood one, where all likelihoods for a given site <i>i</i> are at the <i>i</i> coordinate
 * in the likelihood tensor, but is really faster when computing the likelihood (computeLikelihoods() method).
 * Hence, if you have to compute likelihood many times while holding the tree topology unchanged,
 * you should use patterns. And since this is what you'll have to do in most case (for instance for parameter
 * estimation), we set this as the default method for now. We provide the second method for topology estimation
 * methods (far from achieved!)
 */
class HomogeneousTreeLikelihood : public AbstractHomogeneousTreeLikelihood
{
	protected:

		/**
		 * @brief This contains all likelihood values used for computation.
		 *
		 * <pre>
		 * x[n][i][c][s]
		 *   |------------> Node n (pointer)
		 *      |---------> Site i
		 *         |------> Rate class c
		 *            |---> Ancestral state s
		 * </pre> 
		 * We call this the <i>likelihood array</i> for each node.
		 */
		mutable map<const Node *, VVVdouble> _likelihoods;
	
		/**
		 * @brief This contains all likelihood first order derivatives values used for computation.
		 *
		 * <pre>
		 * x[n][i][c][s]
		 *   |------------> Node n (pointer)
		 *      |---------> Site i
		 *         |------> Rate class c
		 *            |---> Ancestral state s
		 * </pre> 
		 * We call this the <i>dLikelihood array</i> for each node.
		 */
		mutable map<const Node *, VVVdouble> _dLikelihoods;
	
		/**
		 * @brief This contains all likelihood second order derivatives values used for computation.
		 *
		 * <pre>
		 * x[n][i][c][s]
		 *   |------------> Node n (pointer)
		 *      |---------> Site i
		 *         |------> Rate class c
		 *            |---> Ancestral state s
		 * </pre> 
		 * We call this the <i>d2Likelihood array</i> for each node.
		 */
		mutable map<const Node *, VVVdouble> _d2Likelihoods;
		
		/**
		 * @brief This map defines the pattern network.
		 *
		 * Let n1 be a node in the tree, and n11 and n12 its sons.
		 * Providing the likelihood array is known for nodes n11 and n12,
		 * the likelihood array for node n1 and site <i>i</i> (_likelihood[n1][i]) must be computed	
		 * using arrays _patternLinks[n1][n11][i] and _patternLinks[n1][n12][i].
		 * This network is intialized once for all in the constructor of this class.
		 *
		 * The double map contains the position of the site to use (second dimension)
		 * of the likelihoods array.
		 * Pointer are no longer used, since the pattern network is used for both
		 * likelihoods arrays and dLikelihoods arrays.
		 */
		map< const Node *, map< const Node *, vector<unsigned int> > > _patternLinks;
		
		/**
		 * @brief As previous, but for the global container.
		 *
		 * The size of this vector is equal to the number of sites in the container,
		 * each element corresponds to a site in the container and points to the
		 * corresponding column in the likelihood array of the root node.
		 * If the container contains no repeated site, there will be a strict
		 * equivalence between each site and the likelihood array of the root node.
		 * However, if this is not the case, some pointers may point toward the same
		 * element in the likelihood array.
		 */
		vector<unsigned int> _rootPatternLinks;
		
	public:
		HomogeneousTreeLikelihood(
			Tree & tree,
			const SiteContainer & data,
			SubstitutionModel * model,
			DiscreteDistribution * rDist,
			bool verbose = true)
			throw (Exception);
	
		virtual ~HomogeneousTreeLikelihood();
	
	public:

		/**
		 * @name The TreeLikelihood interface.
		 *
		 * Other methods are implemented in the AbstractHomogeneousTreeLikelihood class.
		 *
		 * @{
		 */
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
		VVdouble getLikelihoodForEachSiteForEachRate() const;
		VVdouble getLogLikelihoodForEachSiteForEachRate() const;
		VVdouble getPosteriorProbabilitiesOfEachRate() const;
		Vdouble  getRateWithMaxPostProbOfEachSite() const;
		Vint     getRateClassWithMaxPostProbOfEachSite() const;
		Vdouble  getPosteriorRateOfEachSite() const;
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
		 * @brief Implements the DerivableFirstOrder interface.
		 *
		 * Update the parameter list and call the applyParameters() method.
		 * Then compute the likelihoods at each node (computeLikelihood() method),
		 * compute the first order derivatives and call the getDLogLikelihood() method.
		 *
		 * If a subset of the whole parameter list is passed to the function,
		 * only these parameters are updated and the other remain constant (i.e.
		 * equal to their last value).
		 *
		 * @param parameters The parameter list to pass to the function.
		 */
		double getFirstOrderDerivative(const string & variable) const throw (Exception);

		/**
		 * @brief Implements the DerivableSecondOrder interface.
		 *
		 * Update the parameter list and call the applyParameters() method.
		 * Then compute the likelihoods at each node (computeLikelihood() method),
		 * compute the first and second order derivatives and call
		 * the getD2LogLikelihood() method.
		 *
		 * If a subset of the whole parameter list is passed to the function,
		 * only these parameters are updated and the other remain constant (i.e.
		 * equal to their last value).
		 *
		 * @param parameters The parameter list to pass to the function.
		 */
		double getSecondOrderDerivative(const string & variable) const throw (Exception);
		double getSecondOrderDerivative(const string & variable1, const string & variable2) const throw (Exception) { return 0; } // Not implemented for now.
	
		
	public:	// Specific methods:
	
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
		 * @brief This method initializes the leaves according to a sequence file.
		 * likelihood is set to 1 for the state corresponding to the sequence site,
		 * otherwise it is set to 0.
		 *
		 * All likelihood arrays at each nodes are initialized according to alphabet
		 * size and sequences length, and filled with 1.
		 *
		 * NB: This method is recursive.
		 *
		 * @param node      The node defining the subtree to analyse.
		 * @param sequences The data to be used for initialization.
		 */
		virtual void initTreeLikelihoods(const Node * node, const SiteContainer & sequences) throw (Exception);

		/**
		 * @brief This method initializes the leaves according to a sequence file.
		 *
		 * likelihood is set to 1 for the state corresponding to the sequence site,
		 * otherwise it is set to 0.
		 *
		 * All likelihood arrays at each nodes are initialized according to alphabet
		 * size and sequences length, and filled with 1.
		 *
		 * NB: This method is recursive.
		 *
		 * @param node      The node defining the subtree to analyse.
		 * @param sequences The data to be used for initialization.
		 * @return The shrunk sub-dataset for the subtree defined by <i>node</i>.
		 */
		virtual SiteContainer * initTreeLikelihoodsWithPatterns(const Node * node, const SiteContainer & sequences) throw (Exception);
		
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


#endif	//_HOMOGENEOUSTREELIKELIHOOD_H_
