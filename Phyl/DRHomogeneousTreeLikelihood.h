//
// File: DRHomogeneousTreeLikelihood.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Fri Oct 17 18:14:51 2003
//

#ifndef _DRHOMOGENEOUSTREELIKELIHOOD_H_
#define _DRHOMOGENEOUSTREELIKELIHOOD_H_

#include "AbstractHomogeneousTreeLikelihood.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>
#include <NumCalc/DiscreteDistribution.h>

// From the STL:
#include <map>

using namespace std;

/**
 * @brief This class implement the computation of likelihood for a tree using the double-recursive
 * method.
 *
 * The substitution model is constant over the tree (homogeneous model).
 * A non uniform distribution of rates among the sites is allowed (ASRV models).</p>
 *
 * The Felsenstein recursive algorithm is used for conputation.
 * As in the HomogeneousTreeLikelihood class, a likelihood tensor is defined:
 * <ul>
 * <li>Site</li>
 * <li>Rate class</li>
 * <li>Ancestral state</li>
 * </ul>
 * However in this class, a node will be attached a set of tensor instead of one single tensor,
 * one tensor for each subtree it defines.
 * These tensors are stored into a map of map, with each node as a primary key and each neighbor
 * of the node (sons + father) a a secondary key.
 *
 * No pattern are used here.
 */
class DRHomogeneousTreeLikelihood : public AbstractHomogeneousTreeLikelihood
{
	protected:
		SiteContainer * _shrunkData;

		/**
		 * @brief This contains all likelihood values used for computation.
		 *
		 * <pre>
		 * x[n][b][i][c][s]
		 *   |---------------> Node n (pointer)
		 *      |------------> Neighbor node of n (pointer)
		 *         |---------> Site i
		 *            |------> Rate class c
		 *               |---> Ancestral state s
		 * </pre> 
		 * We call this the <i>likelihood array</i> for each node.
		 */
		mutable map<const Node *, map<const Node *, VVVdouble> > _likelihoods;
	
		/**
		 * @brief This contains all likelihood first order derivatives values used for computation.
		 *
		 * <pre>
		 * x[b][i][c][s]
		 *   |------------> Neighbor node of n (pointer)
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
		 * x[b][i][c][s]
		 *   |------------> Neighbor node of n (pointer)
		 *      |---------> Site i
		 *         |------> Rate class c
		 *            |---> Ancestral state s
		 * </pre> 
		 * We call this the <i>d2Likelihood array</i> for each node.
		 */
		mutable map<const Node *, VVVdouble> _d2Likelihoods;

		mutable map<const Node *, VVdouble> _leavesLikelihoods;

		mutable VVVdouble _rootLikelihoods;
						
		//some values we'll need:
		unsigned int _nbDistinctSites; //the numer of distinct sites in the container
		
	public:
		DRHomogeneousTreeLikelihood(
			Tree & tree,
			const SiteContainer & data,
			SubstitutionModel * model,
			DiscreteDistribution * rDist,
			bool verbose = true)
			throw (Exception);
	
		virtual ~DRHomogeneousTreeLikelihood();
	
	public:

		/**
		 * @name The TreeLikelihood interface.
		 *
		 * Other methods are implemented in the AbstractTreeLikelihood class.
		 *
		 * @{
		 */
		double getLikelihood () const;
		double getLogLikelihood() const;
		double getLikelihoodForASite (unsigned int site) const;
		double getLogLikelihoodForASite(unsigned int site) const;
		double getLikelihoodForASiteForAState (unsigned int site, int state) const;
		double getLogLikelihoodForASiteForAState(unsigned int site, int state) const;
		/** @} */

		void computeTreeLikelihood();

		
		/**
		 * @name The DiscreteRatesAcrossSites interface implementation:
		 *
		 * @{
		 */
		double getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
		double getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
		double getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const;
		double getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const;
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
		 * @name DerivableFirstOrder interface.
		 *
		 * @{
		 */
		double getFirstOrderDerivative(const string & variable) const throw (Exception);
		/** @{ */

		/**
		 * @name DerivableSecondOrder interface.
		 *
		 * @{
		 */
		double getSecondOrderDerivative(const string & variable) const throw (Exception);
		double getSecondOrderDerivative(const string & variable1, const string & variable2) const throw (Exception) { return 0; } // Not implemented for now.
		/** @} */
		
	public:	// Specific methods:
	
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
	
		virtual map<const Node *, VVVdouble> getLikelihoodArraysForNode(const Node * node)
	  { return _likelihoods[node]; }

		virtual VVVdouble getTransitionProbabilitiesForNode(const Node * node)
		{ return _pxy[node]; }

		virtual vector<unsigned int> getRootPatternLinks() const
		{ return _rootPatternLinks; }

		virtual VVVdouble getRootLikelihoods() const
		{ return _rootLikelihoods; }

		virtual unsigned int getNumberOfDistinctSites() const
		{ return _nbDistinctSites; }

		virtual VVVdouble computeLikelihoodAtNode(const Node * node);

		virtual VVdouble getLeafLikelihoods(const Node * node) const
		{ return _leavesLikelihoods[node]; }

		virtual const SiteContainer * getShrunkData() const
		{ return _shrunkData; }
		

	protected:
		
		/**
		 * @brief This method initializes the leaves according to a sequence container.
		 *
		 * Here the container _shrunkData is used.
		 * Likelihood is set to 1 for the state corresponding to the sequence site,
		 * otherwise it is set to 0.
		 *
		 * All likelihood arrays at each nodes are initialized according to alphabet
		 * size and sequences length, and filled with 1.
		 *
		 * NB: This method is recursive.
		 *
		 * @param node      The node defining the subtree to analyse.
		 * @param sequences The sequence container to use.
		 */
		virtual void initTreeLikelihoods(const Node * node, const SequenceContainer & sequences) throw (Exception);

		/**
		 * Initialize the arrays corresponding to each son node for the node passed as argument.
		 * The method is called for each son node and the result stored in the corresponding array.
		 */
		virtual void computeSubtreeLikelihoodPostfix(const Node * node); //Recursive method.
		/**
		 * This method initilize the remaining lieklihood arrays, corresponding to father nodes.
		 * It must be called after the postfix method because it requires that the arrays for
		 * son nodes to be be computed.
		 */
		virtual void computeSubtreeLikelihoodPrefix(const Node * node); //Recursive method.

		virtual void computeRootLikelihood();

		virtual void computeDownSubtreeDLikelihood(const Node *);

		virtual void computeDownSubtreeD2Likelihood(const Node *);
	
		void fireParameterChanged(const ParameterList & params);

		void resetLikelihoodArrays(const Node * node);
	
		/**
		 * @brief This method is mainly for debugging purpose.
		 *
		 * @param node The node at which likelihood values must be displayed.
		 */
		virtual void displayLikelihood(const Node * node);

};


#endif	//_DRHOMOGENEOUSTREELIKELIHOOD_H_
