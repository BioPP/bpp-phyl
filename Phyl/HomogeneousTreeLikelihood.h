//
// File: HomogeneousTreeLikelihood.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Fri Oct 17 18:14:51 2003
//

#ifndef _HOMOGENEOUSTREELIKELIHOOD_H_
#define _HOMOGENEOUSTREELIKELIHOOD_H_

#include "AbstractTreeLikelihood.h"
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
class HomogeneousTreeLikelihood : public AbstractTreeLikelihood
{
	protected:
		SubstitutionModel * _model;
		DiscreteDistribution * _rateDistribution;

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
		mutable map<Node *, VVVdouble> _likelihoods;
	
		/**
		 * @brief This contains all likelihood derivatives values used for computation.
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
		mutable map<Node *, VVVdouble> _dLikelihoods;
	
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
		//map< Node *, map< Node *, vector< VVdouble *> > > _patternLinks;
		map< Node *, map< Node *, vector<int> > > _patternLinks;
		
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
		//vector<VVdouble *> _rootPatternLinks;
		vector<int> _rootPatternLinks;
		
		/**
		 * @brief Pointer toward all nodes in the tree.
		 *
		 * This is used for parameters estimation only, not for likelihood computation.
		 * The order of the nodes in the vector if the order of the named branches and
		 * is initalized once for all in the constructor. It then depends of the tree
		 * topology. This may lead to some problems when we'll act on tree topology...
		 */
		 vector<Node *> _nodes;
			
		//some values we'll need:
		unsigned int _nbSites,   //the number of sites in the container
		             _nbClasses, //the number of rate classes
		             _nbStates,  //the number of states in the alphabet
		             _nbNodes;   //the number of nodes in the tree
		
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
		 * Other methods are implemented in the AbstractTreeLikelihood class.
		 *
		 * @{
		 */
		double getLikelihood () const;
		double getLogLikelihood() const;
		double getLikelihoodForASite (unsigned int site) const;
		double getLogLikelihoodForASite(unsigned int site) const;
		ParameterList getBranchLengthsParameters() const;
		ParameterList getSubstitutionModelParameters() const;
		/** @} */

		// Add this one:
		
		/**
		 * @brief Get the likelihood for a site knowing its rate class.
		 *
		 * @param site      The site index.
		 * @param rateClass The rate class index.
		 * @return The likelihood for the specified site and rate class.
		 */
		virtual double getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
		
		/**
		 * @brief Get the logarithm of the likelihood for a site knowing its rate class.
		 *
		 * @param site      The site index.
		 * @param rateClass The rate class index.
		 * @return The logarithm of the likelihood for the specified site and rate class.
		 */
		virtual double getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
	
		/**
		 * @brief Get the likelihood for each site and each rate class.
		 *
		 * @return A two-dimension vector with all likelihoods.
		 */
		virtual VVdouble getLikelihoodForEachSiteForEachRate() const;
		
		/**
		 * @brief Get the logarithm of the likelihood for each site and each rate class.
		 *
		 * @return A two-dimension vector with all log likelihoods:
		 * <code>V[i][j] =</code> likelihood of site i and rate class j.
		 */
		virtual VVdouble getLogLikelihoodForEachSiteForEachRate() const;
		
		/**
		 * @brief Get the posterior probability for each site of belonging to a
		 * particular rate class.
		 *
		 * @return A two-dimension vector with all posterior probabilities:
		 * <code>V[i][j] =</code> probablity for site i of belonging to rate class j.
		 */
		virtual VVdouble getPosteriorProbabilitiesOfEachRate() const;
		
		/**
		 * @brief Get the posterior rate class (the one with maximum posterior
		 * probability) for each site.
		 *
		 * @return A vector with all rate classes indexes.
		 */
		virtual Vint getPosteriorRateClassOfEachSite() const;
	
		/**
		 * @brief Get the posterior rate (the one with maximum posterior
		 * probability) for each site.
		 *
		 * @return A vector with all rates.
		 */
		virtual Vdouble getPosteriorRateOfEachSite() const;

		/**
		 * @brief Get the substitution model used for the computation.
		 *
		 * @return A const pointer toward the substitution model of this instance.
		 */
		virtual const SubstitutionModel * getSubstitutionModel() const;
		
		/**
		 * @brief Get the substitution model used for the computation.
		 *
		 * @return A pointer toward the substitution model of this instance.
		 */
		virtual SubstitutionModel * getSubstitutionModel();
	
		/**
		 * @brief Get the rate distribution used for the computation.
		 *
		 * @return A const pointer toward the rate distribution of this instance.
		 */
		virtual const DiscreteDistribution * getRateDistribution() const;
		
		/**
		 * @brief Get the rate distribution used for the computation.
		 *
		 * @return A pointer toward the rate distribution of this instance.
		 */
		virtual DiscreteDistribution * getRateDistribution();
			  
		//The Optimizable interface is implemented here:
		
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
		double f(const ParameterList & parameters) const;
		
	public:	// Specific methods:
	
		/**
		 * @brief This builds the <i>parameters</i> list from all paramtrizable objects,
		 * <i>i.e.</i> substitution model, rate distribution and tree.
		 */
		virtual void initParameters();

		/**
		 * @brief This removes a particular parameter from the list.
		 *
		 * This method may be used to not estimate a parameter after having
		 * fixed its value. The previous method reset all calls of thos one.
		 *
		 * @param name The name of the parameter to ignore.
		 */
		virtual void ignoreParameter(const string & name) throw (ParameterNotFoundException);

		virtual ParameterList getRateDistributionParameters() const;
	
		virtual double getDLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;

		virtual double getDLikelihoodForASite(unsigned int site) const;

		virtual double getDLogLikelihoodForASite(unsigned int site) const;
		
		virtual double getDLogLikelihood() const;
		
		virtual void computeTreeDLikelihood(const string & variable);

		virtual void computeDownSubtreeDLikelihood(Node *);
		
		virtual double df(const string & variable, const ParameterList & parameters) const;
	
	
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
		virtual void initTreeLikelihoods(Node * node, const SiteContainer & sequences) throw (Exception);

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
		virtual SiteContainer * initTreeLikelihoodsWithPatterns(Node * node, const SiteContainer & sequences) throw (Exception);
		
		void computeSubtreeLikelihood(Node * node); //Recursive method.

		/**
		 * @brief All parameters are stores in a parameter list.
		 *
		 * This function apply these parameters to the substitution model,
		 * to the rate distribution and to the branch lengths.
		 */
		virtual void applyParameters() throw (Exception);	
	
		/**
		 * @brief This method is mainly for debuggin purpose.
		 *
		 * @param node The node at which likelihood values must be displayed.
		 */
		virtual void displayLikelihood(Node * node);
};


#endif	//_HOMOGENEOUSTREELIKELIHOOD_H_
