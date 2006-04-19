//
// File: HomogeneousTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Oct 17 18:14:51 2003
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
 * @brief Likelihood data structure for a node.
 * 
 * This class is for use with the DRASRTreeParsimonyData class.
 * 
 * Store all conditionnal likelihoods:
 * <pre>
 * x[i][c][s]
 *   |---------> Site i
 *      |------> Rate class c
 *         |---> Ancestral state s
 * </pre> 
 * We call this the <i>likelihood array</i> for each node.
 * In the same way, we store first and second order derivatives.
 *
 * @see DRASRTreeLikelihoodData
 */
class DRASRTreeLikelihoodNodeData :
	public TreeLikelihoodNodeData
{
	protected:
		mutable VVVdouble _nodeLikelihoods;
		mutable VVVdouble _nodeDLikelihoods;
		mutable VVVdouble _nodeD2Likelihoods;
		const Node * _node;

	public:
		const Node * getNode() const { return _node; }
		void setNode(const Node * node) { _node = node; }

		VVVdouble & getLikelihoodArray() { return _nodeLikelihoods; }
		const VVVdouble & getLikelihoodArray() const { return _nodeLikelihoods; }
		
		VVVdouble & getDLikelihoodArray() { return _nodeDLikelihoods; }
		const VVVdouble & getDLikelihoodArray() const { return _nodeDLikelihoods; }

		VVVdouble & getD2LikelihoodArray() { return _nodeD2Likelihoods; }
		const VVVdouble & getD2LikelihoodArray() const { return _nodeD2Likelihoods; }
};

/**
 * @brief discrete Rate Across Sites, (simple) Recursive likelihood data structure.
 */
class DRASRTreeLikelihoodData :
	public virtual AbstractTreeLikelihoodData
{
	protected:

		/**
		 * @brief This contains all likelihood values used for computation.
		 *
		 */
		mutable map<const Node *, DRASRTreeLikelihoodNodeData> _nodeData;
			
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
		 */
		mutable map< const Node *, map< const Node *, vector<unsigned int> > > _patternLinks;
		SiteContainer * _shrunkData;
		unsigned int _nbSites; 
		unsigned int _nbStates;
		unsigned int _nbClasses;
		unsigned int _nbDistinctSites; 

	public:
		DRASRTreeLikelihoodData(TreeTemplate<Node> & tree, unsigned int nbClasses) : _nbClasses(nbClasses)
    {
      _tree = &tree;
      _shrunkData = NULL;
    }
		virtual ~DRASRTreeLikelihoodData() { delete _shrunkData; }

	public:
		DRASRTreeLikelihoodNodeData & getNodeData(const Node * node)
		{ 
			return _nodeData[node];
		}
		const DRASRTreeLikelihoodNodeData & getNodeData(const Node * node) const
		{ 
			return _nodeData[node];
		}
		unsigned int getArrayPosition(const Node* parent, const Node* son, unsigned int currentPosition) const
		{
			return _patternLinks[parent][son][currentPosition];
		}
		unsigned int getRootArrayPosition(unsigned int currentPosition) const
		{
			return _rootPatternLinks[currentPosition];
		}
		const vector<unsigned int> & getArrayPositions(const Node* parent, const Node* son) const
		{
			return _patternLinks[parent][son];
		}
		vector<unsigned int> & getArrayPositions(const Node* parent, const Node* son)
		{
			return _patternLinks[parent][son];
		}
		unsigned int getArrayPosition(const Node* parent, const Node* son, unsigned int currentPosition)
		{
			return _patternLinks[parent][son][currentPosition];
		}


		VVVdouble & getLikelihoodArray(const Node *node)
		{
			return _nodeData[node].getLikelihoodArray();
		}
		
		VVVdouble & getDLikelihoodArray(const Node *node)
		{
			return _nodeData[node].getDLikelihoodArray();
		}
		
		VVVdouble & getD2LikelihoodArray(const Node *node)
		{
			return _nodeData[node].getD2LikelihoodArray();
		}

		unsigned int getNumberOfDistinctSites() const { return _nbDistinctSites; }
		unsigned int getNumberOfSites() const { return _nbSites; }
		unsigned int getNumberOfStates() const { return _nbStates; }
		unsigned int getNumberOfClasses() const { return _nbClasses; }
		
		void initLikelihoods(const SiteContainer & sites, const SubstitutionModel & model) throw (Exception);

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
     * @param model     The model to use.
		 */
		virtual void initLikelihoods(const Node * node, const SiteContainer & sequences, const SubstitutionModel & model) throw (Exception);

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
     * @param model     The model to use.
		 * @return The shrunk sub-dataset for the subtree defined by <i>node</i>.
		 */
		virtual SiteContainer * initLikelihoodsWithPatterns(const Node * node, const SiteContainer & sequences, const SubstitutionModel & model) throw (Exception);
		
};

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
 * estimation), we set this as the default method for now.
 * The second method is for testing purpose only.
 *
 * For topology estimation, consider using the DRHomogeneousTreeLikelihood class.
 */
class HomogeneousTreeLikelihood :
	public virtual AbstractHomogeneousTreeLikelihood
{
	protected:

		mutable DRASRTreeLikelihoodData *_likelihoodData;

	public:
    /**
     * @brief Build a new HomogeneousTreeLikelihood object.
     *
     * @param tree The tree to use.
     * @param model The substitution model to use.
     * @param rDist The rate across sites distribution to use.
     * @param checkRooted Tell if we have to check for the tree to be unrooted.
     * If true, any rooted tree will be unrooted before likelihood computation.
     * @param verbose Should I display some info?
     * @throw Exception in an error occured.
     */
		HomogeneousTreeLikelihood(
			TreeTemplate<Node> * tree,
			SubstitutionModel * model,
			DiscreteDistribution * rDist,
      bool checkRooted = true,
			bool verbose = true)
			throw (Exception);
	
    /**
     * @brief Build a new HomogeneousTreeLikelihood object.
     *
     * @param tree The tree to use.
     * @param data Sequences to use.
     * @param model The substitution model to use.
     * @param rDist The rate across sites distribution to use.
     * @param checkRooted Tell if we have to check for the tree to be unrooted.
     * If true, any rooted tree will be unrooted before likelihood computation.
     * @param verbose Should I display some info?
     * @throw Exception in an error occured.
     */
		HomogeneousTreeLikelihood(
			TreeTemplate<Node> * tree,
			const SiteContainer & data,
			SubstitutionModel * model,
			DiscreteDistribution * rDist,
      bool checkRooted = true,
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
    void setData(const SiteContainer & sites) throw (Exception);
  	double getLikelihood () const;
		double getLogLikelihood() const;
		double getLikelihoodForASite (unsigned int site) const;
		double getLogLikelihoodForASite(unsigned int site) const;
		double getLikelihoodForASiteForAState (unsigned int site, int state) const;
		double getLogLikelihoodForASiteForAState(unsigned int site, int state) const;
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
		VVdouble getPosteriorProbabilitiesOfEachRate() const;
		Vdouble  getRateWithMaxPostProbOfEachSite() const;
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


#endif	//_HOMOGENEOUSTREELIKELIHOOD_H_

