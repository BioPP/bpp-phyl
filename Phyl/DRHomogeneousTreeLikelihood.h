//
// File: DRHomogeneousTreeLikelihood.h
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
 * @brief Likelihood data structure for a leaf.
 * 
 * This class is for use with the DRASDRTreeParsimonyData class.
 * 
 * Store the likelihoods arrays associated to a leaf.
 * 
 * @see DRASDRTreeLikelihoodData
 */
class DRASDRTreeLikelihoodLeafData :
	public TreeLikelihoodNodeData
{
	protected:
		mutable VVdouble _leafLikelihood;
		const Node * _leaf;

	public:
		const Node * getNode() const { return _leaf; }
		void setNode(const Node * node) { _leaf = node; }

		VVdouble & getLikelihoodArray()	{	return _leafLikelihood;	}
};

/**
 * @brief Likelihood data structure for a node.
 * 
 * This class is for use with the DRASDRTreeParsimonyData class.
 * 
 * Store for each neighbor node an array with conditionnal likelihoods.
 *
 * @see DRASDRTreeLikelihoodData
 */
class DRASDRTreeLikelihoodNodeData :
	public TreeLikelihoodNodeData
{
	protected:
		/**
		 * @brief This contains all likelihood values used for computation.
		 *
		 * <pre>
		 * x[b][i][c][s]
		 *   |------------> Neighbor node of n (pointer)
 		 *      |---------> Site i
		 *         |------> Rate class c
		 *            |---> Ancestral state s
		 * </pre>
		 * We call this the <i>likelihood array</i> for each node.
		 */

		mutable map<const Node *, VVVdouble > _nodeLikelihoods;
		/**
		 * @brief This contains all likelihood first order derivatives values used for computation.
		 *
		 * <pre>
		 * x[i]
		 *   |---------> Site i
		 * </pre> 
		 * We call this the <i>dLikelihood array</i> for each node.
		 */
		mutable Vdouble _nodeDLikelihoods;
	
		/**
		 * @brief This contains all likelihood second order derivatives values used for computation.
		 *
		 * <pre>
		 * x[i]
		     |---------> Site i
		 * </pre> 
		 * We call this the <i>d2Likelihood array</i> for each node.
		 */
		mutable Vdouble _nodeD2Likelihoods;
		
		const Node * _node;

	public:
		DRASDRTreeLikelihoodNodeData() {}
		~DRASDRTreeLikelihoodNodeData() {}

	public:
		const Node * getNode() const { return _node; }
		void setNode(const Node * node) { _node = node; }

		const map<const Node *, VVVdouble> & getLikelihoodArrays() const { return _nodeLikelihoods; }
		map<const Node *, VVVdouble> & getLikelihoodArrays() { return _nodeLikelihoods; }
		VVVdouble & getLikelihoodArrayForNeighbor(const Node * neighbor)
		{
			return _nodeLikelihoods[neighbor];
		}
		const VVVdouble & getLikelihoodArrayForNeighbor(const Node * neighbor) const
		{
			return _nodeLikelihoods[neighbor];
		}
		Vdouble & getDLikelihoodArray() { return _nodeDLikelihoods;	}
		const Vdouble & getDLikelihoodArray() const	{	return _nodeDLikelihoods;	}
		Vdouble & getD2LikelihoodArray()	{	return _nodeD2Likelihoods; }
		const Vdouble & getD2LikelihoodArrayForNeighbor() const	{ return _nodeD2Likelihoods; }

		bool isNeighbor(const Node * neighbor) const
		{
			return _nodeLikelihoods.find(neighbor) != _nodeLikelihoods.end();
		}

		void eraseNeighborArrays()
		{
			_nodeLikelihoods.erase(_nodeLikelihoods.begin(), _nodeLikelihoods.end());
			_nodeDLikelihoods.erase(_nodeDLikelihoods.begin(), _nodeDLikelihoods.end());
			_nodeD2Likelihoods.erase(_nodeD2Likelihoods.begin(), _nodeD2Likelihoods.end());
		}
};

/**
 * @brief Likelihood data structure for rate across sites models, using a double-recursive algorithm.
 */
class DRASDRTreeLikelihoodData :
	public virtual AbstractTreeLikelihoodData
{
	protected:

		mutable map<const Node *, DRASDRTreeLikelihoodNodeData> _nodeData;
		mutable map<const Node *, DRASDRTreeLikelihoodLeafData> _leafData;
		mutable VVVdouble _rootLikelihoods;
		mutable VVdouble  _rootLikelihoodsS;
		mutable Vdouble   _rootLikelihoodsSR;

		SiteContainer * _shrunkData;
		unsigned int _nbSites; 
		unsigned int _nbStates;
		unsigned int _nbClasses;
		unsigned int _nbDistinctSites; 

	public:
		DRASDRTreeLikelihoodData(TreeTemplate<Node> & tree, unsigned int nbClasses) : _nbClasses(nbClasses) { _tree = &tree; }
		virtual ~DRASDRTreeLikelihoodData() { delete _shrunkData; }

	public:
		DRASDRTreeLikelihoodNodeData & getNodeData(const Node * node)
		{ 
			return _nodeData[node];
		}
		const DRASDRTreeLikelihoodNodeData & getNodeData(const Node * node) const
		{ 
			return _nodeData[node];
		}
		DRASDRTreeLikelihoodLeafData & getLeafData(const Node * node)
		{ 
			return _leafData[node];
		}
		const DRASDRTreeLikelihoodLeafData & getLeafData(const Node * node) const
		{ 
			return _leafData[node];
		}
		unsigned int getArrayPosition(const Node* parent, const Node* son, unsigned int currentPosition) const
		{
			return currentPosition;
		}

		const map<const Node *, VVVdouble> & getLikelihoodArrays(const Node *node) const 
		{
			return _nodeData[node].getLikelihoodArrays();
		}
		map<const Node *, VVVdouble> & getLikelihoodArrays(const Node *node)
		{
			return _nodeData[node].getLikelihoodArrays();
		}

		VVVdouble & getLikelihoodArray(const Node *parent, const Node *neighbor)
		{
			return _nodeData[parent].getLikelihoodArrayForNeighbor(neighbor);
		}
		const VVVdouble & getLikelihoodArray(const Node *parent, const Node *neighbor) const
		{
			return _nodeData[parent].getLikelihoodArrayForNeighbor(neighbor);
		}
		
		Vdouble & getDLikelihoodArray(const Node * node)
		{
			return _nodeData[node].getDLikelihoodArray();
		}
		
		const Vdouble & getDLikelihoodArray(const Node * node) const
		{
			return _nodeData[node].getDLikelihoodArray();
		}
		
		Vdouble & getD2LikelihoodArray(const Node * node)
		{
			return _nodeData[node].getD2LikelihoodArray();
		}

		const Vdouble & getD2LikelihoodArray(const Node * node) const
		{
			return _nodeData[node].getD2LikelihoodArray();
		}

		VVdouble & getLeafLikelihoods(const Node * node)
		{
			return _leafData[node].getLikelihoodArray();
		}
		const VVdouble & getLeafLikelihoods(const Node * node) const
		{
			return _leafData[node].getLikelihoodArray();
		}
		VVVdouble & getRootLikelihoodArray() { return _rootLikelihoods; }
		VVdouble  & getRootSiteLikelihoodArray() { return _rootLikelihoodsS; }
		Vdouble   & getRootRateSiteLikelihoodArray() { return _rootLikelihoodsSR; }

		unsigned int getNumberOfDistinctSites() const { return _nbDistinctSites; }
		unsigned int getNumberOfSites() const { return _nbSites; }
		unsigned int getNumberOfStates() const { return _nbStates; }
		unsigned int getNumberOfClasses() const { return _nbClasses; }

		const SiteContainer * getShrunkData() const { return _shrunkData; }
		
		void init(const SiteContainer & sites, const SubstitutionModel & model) throw (Exception);
		//void reInit() throw (Exception);

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
		 * @param model     The model, used for initializing leaves' likelihoods.
		 */
		void initTreeLikelihoods(const Node * node, const SequenceContainer & sequences, const SubstitutionModel & model) throw (Exception);
		//void reInit(const Node * node) throw (Exception);

};



/**
 * @brief This class implements the likelihood computation for a tree using the double-recursive
 * algorithm.
 *
 * The substitution model is the same over the tree (homogeneous model).
 * A non-uniform distribution of rates among the sites is allowed (ASRV models).</p>
 *
 * The Felsenstein recursive algorithm is used for conputation.
 * As in the HomogeneousTreeLikelihood class, a likelihood tensor is defined:
 * 
 * -Site
 * -Rate class
 * -Ancestral state
 * 
 * However, in this class, a node will be attached a set of tensor instead of one single tensor,
 * one tensor for each subtree it defines.
 * These tensors are stored into a map of map, with each node as a primary key and each neighbor
 * of the node (sons + father) a a secondary key.
 *
 * All nodes share the same site patterns.
 */
class DRHomogeneousTreeLikelihood : public virtual AbstractHomogeneousTreeLikelihood
{
	protected:
		mutable DRASDRTreeLikelihoodData *_likelihoodData;
		
		//some values we'll need:
		unsigned int _nbDistinctSites; //the number of distinct sites in the container
		
	public:
		DRHomogeneousTreeLikelihood(
			TreeTemplate<Node> & tree,
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

		DRASDRTreeLikelihoodData * getLikelihoodData() { return _likelihoodData; }
		const DRASDRTreeLikelihoodData * getLikelihoodData() const { return _likelihoodData; }
	
		virtual VVVdouble computeLikelihoodAtNode(const Node * node) const;

		/**
		 * @brief Retrieves all Pij(t) for a particular node.
		 *
		 * These intermediate results may be used by other methods.
		 */
		virtual const VVVdouble & getTransitionProbabilitiesForNode(const Node * node) const { return _pxy[node]; }
		
	protected:
	
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

		virtual void computeTreeDLikelihoodAtNode(const Node * node);
		virtual void computeTreeDLikelihoods();
		
		virtual void computeTreeD2LikelihoodAtNode(const Node * node);
		virtual void computeTreeD2Likelihoods();

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

