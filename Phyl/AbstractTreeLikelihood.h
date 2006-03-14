//
// File: AbstractTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Oct 17 17:57:21 2003
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

#ifndef _ABSTRACTTREELIKELIHOOD_H_
#define _ABSTRACTTREELIKELIHOOD_H_

#include "TreeLikelihood.h"
#include "Tree.h"
#include "TreeTemplate.h"

//From SeqLib:
#include <Seq/SiteContainer.h>

/**
 * @brief TreeLikelihood partial data structure.
 *
 * Stores inner computation for a given node.
 *
 * @see TreeLikelihoodData
 */
class TreeLikelihoodNodeData
{
	public:
		TreeLikelihoodNodeData() {}
		virtual ~TreeLikelihoodNodeData() {}

	public:
		/**
		 * @brief Get the node associated to this data structure.
		 *
		 * @return The node associated to this structure.
		 */
		virtual const Node * getNode() const = 0;
};

/**
 * @brief TreeLikelihood data structure.
 *
 * Stores all the inner computations:
 * - conditionnal likelihoods for each node,
 * - correspondance between sites in the dataset and array indices.
 * 
 * @see TreeLikelihoodNodeData
 */
class TreeLikelihoodData
{
	public:
		TreeLikelihoodData() {}
		virtual ~TreeLikelihoodData() {}

	public:
		virtual const Alphabet * getAlphabet() const = 0;
		virtual const TreeTemplate<Node> * getTree() const = 0;  
		virtual TreeTemplate<Node> * getTree() = 0;
		virtual unsigned int getArrayPosition(const Node* parent, const Node* son, unsigned int currentPosition) const = 0;
		virtual unsigned int getRootArrayPosition(const unsigned int site) const = 0;
		virtual       TreeLikelihoodNodeData & getNodeData(const Node * node) = 0;
		virtual const TreeLikelihoodNodeData & getNodeData(const Node * node) const = 0;
    
};


/**
 * @brief Partial implementation of the TreeLikelihoodData interface.
 *
 * This data structure provides a simple compression, by performing and storing computations
 * only one time per identical sites.
 *
 * The compression is achieved by the TreeLikelihood object.
 * The correspondance between sites in the dataset and the arrays in the structures is given
 * by the _rootPatternLinks array: the array indice for site @f$i@f$ if given by:
 * @code
 * _rootPatternLinks[i]
 * @endcode
 *
 * Finally, the _rootWeights array gives for each array position, the number of sites with this
 * pattern.
 * The global likelihood is then given by the product of all likelihoods for each array position,
 * weighted by the corresponding number of sites.
 */
class AbstractTreeLikelihoodData :
	public virtual TreeLikelihoodData
{
	public:
		AbstractTreeLikelihoodData() {}
		virtual ~AbstractTreeLikelihoodData() {}

	protected:
		/**
		 * @brief Links between sites and patterns.
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

		/**
		 * @brief The frequency of each site.
		 */
		vector<unsigned int> _rootWeights;

		TreeTemplate<Node> * _tree;

		const Alphabet * _alphabet;

	public:
		vector<unsigned int> & getRootArrayPositions() { return _rootPatternLinks; }
		const vector<unsigned int> & getRootArrayPositions() const { return _rootPatternLinks; }
		unsigned int getRootArrayPosition(const unsigned int site) const
		{
			return _rootPatternLinks[site];
		}
		unsigned int getWeight(unsigned int pos) const
		{
			return _rootWeights[pos];
		}
		const vector<unsigned int> & getWeights() const
		{ 
			return _rootWeights;
		}

		const Alphabet * getAlphabet() const { return _alphabet; }
		const TreeTemplate<Node> * getTree() const { return _tree; }  
		      TreeTemplate<Node> * getTree()       { return _tree; }
};

/**
 * @brief Partial implementatoin of the TreeLikelihood interface. 
 *
 * This class implements a few methods useful for most of likelihood
 * computation methods.
 * 
 * - The Parametrizable interface;
 * - The getTree() method;
 * 
 * It also adds an abstract method for recursive computations.
 */
class AbstractTreeLikelihood :
	public virtual TreeLikelihood,
	public virtual AbstractParametrizable
{
	protected:
		//const Alphabet * _alphabet;
		const SiteContainer * _data;
		mutable TreeTemplate<Node> *  _tree;
		bool _computeDerivatives;

	public:
		AbstractTreeLikelihood() {}
		virtual ~AbstractTreeLikelihood() {}
	
	public:
		/**
		 * @name The TreeLikelihood interface.
		 *
		 * @{
		 */
		const SiteContainer * getData() const { return _data; }
		const Alphabet * getAlphabet() const { return _data -> getAlphabet(); }	
		Vdouble getLikelihoodForEachSite()                 const;
		Vdouble getLogLikelihoodForEachSite()              const;
		VVdouble getLikelihoodForEachSiteForEachState()    const;
		VVdouble getLogLikelihoodForEachSiteForEachState() const;
		unsigned int getNumberOfSites() const { return _data -> getNumberOfSites(); }
		unsigned int getNumberOfStates() const { return _data -> getAlphabet() -> getSize(); }
		Tree * getTree() const { return _tree; }
		void setComputeDerivatives(bool yn) { _computeDerivatives = yn; }
		bool computeDerivatives() const { return _computeDerivatives; }
		/** @} */

	protected:
		
		/**
		 * @brief Recompute _pxy, _dpxy and _d2pxy arrays, and derivatives if needed.
		 *
		 * This method is called when some parameter has changed.
		 *
		 * @param params The parameters that changed.
		 */
		virtual void fireParameterChanged(const ParameterList & params) = 0;
		
};


#endif	//_ABSTRACTTREELIKELIHOOD_H_

