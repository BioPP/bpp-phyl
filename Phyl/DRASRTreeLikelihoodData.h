//
// File: DRASRTreeLikelihoodData.h
// Created by: Julien Dutheil
// Created on: Sat Dec 30 14:20 2006
// From file HomogeneousTreeLikelihood.h
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

#ifndef _DRASRHOMOGENEOUSTREELIKELIHOODDATA_H_
#define _DRASRHOMOGENEOUSTREELIKELIHOODDATA_H_

#include "AbstractTreeLikelihoodData.h"
#include "SubstitutionModel.h"
#include "SitePatterns.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>

// From the STL:
#include <map>

using namespace std;

namespace bpp
{

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
    DRASRTreeLikelihoodNodeData * clone() const { return new DRASRTreeLikelihoodNodeData(*this); }

	public:
		const Node * getNode() const { return _node; }
		void setNode(const Node & node) { _node = &node; }

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
	public AbstractTreeLikelihoodData
{
	protected:

		/**
		 * @brief This contains all likelihood values used for computation.
		 *
		 */
		mutable map<int, DRASRTreeLikelihoodNodeData> _nodeData;
			
		/**
		 * @brief This map defines the pattern network.
		 *
		 * Let n1 be the id of a node in the tree, and n11 and n12 the ids of its sons.
		 * Providing the likelihood array is known for nodes n11 and n12,
		 * the likelihood array for node n1 and site <i>i</i> (_likelihood[n1][i]) must be computed	
		 * using arrays _patternLinks[n1][n11][i] and _patternLinks[n1][n12][i].
		 * This network is intialized once for all in the constructor of this class.
		 *
		 * The double map contains the position of the site to use (second dimension)
		 * of the likelihoods array.
		 */
		mutable map<int, map<int, vector<unsigned int> > > _patternLinks;
		SiteContainer * _shrunkData;
		unsigned int _nbSites; 
		unsigned int _nbStates;
		unsigned int _nbClasses;
		unsigned int _nbDistinctSites; 
    bool _usePatterns;

	public:
		DRASRTreeLikelihoodData(TreeTemplate<Node> & tree, unsigned int nbClasses, bool usePatterns = true):
      _nodeData(), _patternLinks(), _shrunkData(NULL), _nbSites(0), _nbStates(0),
      _nbClasses(nbClasses), _nbDistinctSites(0), _usePatterns(usePatterns)
    {
      _tree = &tree;
    }

		DRASRTreeLikelihoodData(const DRASRTreeLikelihoodData & data):
      AbstractTreeLikelihoodData(data),
      _nodeData(data._nodeData),
      _patternLinks(data._patternLinks),
      _shrunkData(NULL),
      _nbSites(data._nbSites), _nbStates(data._nbStates),
      _nbClasses(data._nbClasses), _nbDistinctSites(data._nbDistinctSites),
      _usePatterns(data._usePatterns)
    {
      _tree              = data._tree;
      if(data._shrunkData)
        _shrunkData      = dynamic_cast<SiteContainer *>(data._shrunkData->clone());
    }

		DRASRTreeLikelihoodData& operator=(const DRASRTreeLikelihoodData & data)
    {
      AbstractTreeLikelihoodData::operator=(data);
      _nodeData          = data._nodeData;
      _patternLinks      = data._patternLinks;
      _nbSites           = data._nbSites;
      _nbStates          = data._nbStates;
      _nbClasses         = data._nbClasses;
      _nbDistinctSites   = data._nbDistinctSites;
      _tree              = data._tree;
      if(data._shrunkData)
        _shrunkData      = dynamic_cast<SiteContainer *>(data._shrunkData->clone());
      else
        _shrunkData      = NULL;
      _usePatterns       = data._usePatterns;
      return *this;
    }

    virtual ~DRASRTreeLikelihoodData() { delete _shrunkData; }

    DRASRTreeLikelihoodData * clone() const { return new DRASRTreeLikelihoodData(*this); }

	public:
    void setTree(TreeTemplate<Node> & tree)
    { 
      _tree = &tree;
      for(map<int, DRASRTreeLikelihoodNodeData>::iterator it = _nodeData.begin(); it != _nodeData.end(); it++)
      {
        int id = it->second.getNode()->getId();
        it->second.setNode(*_tree->getNode(id));
      }
    }

		DRASRTreeLikelihoodNodeData & getNodeData(const Node * node)
		{ 
			return _nodeData[node->getId()];
		}
		const DRASRTreeLikelihoodNodeData & getNodeData(const Node * node) const
		{ 
			return _nodeData[node->getId()];
		}
		unsigned int getArrayPosition(const Node* parent, const Node* son, unsigned int currentPosition) const
		{
			return _patternLinks[parent->getId()][son->getId()][currentPosition];
		}
		unsigned int getRootArrayPosition(unsigned int currentPosition) const
		{
			return _rootPatternLinks[currentPosition];
		}
		const vector<unsigned int> & getArrayPositions(const Node* parent, const Node* son) const
		{
			return _patternLinks[parent->getId()][son->getId()];
		}
		vector<unsigned int> & getArrayPositions(const Node* parent, const Node* son)
		{
			return _patternLinks[parent->getId()][son->getId()];
		}
		unsigned int getArrayPosition(const Node* parent, const Node* son, unsigned int currentPosition)
		{
			return _patternLinks[parent->getId()][son->getId()][currentPosition];
		}

		VVVdouble & getLikelihoodArray(const Node *node)
		{
			return _nodeData[node->getId()].getLikelihoodArray();
		}
		
		VVVdouble & getDLikelihoodArray(const Node *node)
		{
			return _nodeData[node->getId()].getDLikelihoodArray();
		}
		
		VVVdouble & getD2LikelihoodArray(const Node *node)
		{
			return _nodeData[node->getId()].getD2LikelihoodArray();
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
		 * @return The shrunk sub-dataset + indices for the subtree defined by <i>node</i>.
		 */
		virtual SitePatterns * initLikelihoodsWithPatterns(const Node * node, const SiteContainer & sequences, const SubstitutionModel & model) throw (Exception);
		
};

} //end of namespace bpp.

#endif //_DRASRHOMOGENEOUSTREELIKELIHOODDATA_H_

