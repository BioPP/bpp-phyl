//
// File: DRASDRTreeLikelihoodData.h
// Created by: Julien Dutheil
// Created on: Sat Dec 30 14:20 2006
// From file DRHomogeneousTreeLikelihood.h
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

#ifndef _DRASDRHOMOGENEOUSTREELIKELIHOODDATA_H_
#define _DRASDRHOMOGENEOUSTREELIKELIHOODDATA_H_

#include "AbstractTreeLikelihoodData.h"
#include "SubstitutionModel.h"
#include "PatternTools.h"
#include "SitePatterns.h"

//From SeqLib:
#include <Seq/AlignedSequenceContainer.h>

// From the STL:
#include <map>

using namespace std;

namespace bpp
{

/**
 * @brief Likelihood data structure for a leaf.
 * 
 * This class is for use with the DRASDRTreeLikelihoodData class.
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
#ifndef NO_VIRTUAL_COV
    DRASDRTreeLikelihoodLeafData*
#else
    Clonable*
#endif
    clone() const
    { 
      return new DRASDRTreeLikelihoodLeafData(*this);
    }

  public:
    const Node * getNode() const { return _leaf; }
    void setNode(const Node & node) { _leaf = &node; }

    VVdouble & getLikelihoodArray()  { return _leafLikelihood;  }
};

/**
 * @brief Likelihood data structure for a node.
 * 
 * This class is for use with the DRASDRTreeLikelihoodData class.
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
     *   |------------> Neighbor node of n (id)
      *     |---------> Site i
     *         |------> Rate class c
     *            |---> Ancestral state s
     * </pre>
     * We call this the <i>likelihood array</i> for each node.
     */

    mutable map<int, VVVdouble > _nodeLikelihoods;
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
    //DRASDRTreeLikelihoodNodeData() {}
    //virtual ~DRASDRTreeLikelihoodNodeData() {}

#ifndef NO_VIRTUAL_COV
    DRASDRTreeLikelihoodNodeData*
#else 
    Clonable*
#endif
    clone() const
    { 
      return new DRASDRTreeLikelihoodNodeData(*this);
    }

  public:
    const Node * getNode() const { return _node; }
    
    void setNode(const Node & node) { _node = &node; }

    const map<int, VVVdouble> & getLikelihoodArrays() const { return _nodeLikelihoods; }
    
    map<int, VVVdouble> & getLikelihoodArrays() { return _nodeLikelihoods; }
    
    VVVdouble & getLikelihoodArrayForNeighbor(int neighborId)
    {
      return _nodeLikelihoods[neighborId];
    }
    
    const VVVdouble & getLikelihoodArrayForNeighbor(int neighborId) const
    {
      return _nodeLikelihoods[neighborId];
    }
    
    Vdouble & getDLikelihoodArray() { return _nodeDLikelihoods;  }
    
    const Vdouble & getDLikelihoodArray() const  {  return _nodeDLikelihoods;  }
    
    Vdouble & getD2LikelihoodArray()  {  return _nodeD2Likelihoods; }
    
    const Vdouble & getD2LikelihoodArrayForNeighbor() const  { return _nodeD2Likelihoods; }

    bool isNeighbor(int neighborId) const
    {
      return _nodeLikelihoods.find(neighborId) != _nodeLikelihoods.end();
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
  public AbstractTreeLikelihoodData
{
  protected:

    mutable map<int, DRASDRTreeLikelihoodNodeData> _nodeData;
    mutable map<int, DRASDRTreeLikelihoodLeafData> _leafData;
    mutable VVVdouble _rootLikelihoods;
    mutable VVdouble  _rootLikelihoodsS;
    mutable Vdouble   _rootLikelihoodsSR;

    SiteContainer * _shrunkData;
    unsigned int _nbSites; 
    unsigned int _nbStates;
    unsigned int _nbClasses;
    unsigned int _nbDistinctSites; 

  public:
    DRASDRTreeLikelihoodData(TreeTemplate<Node> & tree, unsigned int nbClasses):
      _nodeData(), _leafData(), _rootLikelihoods(), _rootLikelihoodsS(), _rootLikelihoodsSR(),
      _shrunkData(NULL), _nbSites(0), _nbStates(0), _nbClasses(nbClasses), _nbDistinctSites(0)
    {
      _tree = &tree;
    }

    DRASDRTreeLikelihoodData(const DRASDRTreeLikelihoodData & data):
      AbstractTreeLikelihoodData(data),
      _nodeData(data._nodeData), _leafData(data._leafData),
      _rootLikelihoods(data._rootLikelihoods),
      _rootLikelihoodsS(data._rootLikelihoodsS),
      _rootLikelihoodsSR(data._rootLikelihoodsSR),
      _shrunkData(NULL),
      _nbSites(data._nbSites), _nbStates(data._nbStates),
      _nbClasses(data._nbClasses), _nbDistinctSites(data._nbDistinctSites)
    {
      _tree         = data._tree;
      if(data._shrunkData)
        _shrunkData = dynamic_cast<SiteContainer *>(data._shrunkData->clone());
    }

    DRASDRTreeLikelihoodData & operator=(const DRASDRTreeLikelihoodData & data)
    {
      AbstractTreeLikelihoodData::operator=(data);
      _nodeData          = data._nodeData;
      _leafData          = data._leafData;
      _rootLikelihoods   = data._rootLikelihoods;
      _rootLikelihoodsS  = data._rootLikelihoodsS;
      _rootLikelihoodsSR = data._rootLikelihoodsSR;
      _nbSites           = data._nbSites;
      _nbStates          = data._nbStates;
      _nbClasses         = data._nbClasses;
      _nbDistinctSites   = data._nbDistinctSites;
      _tree              = data._tree;
      if(data._shrunkData)
        _shrunkData      = dynamic_cast<SiteContainer *>(data._shrunkData->clone());
      else
        _shrunkData      = NULL;
      return *this;
    }

    virtual ~DRASDRTreeLikelihoodData() { delete _shrunkData; }

#ifndef NO_VIRTUAL_COV
    DRASDRTreeLikelihoodData*
#else
    Clonable*
#endif
    clone() const { return new DRASDRTreeLikelihoodData(*this); }

  public:
    void setTree(TreeTemplate<Node> & tree)
    { 
      _tree = &tree;
      for(map<int, DRASDRTreeLikelihoodNodeData>::iterator it = _nodeData.begin(); it != _nodeData.end(); it++)
      {
        int id = it->second.getNode()->getId();
        it->second.setNode(*_tree->getNode(id));
      }
      for(map<int, DRASDRTreeLikelihoodLeafData>::iterator it = _leafData.begin(); it != _leafData.end(); it++)
      {
        int id = it->second.getNode()->getId();
        it->second.setNode(*_tree->getNode(id));
      }
    }

    DRASDRTreeLikelihoodNodeData & getNodeData(int nodeId)
    { 
      return _nodeData[nodeId];
    }
    
    const DRASDRTreeLikelihoodNodeData & getNodeData(int nodeId) const
    { 
      return _nodeData[nodeId];
    }
    
    DRASDRTreeLikelihoodLeafData & getLeafData(int nodeId)
    { 
      return _leafData[nodeId];
    }
    
    const DRASDRTreeLikelihoodLeafData & getLeafData(int nodeId) const
    { 
      return _leafData[nodeId];
    }
    
    unsigned int getArrayPosition(int parentId, int sonId, unsigned int currentPosition) const
    {
      return currentPosition;
    }

    const map<int, VVVdouble> & getLikelihoodArrays(int nodeId) const 
    {
      return _nodeData[nodeId].getLikelihoodArrays();
    }
    
    map<int, VVVdouble> & getLikelihoodArrays(int nodeId)
    {
      return _nodeData[nodeId].getLikelihoodArrays();
    }

    VVVdouble & getLikelihoodArray(int parentId, int neighborId)
    {
      return _nodeData[parentId].getLikelihoodArrayForNeighbor(neighborId);
    }
    
    const VVVdouble & getLikelihoodArray(int parentId, int neighborId) const
    {
      return _nodeData[parentId].getLikelihoodArrayForNeighbor(neighborId);
    }
    
    Vdouble & getDLikelihoodArray(int nodeId)
    {
      return _nodeData[nodeId].getDLikelihoodArray();
    }
    
    const Vdouble & getDLikelihoodArray(int nodeId) const
    {
      return _nodeData[nodeId].getDLikelihoodArray();
    }
    
    Vdouble & getD2LikelihoodArray(int nodeId)
    {
      return _nodeData[nodeId].getD2LikelihoodArray();
    }

    const Vdouble & getD2LikelihoodArray(int nodeId) const
    {
      return _nodeData[nodeId].getD2LikelihoodArray();
    }

    VVdouble & getLeafLikelihoods(int nodeId)
    {
      return _leafData[nodeId].getLikelihoodArray();
    }
    
    const VVdouble & getLeafLikelihoods(int nodeId) const
    {
      return _leafData[nodeId].getLikelihoodArray();
    }
    
    VVVdouble & getRootLikelihoodArray() { return _rootLikelihoods; }
    const VVVdouble & getRootLikelihoodArray() const { return _rootLikelihoods; }
    
    VVdouble  & getRootSiteLikelihoodArray() { return _rootLikelihoodsS; }
    const VVdouble  & getRootSiteLikelihoodArray() const { return _rootLikelihoodsS; }
    
    Vdouble   & getRootRateSiteLikelihoodArray() { return _rootLikelihoodsSR; }
    const Vdouble   & getRootRateSiteLikelihoodArray() const { return _rootLikelihoodsSR; }

    unsigned int getNumberOfDistinctSites() const { return _nbDistinctSites; }
    
    unsigned int getNumberOfSites() const { return _nbSites; }
    
    unsigned int getNumberOfStates() const { return _nbStates; }
    
    unsigned int getNumberOfClasses() const { return _nbClasses; }

    const SiteContainer * getShrunkData() const { return _shrunkData; }
    
    /**
     * @brief Resize and initialize all likelihood arrays according to the given data set and substitution model.
     *
     * @param sites The sequences to use as data.
     * @param model The substitution model to use.
     * @throw Exception if an error occures.
     */
    void initLikelihoods(const SiteContainer & sites, const SubstitutionModel & model) throw (Exception);
    
    /**
     * @brief Rebuild likelihood arrays at inner nodes.
     *
     * This method is to be called when the topology of the tree has changed.
     * Node arrays relationship are rebuilt according to the new topology of the tree.
     * The leaves likelihood remain unchanged, so as for the first and second order derivatives.
     */
    void reInit() throw (Exception);
    
#
    void reInit(const Node * node) throw (Exception);

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
     * @param node  The node defining the subtree to analyse.
     * @param sites The sequence container to use.
     * @param model The model, used for initializing leaves' likelihoods.
     */
    void initLikelihoods(const Node * node, const SiteContainer & sites, const SubstitutionModel & model) throw (Exception);

};

} //end of namespace bpp.

#endif //_DRASDRHOMOGENEOUSTREELIKELIHOODDATA_H_

