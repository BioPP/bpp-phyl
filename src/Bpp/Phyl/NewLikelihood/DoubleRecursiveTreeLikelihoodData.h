//
// File: DoubleRecursiveTreeLikelihoodData.h
// Created by: Laurent Guéguen
// Created on: Sat Dec 30 14:20 2006
// From file HomogeneousTreeLikelihood.h
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _DOUBLERECURSIVETREELIKELIHOODDATA_H_
#define _DOUBLERECURSIVETREELIKELIHOODDATA_H_

#include "AbstractTreeLikelihoodData.h"
#include "SubstitutionProcess.h"
#include "../SitePatterns.h"

#include <Bpp/Numeric/VectorTools.h>

// From the STL:
#include <map>

namespace bpp
{
  namespace newlik
  {

/**
 * @brief Likelihood data structure for a leaf.
 * 
 * This class is for use with the DoubleRecursiveTreeLikelihoodData
 * class.
 * 
 * Store the likelihoods arrays associated to a leaf.
 * 
 * @see DoubleRecursiveTreeLikelihoodData
 */
    
    class DoubleRecursiveTreeLikelihoodLeafData :
      public virtual TreeLikelihoodNodeData
    {
    private:
      mutable VVdouble leafLikelihood_;
      int leafId_;

    public:
      DoubleRecursiveTreeLikelihoodLeafData() : leafLikelihood_(), leafId_(-1) {}

      DoubleRecursiveTreeLikelihoodLeafData(const DoubleRecursiveTreeLikelihoodLeafData& data) :
        leafLikelihood_(data.leafLikelihood_), leafId_(data.leafId_) {}
    
      DoubleRecursiveTreeLikelihoodLeafData& operator=(const DoubleRecursiveTreeLikelihoodLeafData& data)
      {
        leafLikelihood_ = data.leafLikelihood_;
        leafId_           = data.leafId_;
        return *this;
      }

#ifndef NO_VIRTUAL_COV
      DoubleRecursiveTreeLikelihoodLeafData*
#else
      Clonable*
#endif
      clone() const
      { 
        return new DoubleRecursiveTreeLikelihoodLeafData(*this);
      }

    public:
      int getNodeId() const { return leafId_; }
      void setNodeId(int leafId) { leafId_ = leafId; }

      VVdouble& getLikelihoodArray()  { return leafLikelihood_;  }
    };


/**
 * @brief Likelihood data structure for a node.
 * 
 * This class is for use with the DoubleRecursiveTreeLikelihoodData class.
 * 
 * Store all conditional likelihoods:
 * <pre>
 * x[b][c][i][s]
 *   |------------> Neighbor node of n (id)
 *      |---------> Class c
 *         |------> Site i
 *            |---> Ancestral state s
 * </pre> 
 * We call this the <i>likelihood array</i> for each node.
 *
 *
 * @see DoubleRecursiveTreeLikelihoodData
 */

    class DoubleRecursiveTreeLikelihoodNodeData :
      public virtual TreeLikelihoodNodeData
    {
    private:
      mutable std::map<int, VVVdouble> nodeLikelihoods_;

      /**
       * @brief This contains all likelihood first order derivatives values used for computation.
       *
       * <pre>
       * x[c][i]
       *   |---------> Class c
       *      |------> Site i
       * </pre> 
       * We call this the <i>dLikelihood array</i> for each node.
       */

      mutable VVdouble nodeDLogLikelihoods_;
  
      /**
       * @brief This contains all likelihood second order derivatives values used for computation.
       *
       * <pre>
       * x[c][i]
       *   |---------> Class c
       *      |------> Site i
       * </pre> 
       * We call this the <i>d2Likelihood array</i> for each node.
       */

      mutable VVdouble nodeD2LogLikelihoods_;
    
      int nodeId_;

    public:
      DoubleRecursiveTreeLikelihoodNodeData() :
        nodeLikelihoods_(),
        nodeDLogLikelihoods_(),
        nodeD2LogLikelihoods_(),
        nodeId_(-1)
      {}
    
      DoubleRecursiveTreeLikelihoodNodeData(const DoubleRecursiveTreeLikelihoodNodeData& data) :
        nodeLikelihoods_(data.nodeLikelihoods_),
        nodeDLogLikelihoods_(data.nodeDLogLikelihoods_),
        nodeD2LogLikelihoods_(data.nodeD2LogLikelihoods_),
        nodeId_(data.nodeId_)
      {}
    
      DoubleRecursiveTreeLikelihoodNodeData& operator=(const DoubleRecursiveTreeLikelihoodNodeData& data)
      {
        nodeLikelihoods_   = data.nodeLikelihoods_;
        nodeDLogLikelihoods_  = data.nodeDLogLikelihoods_;
        nodeD2LogLikelihoods_ = data.nodeD2LogLikelihoods_;
        nodeId_            = data.nodeId_;
        return *this;
      }
 
#ifndef NO_VIRTUAL_COV
      DoubleRecursiveTreeLikelihoodNodeData*
#else
      Clonable*
#endif
      clone() const
      {
        return new DoubleRecursiveTreeLikelihoodNodeData(*this);
      }

    public:
      int getNodeId() const { return nodeId_; }
      void setNodeId(int nodeId) { nodeId_ = nodeId; }

      std::map<int, VVVdouble>& getLikelihoodArrays() const { return nodeLikelihoods_; }

      std::map<int, VVVdouble>& getLikelihoodArrays() { return nodeLikelihoods_; }

      VVVdouble& getLikelihoodArrayForNeighbor(int neighborId)
      {
        return nodeLikelihoods_[neighborId];
      }

      const VVVdouble& getLikelihoodArrayForNeighbor(int neighborId) const
      {
        return nodeLikelihoods_[neighborId];
      }
      
      VVdouble& getDLogLikelihoodArray() { return nodeDLogLikelihoods_; }
      const VVdouble& getDLogLikelihoodArray() const { return nodeDLogLikelihoods_; }

      VVdouble& getD2LogLikelihoodArray() { return nodeD2LogLikelihoods_; }
      const VVdouble& getD2LogLikelihoodArray() const { return nodeD2LogLikelihoods_; }

      void eraseNeighborArrays()
      {
        nodeLikelihoods_.erase(nodeLikelihoods_.begin(), nodeLikelihoods_.end());
        nodeDLogLikelihoods_.erase(nodeDLogLikelihoods_.begin(), nodeDLogLikelihoods_.end());
        nodeD2LogLikelihoods_.erase(nodeD2LogLikelihoods_.begin(), nodeD2LogLikelihoods_.end());
      }

      /*
       * @brief Sets all neighbor likelihoods to 1.
       *
       */
      
      void resetLikelihoodArrays()
      {
        std::map<int, VVVdouble >::iterator it;
        for (it=nodeLikelihoods_.begin(); it!=nodeLikelihoods_.end(); it++)
          TreeLikelihoodData::resetLikelihoodArray(it->second);
      }

    };

/**
 * @brief Likelihood data structure suporting simple recursion.
 */
    class DoubleRecursiveTreeLikelihoodData :
      public newlik::AbstractTreeLikelihoodData
    {
    private:
      /**
       * @brief This contains all likelihood values used for computation.
       *
       */

      mutable std::map<int, DoubleRecursiveTreeLikelihoodNodeData> nodeData_;
      mutable std::map<int, DoubleRecursiveTreeLikelihoodLeafData> leafData_;

      /*
       * Store conditional likelihoods:
       * <pre>
       * x[c][i][s]
       *   |---------> Class c
       *      |------> Site i
       *         |---> Ancestral state s
       * </pre> 
       * We call this the <i>likelihood array</i> for the root.
       */
      
      mutable VVVdouble rootLikelihoods_;

      /*
       * Store conditional likelihoods for all classes and sites
       * <pre>
       * x[c][i]
       *   |------> Class c
       *      |---> Site i
       * </pre> 
       * We call this the <i>likelihood array</i> for the root.
       */
      
      mutable VVdouble  rootLikelihoodsS_;

      /*
       * Store conditional likelihoods for all sites:
       * <pre>
       * x[i]
       *   |------> Site i
       * </pre> 
       * We call this the <i>likelihood array</i> for the root.
       */

      mutable Vdouble   rootLikelihoodsSC_;

      std::auto_ptr<SiteContainer> shrunkData_;
      size_t nbSites_; 
      size_t nbStates_;
      size_t nbClasses_;
      size_t nbDistinctSites_; 

      const TreeTemplate<Node>* tree_;

    public:
      DoubleRecursiveTreeLikelihoodData(size_t nbClasses) :
        AbstractTreeLikelihoodData(),
        nodeData_(), leafData_(), rootLikelihoods_(), rootLikelihoodsS_(), rootLikelihoodsSC_(), shrunkData_(0), nbSites_(0), nbStates_(0),
        nbClasses_(nbClasses), nbDistinctSites_(0), tree_(0)
      {}

      DoubleRecursiveTreeLikelihoodData(const DoubleRecursiveTreeLikelihoodData& data):
        AbstractTreeLikelihoodData(data),
        nodeData_(data.nodeData_), leafData_(data.leafData_),
        rootLikelihoods_(data.rootLikelihoods_),
        rootLikelihoodsS_(data.rootLikelihoodsS_),
        rootLikelihoodsSC_(data.rootLikelihoodsSC_),
        shrunkData_(0),
        nbSites_(data.nbSites_), nbStates_(data.nbStates_),
        nbClasses_(data.nbClasses_), nbDistinctSites_(data.nbDistinctSites_), tree_(data.tree_)
      {
        if (data.shrunkData_.get())
          shrunkData_.reset(dynamic_cast<SiteContainer*>(data.shrunkData_->clone()));
      }

      DoubleRecursiveTreeLikelihoodData& operator=(const DoubleRecursiveTreeLikelihoodData & data)
      {
        AbstractTreeLikelihoodData::operator=(data);
        nodeData_          = data.nodeData_;
        leafData_          = data.leafData_;
        rootLikelihoods_   = data.rootLikelihoods_;
        rootLikelihoodsS_  = data.rootLikelihoodsS_;
        rootLikelihoodsSC_ = data.rootLikelihoodsSC_;
        nbSites_           = data.nbSites_;
        nbStates_          = data.nbStates_;
        nbClasses_         = data.nbClasses_;
        nbDistinctSites_   = data.nbDistinctSites_;
        if (data.shrunkData_.get())
          shrunkData_.reset(dynamic_cast<SiteContainer*>(data.shrunkData_->clone()));
        else
          shrunkData_.reset();
        tree_              = data.tree_;
        
        return *this;
      }

      virtual ~DoubleRecursiveTreeLikelihoodData() {}

      DoubleRecursiveTreeLikelihoodData* clone() const { return new DoubleRecursiveTreeLikelihoodData(*this); }

    public:

      /*
       *@brief Set the tree associated to the data.
       *
       * All node data will be actualized accordingly by calling the setNode() method on the corresponding nodes.
       * @warning: the old tree and the new tree must be two clones! And particularly, they have to share the
       * same topology and nodes id.
       *
       * @param tree The tree to be associated to this data.
       */
      
      // void setTree(const TreeTemplate<Node>* tree)
      // { 
      //   tree_ = tree;
      //   for (std::map<int, DoubleRecursiveTreeLikelihoodNodeData>::iterator it = nodeData_.begin(); it != nodeData_.end(); it++)
      //   {
      //     int id = it->second.getNodeId();
      //     it->second.setNodeId(tree_->getNode(id));
      //   }
      //   for (std::map<int, DoubleRecursiveTreeLikelihoodLeafData>::iterator it = leafData_.begin(); it != leafData_.end(); it++)
      //   {
      //     int id = it->second.getNode()->getId();
      //     it->second.setNode(tree_->getNode(id));
      //   }
      // }

      DoubleRecursiveTreeLikelihoodNodeData& getNodeData(int nodeId)
      { 
        return nodeData_[nodeId];
      }
      
      const DoubleRecursiveTreeLikelihoodNodeData& getNodeData(int nodeId) const
      { 
        return nodeData_[nodeId];
      }

      DoubleRecursiveTreeLikelihoodLeafData& getLeafData(int nodeId)
      { 
        return leafData_[nodeId];
      }
    
      const DoubleRecursiveTreeLikelihoodLeafData& getLeafData(int nodeId) const
      { 
        return leafData_[nodeId];
      }
      size_t getArrayPosition(int parentId, int sonId, size_t currentPosition) const
      {
        return currentPosition;
      }
      const std::map<int, VVVdouble>& getLikelihoodArrays(int nodeId) const 
      {
        return nodeData_[nodeId].getLikelihoodArrays();
      }
    
      std::map<int, VVVdouble>& getLikelihoodArrays(int nodeId)
      {
        return nodeData_[nodeId].getLikelihoodArrays();
      }

      VVVdouble& getLikelihoodArray(int parentId, int neighborId)
      {
        return nodeData_[parentId].getLikelihoodArrayForNeighbor(neighborId);
      }
    
      const VVVdouble& getLikelihoodArray(int parentId, int neighborId) const
      {
        return nodeData_[parentId].getLikelihoodArrayForNeighbor(neighborId);
      }
    
      VVdouble& getDLogLikelihoodArray(int nodeId)
      {
        return nodeData_[nodeId].getDLogLikelihoodArray();
      }
    
      const VVdouble& getDLogLikelihoodArray(int nodeId) const
      {
        return nodeData_[nodeId].getDLogLikelihoodArray();
      }
      
      VVdouble& getD2LogLikelihoodArray(int nodeId)
      {
        return nodeData_[nodeId].getD2LogLikelihoodArray();
      }

      const VVdouble& getD2LogLikelihoodArray(int nodeId) const
      {
        return nodeData_[nodeId].getD2LogLikelihoodArray();
      }

      VVdouble& getLeafLikelihoods(int nodeId)
      {
        return leafData_[nodeId].getLikelihoodArray();
      }
    
      const VVdouble& getLeafLikelihoods(int nodeId) const
      {
        return leafData_[nodeId].getLikelihoodArray();
      }
    
      VVVdouble& getRootLikelihoodArray() { return rootLikelihoods_; }
      const VVVdouble & getRootLikelihoodArray() const { return rootLikelihoods_; }
    

      VVdouble& getRootStateLikelihoodArray() { return rootLikelihoodsS_; }
      const VVdouble& getRootStateLikelihoodArray() const { return rootLikelihoodsS_; }
    
      Vdouble& getRootStateClassLikelihoodArray() { return rootLikelihoodsSC_; }
      const Vdouble& getRootStateClassLikelihoodArray() const { return rootLikelihoodsSC_; }

      size_t getNumberOfDistinctSites() const { return nbDistinctSites_; }
      size_t getNumberOfSites() const { return nbSites_; }
      size_t getNumberOfStates() const { return nbStates_; }
      size_t getNumberOfClasses() const { return nbClasses_; }
    
      /**
       * @brief Resize and initialize all likelihood arrays according to the given data set and substitution process.
       *
       * @param sites The sequences to use as data.
       * @param process The substitution process to use.
       * @throw Exception if an error occures.
       */

      void initLikelihoods(const SiteContainer& sites, const SubstitutionProcess& process) throw (Exception);

      /**
       * @brief Rebuild likelihood arrays at inner nodes.
       *
       * This method is to be called when the topology of the tree has changed.
       * Node arrays relationship are rebuilt according to the new topology of the tree.
       * The leaves likelihood remain unchanged, so as for the first and second order derivatives.
       */
      void reInit() throw (Exception);
    
      void reInit(const Node* node) throw (Exception);


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
       * @param process   The substitution process to use.
       */

      void initLikelihoods_(const Node* node, const SiteContainer& sequences, const SubstitutionProcess& process) throw (Exception);

    };

  } //end of namespace newlik.
} //end of namespace bpp.

#endif //_DOUBLERECURSIVETREELIKELIHOODDATA_H_

