//
// File: AbstractLikelihoodNode.h
// Created by: Laurent Guéguen
// Created on: mardi 23 juin 2015, à 00h 26
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

#ifndef _ABSTRACT_RECURSIVE_LIKELIHOOD_NODE_H_
#define _ABSTRACT_RECURSIVE_LIKELIHOOD_NODE_H_

#include "LikelihoodNode.h"
#include "ComputingNode.h"

// From the STL:
#include <map>

namespace bpp
{

/**
 * @brief Basic Likelihood data structure for a node.
 * 
 * This class is for use with the AbstractTreeLikelihoodData class.
 * 
 * Store all conditionnal likelihoods:
 * <pre>
 * x[i][s]
 *   |------> Site i
 *      |---> Ancestral state s
 * </pre> 
 * We call this the <i>likelihood array</i> for each node.
 * In the same way, we store first and second order derivatives.
 *
 * @see AbstractTreeLikelihoodData
 */
  
  class AbstractLikelihoodNode :
    public LikelihoodNode
  {
  protected:

    /*
     * @brief likelihoods and derivatives
     *
     */
    
    mutable VVdouble nodeLikelihoods_;
    mutable VVdouble nodeDLikelihoods_;
    mutable VVdouble nodeD2Likelihoods_;

  public:
    AbstractLikelihoodNode() :
      LikelihoodNode(),
      nodeLikelihoods_(),
      nodeDLikelihoods_(),
      nodeD2Likelihoods_()
    {}
    
    AbstractLikelihoodNode(const Node& np):
      LikelihoodNode(np),
      nodeLikelihoods_(),
      nodeDLikelihoods_(),
      nodeD2Likelihoods_()
    {}

    AbstractLikelihoodNode(int num, std::string st):
      LikelihoodNode(num, st),
      nodeLikelihoods_(),
      nodeDLikelihoods_(),
      nodeD2Likelihoods_()
    {}

    AbstractLikelihoodNode(const AbstractLikelihoodNode& data) :
      LikelihoodNode(data),
      nodeLikelihoods_(data.nodeLikelihoods_),
      nodeDLikelihoods_(data.nodeDLikelihoods_),
      nodeD2Likelihoods_(data.nodeD2Likelihoods_)
    {}
    
    AbstractLikelihoodNode& operator=(const AbstractLikelihoodNode& data)
    {
      LikelihoodNode::operator=(data);
      
      nodeLikelihoods_   = data.nodeLikelihoods_;
      nodeDLikelihoods_  = data.nodeDLikelihoods_;
      nodeD2Likelihoods_ = data.nodeD2Likelihoods_;

      return *this;
    }
 
    virtual AbstractLikelihoodNode* clone() const
    {
      return new AbstractLikelihoodNode(*this);
    }

  public:
    VVdouble& getLikelihoodArray() { return nodeLikelihoods_; }
    const VVdouble& getLikelihoodArray() const { return nodeLikelihoods_; }
    
    VVdouble& getDLikelihoodArray() { return nodeDLikelihoods_; }
    const VVdouble& getDLikelihoodArray() const { return nodeDLikelihoods_; }
    
    VVdouble& getD2LikelihoodArray() { return nodeD2Likelihoods_; }
    const VVdouble& getD2LikelihoodArray() const { return nodeD2Likelihoods_; }

    virtual void resetLikelihoods(size_t nbSites, size_t nbStates)
    {
      nodeLikelihoods_.resize(nbSites);
      nodeDLikelihoods_.resize(nbSites);
      nodeD2Likelihoods_.resize(nbSites);
      
      for (size_t i = 0; i < nbSites; i++)
      {
        Vdouble* likelihoods_i = &(nodeLikelihoods_[i]);
        Vdouble* dLikelihoods_i = &(nodeDLikelihoods_[i]);
        Vdouble* d2Likelihoods_i = &(nodeD2Likelihoods_[i]);
        likelihoods_i->resize(nbStates);
        dLikelihoods_i->resize(nbStates);
        d2Likelihoods_i->resize(nbStates);
        for (size_t s = 0; s < nbStates; s++)
        {
          (*likelihoods_i)[s] = 1; // All likelihoods_ are initialized to 1.
          (*dLikelihoods_i)[s] = 0; // All dLikelihoods_ are initialized to 0.
          (*d2Likelihoods_i)[s] = 0; // All d2Likelihoods_ are initialized to 0.
        }
      }
    }
    

    virtual void resetLikelihoods(unsigned char DX)
    {
      size_t nSites=nodeLikelihoods_.size();
      size_t nbStates=nodeLikelihoods_[0].size();

      double init=(DX==ComputingNode::D0?1:0);
      
      VVdouble* array;
      if (DX==ComputingNode::D0)
        array=&nodeLikelihoods_;
      else
        if (DX==ComputingNode::D1)
          array=&nodeDLikelihoods_;
        else
          if (DX==ComputingNode::D2)
            array=&nodeD2Likelihoods_;
          else
            throw Exception("AbstractLikelihoodNode::resetLikelihoods: unknown function modifier " + TextTools::toString(DX));
        
      
      for (size_t i = 0; i < nSites; i++)
      {
        Vdouble* array_i = &(*array)[i];

        for (size_t s = 0; s < nbStates; s++)
          (*array_i)[s] = init; 

      }
    }
    

    double& operator()(size_t nSite, size_t nState)
    {
      return nodeLikelihoods_[nSite][nState];
    }

    double operator()(size_t nSite, size_t nState) const
    {
      return nodeLikelihoods_[nSite][nState];
    }
    
    Vdouble& operator[](size_t nSite)
    {
      return nodeLikelihoods_[nSite];
    }

    virtual void computeUpwardPartialLikelihoods(const ComputingNode& lNode, unsigned char DX, Vint* vBrid= NULL)
    {
      resetLikelihoods(DX);
    }

    virtual void computeUpwardPartialLikelihoods(const ComputingNode& lNode, const std::vector<const std::vector<size_t>* >& vPatterns, unsigned char DX, Vint* vBrid= NULL)
    {
      resetLikelihoods(DX);
    }

  };

  
} //end of namespace bpp.

#endif //_ABSTRACT_RECURSIVE_LIKELIHOOD_NODE_H_

