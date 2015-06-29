//
// File: SingleRecursiveLikelihoodNode.h
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

#ifndef _SINGLE_RECURSIVE_LIKELIHOOD_NODE_H_
#define _SINGLE_RECURSIVE_LIKELIHOOD_NODE_H_

#include "AbstractLikelihoodNode.h"

#include <Bpp/Numeric/VectorTools.h>

// From the STL:
#include <map>

namespace bpp
{

/**
 * @brief Likelihood data structure for a node.
 * 
 * This class is for use with the SingleRecursiveTreeLikelihoodData class.
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
 * @see SingleRecursiveTreeLikelihoodData
 */
  
  class SingleRecursiveLikelihoodNode :
    public AbstractLikelihoodNode
  {
  private:

    /*
     * @brief sum of the likelihoods under the node multiplied by the
     * transition probabilities from the father, and its derivatives.
     *
     *  \sum_y P_{x,y}(t) . L(D_n|y)
     */
    
    mutable VVdouble node_fatherLikelihoods_;
    mutable VVdouble node_fatherDLikelihoods_;
    mutable VVdouble node_fatherD2Likelihoods_;

  public:
    SingleRecursiveLikelihoodNode() :
      AbstractLikelihoodNode(),
      node_fatherLikelihoods_(),
      node_fatherDLikelihoods_(),
      node_fatherD2Likelihoods_()
    {}
    
    SingleRecursiveLikelihoodNode(const Node& np):
      AbstractLikelihoodNode(np),
      node_fatherLikelihoods_(),
      node_fatherDLikelihoods_(),
      node_fatherD2Likelihoods_()
    {}

    SingleRecursiveLikelihoodNode(int num, std::string st):
      AbstractLikelihoodNode(num, st),
      node_fatherLikelihoods_(),
      node_fatherDLikelihoods_(),
      node_fatherD2Likelihoods_()
    {}
    
    SingleRecursiveLikelihoodNode(const SingleRecursiveLikelihoodNode& data) :
      AbstractLikelihoodNode(data),
      node_fatherLikelihoods_(data.node_fatherLikelihoods_),
      node_fatherDLikelihoods_(data.node_fatherDLikelihoods_),
      node_fatherD2Likelihoods_(data.node_fatherD2Likelihoods_)
    {}
    
    SingleRecursiveLikelihoodNode& operator=(const SingleRecursiveLikelihoodNode& data)
    {
      AbstractLikelihoodNode::operator=(data);
      
      node_fatherLikelihoods_   = data.node_fatherLikelihoods_;
      node_fatherDLikelihoods_  = data.node_fatherDLikelihoods_;
      node_fatherD2Likelihoods_ = data.node_fatherD2Likelihoods_;
      return *this;
    }
 
    SingleRecursiveLikelihoodNode* clone() const
    {
      return new SingleRecursiveLikelihoodNode(*this);
    }

    VVdouble& getToFatherLikelihoodArray() { return node_fatherLikelihoods_; }
    const VVdouble& getToFatherLikelihoodArray() const { return node_fatherLikelihoods_; }
    
    VVdouble& getToFatherDLikelihoodArray() { return node_fatherDLikelihoods_; }
    const VVdouble& getToFatherDLikelihoodArray() const { return node_fatherDLikelihoods_; }
    
    VVdouble& getToFatherD2LikelihoodArray() { return node_fatherD2Likelihoods_; }
    const VVdouble& getToFatherD2LikelihoodArray() const { return node_fatherD2Likelihoods_; }


    void resetLikelihoods(size_t nbSites, size_t nbStates)
    {
      AbstractLikelihoodNode::resetLikelihoods(nbSites, nbStates);
      node_fatherLikelihoods_.resize(nbSites);
      node_fatherDLikelihoods_.resize(nbSites);
      node_fatherD2Likelihoods_.resize(nbSites);
      
      for (size_t i = 0; i < nbSites; i++)
      {
        Vdouble* _fatherlikelihoods_i = &(node_fatherLikelihoods_[i]);
        Vdouble* _fatherdLikelihoods_i = &(node_fatherDLikelihoods_[i]);
        Vdouble* _fatherd2Likelihoods_i = &(node_fatherD2Likelihoods_[i]);
        _fatherlikelihoods_i->resize(nbStates);
        _fatherdLikelihoods_i->resize(nbStates);
        _fatherd2Likelihoods_i->resize(nbStates);
        for (size_t s = 0; s < nbStates; s++)
        {
          (*_fatherlikelihoods_i)[s] = 0; // All _fatherlikelihoods_ are initialized to 0.
          (*_fatherdLikelihoods_i)[s] = 0; // All _fatherdLikelihoods_ are initialized to 0.
          (*_fatherd2Likelihoods_i)[s] = 0; // All _fatherd2Likelihoods_ are initialized to 0.
        }
      }
    }

    void resetLikelihoods(unsigned char DX)
    {
      AbstractLikelihoodNode::resetLikelihoods(DX);
      
      size_t nSites=nodeLikelihoods_.size();
      size_t nbStates=nodeLikelihoods_[0].size();

      double init=0;
      
      VVdouble* array;
      if (DX==ComputingNode::D0)
        array=&node_fatherLikelihoods_;
      else
        if (DX==ComputingNode::D1)
          array=&node_fatherDLikelihoods_;
        else
          if (DX==ComputingNode::D2)
            array=&node_fatherD2Likelihoods_;
          else
            throw Exception("SingleRecursiveLikelihoodNode::resetLikelihoods: unknown function modifier " + TextTools::toString(DX));
        
      
      for (size_t i = 0; i < nSites; i++)
      {
        Vdouble* array_i = &(*array)[i];

        for (size_t s = 0; s < nbStates; s++)
          (*array_i)[s] = init; 
      }
    }

    
    /*
     * @brief
     *
     * compute the PartialLikelihoods from the sons and perform as if
     * all sons branches are under DX
     *
     */
     
    void computeUpwardPartialLikelihoods(const ComputingNode& cNode, unsigned char DX, Vint* vBrid= NULL)
    {
      resetLikelihoods(DX);

      size_t nbNodes=getNumberOfSons();
      
      if (DX==ComputingNode::D0)
        for (size_t l = 0; l < nbNodes; l++){
          SingleRecursiveLikelihoodNode* son=dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l));

          cNode.getSon(l)->setUpwardPartialLikelihoods(&(son->getToFatherLikelihoodArray()),
                                                        &(son->getLikelihoodArray()), ComputingNode::D0);
          nodeLikelihoods_*=son->getToFatherLikelihoodArray();
          
        }
      else if (DX==ComputingNode::D1){
        for (size_t l = 0; l < nbNodes; l++){
          const ComputingNode* cSon = cNode.getSon(l);
          SingleRecursiveLikelihoodNode* son=dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l));

          cSon->setUpwardPartialLikelihoods(&(son->getToFatherDLikelihoodArray()),
                                            &(son->getDLikelihoodArray()), ComputingNode::D0);

          if (vBrid && VectorTools::contains(*vBrid,cSon->getId()))
            cSon->addUpwardPartialLikelihoods(&(son->getToFatherDLikelihoodArray()),
                                              &(son->getLikelihoodArray()), ComputingNode::D1);
        }
        for (size_t l = 0; l < nbNodes; l++){
          getToFatherDLikelihoodArray()= dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l))->getToFatherDLikelihoodArray(); // as
                                                                                                                                 // temp//variable
          //!!!!!! Check if real copy or not
          
          for (size_t k = 0; k < nbNodes; k++) 
            if (k!=l)
              getToFatherDLikelihoodArray() *= dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(k))->getToFatherLikelihoodArray();
          nodeDLikelihoods_+=getToFatherDLikelihoodArray();
        }
      }
      else if (DX==ComputingNode::D2){
        for (size_t l = 0; l < nbNodes; l++){
          const ComputingNode* cSon = cNode.getSon(l);
          SingleRecursiveLikelihoodNode* son=dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l));

          cSon->setUpwardPartialLikelihoods(&(son->getToFatherD2LikelihoodArray()),
                                            &(son->getD2LikelihoodArray()), ComputingNode::D0);

          if (vBrid && VectorTools::contains(*vBrid,cSon->getId())) {
            cSon->setUpwardPartialLikelihoods(&getToFatherD2LikelihoodArray(), //as
                                              //temp variable
                                              &(son->getDLikelihoodArray()), ComputingNode::D1);
            getToFatherD2LikelihoodArray()*=2;
            son->getToFatherD2LikelihoodArray()+=getToFatherD2LikelihoodArray();
            cSon->addUpwardPartialLikelihoods(&son->getToFatherD2LikelihoodArray(),
                                              &(son->getLikelihoodArray()), ComputingNode::D2);
          }
        }
        for (size_t l = 0; l < nbNodes; l++){
          getToFatherD2LikelihoodArray()= dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l))->getToFatherD2LikelihoodArray();
          
          for (size_t k = 0; k < nbNodes; k++) 
            if (k!=l)
              getToFatherD2LikelihoodArray() *= dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(k))->getToFatherLikelihoodArray();
          nodeD2Likelihoods_+=getToFatherD2LikelihoodArray();
        }
        for (size_t l = 0; l < nbNodes; l++){
          getToFatherD2LikelihoodArray()= dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l))->getToFatherDLikelihoodArray();
          for (size_t l2 = l+1; l2 < nbNodes; l2++){
            getToFatherD2LikelihoodArray() *= dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l2))->getToFatherDLikelihoodArray();
            
            for (size_t k = 0; k < nbNodes; k++) 
              if ((k!=l) && (k!=l2))
                getToFatherD2LikelihoodArray() *= dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(k))->getToFatherLikelihoodArray();
            getToFatherD2LikelihoodArray()*=2;
            
            nodeD2Likelihoods_+=getToFatherD2LikelihoodArray();
          }
        }
      }
    }
    
    void computeUpwardPartialLikelihoods(const ComputingNode& cNode, const std::vector<const std::vector<size_t>* >& vPatterns, unsigned char DX,  Vint* vBrid= NULL)
    {
      resetLikelihoods(DX);
      
      size_t nbNodes=getNumberOfSons();

      if (DX==ComputingNode::D0)
        for (size_t l = 0; l < nbNodes; l++){
          SingleRecursiveLikelihoodNode* son=dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l));

          cNode.getSon(l)->setUpwardPartialLikelihoods(&(son->getToFatherLikelihoodArray()),
                                                        &(son->getLikelihoodArray()),
                                                        ComputingNode::D0);

          ComputingNode::multiplyPartialLikelihoods(&nodeLikelihoods_,&son->getToFatherLikelihoodArray(), *vPatterns[l]);
        }
      
      else if (DX==ComputingNode::D1){
        for (size_t l = 0; l < nbNodes; l++){
          const ComputingNode* cSon = cNode.getSon(l);
          SingleRecursiveLikelihoodNode* son=dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l));

          cSon->setUpwardPartialLikelihoods(&(son->getToFatherDLikelihoodArray()),
                                            &(son->getDLikelihoodArray()), ComputingNode::D0);

          if (vBrid && VectorTools::contains(*vBrid,cSon->getId()))
            cSon->addUpwardPartialLikelihoods(&(son->getToFatherDLikelihoodArray()),
                                              &(son->getLikelihoodArray()), ComputingNode::D1);
        }
        for (size_t l = 0; l < nbNodes; l++){
          SingleRecursiveLikelihoodNode* son=dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l));
          ComputingNode::setPartialLikelihoods(&getToFatherDLikelihoodArray(),
                                               &son->getToFatherDLikelihoodArray(),
                                               *vPatterns[l]);
          // as temp variable
          
          for (size_t k = 0; k < nbNodes; k++) 
            if (k!=l){
              SingleRecursiveLikelihoodNode* sonk=dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(k));
              ComputingNode::multiplyPartialLikelihoods(&getToFatherDLikelihoodArray(),
                                                        &sonk->getToFatherLikelihoodArray(),
                                                        *vPatterns[k]);
            }
            
          nodeDLikelihoods_+=getToFatherDLikelihoodArray();
        }
      }
      else if (DX==ComputingNode::D2){
        for (size_t l = 0; l < nbNodes; l++){
          const ComputingNode* cSon = cNode.getSon(l);
          SingleRecursiveLikelihoodNode* son=dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l));

          cSon->setUpwardPartialLikelihoods(&(son->getToFatherD2LikelihoodArray()),
                                            &(son->getD2LikelihoodArray()), ComputingNode::D0);

          if (vBrid && VectorTools::contains(*vBrid,cSon->getId())) {
            cSon->setUpwardPartialLikelihoods(&getToFatherD2LikelihoodArray(), //as
                                              //temp variable
                                              &(son->getDLikelihoodArray()), ComputingNode::D1);
            getToFatherD2LikelihoodArray()*=2;
            son->getToFatherD2LikelihoodArray()+=getToFatherD2LikelihoodArray();
            cSon->addUpwardPartialLikelihoods(&son->getToFatherD2LikelihoodArray(),
                                              &(son->getLikelihoodArray()), ComputingNode::D2);
          }
        }
        for (size_t l = 0; l < nbNodes; l++){
          SingleRecursiveLikelihoodNode* son=dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l));
          ComputingNode::setPartialLikelihoods(&getToFatherD2LikelihoodArray(),
                                               &son->getToFatherD2LikelihoodArray(),
                                               *vPatterns[l]);
            
          for (size_t k = 0; k < nbNodes; k++) 
            if (k!=l){
              SingleRecursiveLikelihoodNode* sonk=dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(k));

              ComputingNode::multiplyPartialLikelihoods(&getToFatherD2LikelihoodArray(),
                                                        &(sonk->getToFatherLikelihoodArray()),
                                                        *vPatterns[k]);
            }
            
          nodeD2Likelihoods_+=getToFatherD2LikelihoodArray();
        }
        for (size_t l = 0; l < nbNodes; l++){
          SingleRecursiveLikelihoodNode* son=dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l));

          ComputingNode::setPartialLikelihoods(&getToFatherD2LikelihoodArray(),
                                               &(son->getToFatherDLikelihoodArray()),
                                               *vPatterns[l]);

          for (size_t l2 = l+1; l2 < nbNodes; l2++){
            SingleRecursiveLikelihoodNode* son2=dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(l2));
            ComputingNode::multiplyPartialLikelihoods(&getToFatherD2LikelihoodArray(),
                                                      &(son2->getToFatherDLikelihoodArray()),
                                                      *vPatterns[l2]);
            
            for (size_t k = 0; k < nbNodes; k++) 
              if ((k!=l) && (k!=l2)){
                SingleRecursiveLikelihoodNode* sonk=dynamic_cast<SingleRecursiveLikelihoodNode*>(getSon(k));
                  
                ComputingNode::multiplyPartialLikelihoods(&getToFatherD2LikelihoodArray(),
                                                          &(sonk->getToFatherLikelihoodArray()),
                                                          *vPatterns[k]);
              }
              
            getToFatherD2LikelihoodArray()*=2;
            
            nodeD2Likelihoods_+=getToFatherD2LikelihoodArray();
          }
        }
      }

    }

    
  };

  
} //end of namespace bpp.

#endif //_SINGLE_RECURSIVE_LIKELIHOOD_NODE_H_

