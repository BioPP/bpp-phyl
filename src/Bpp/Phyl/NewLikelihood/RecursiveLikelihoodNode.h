//
// File: RecursiveLikelihoodNode.h
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

#ifndef _RECURSIVE_LIKELIHOOD_NODE_H_
#define _RECURSIVE_LIKELIHOOD_NODE_H_

#include "AbstractLikelihoodNode.h"

#include <Bpp/Numeric/VectorTools.h>

// From the STL:
#include <map>

namespace bpp
{

/**
 * @brief Likelihood data structure for a node.
 * 
 * This class is for use with the RecursiveTreeLikelihoodData class.
 * 
 * Store all conditional likelihoods:
 * <pre>
 * x[i][s]
 *   |------> Site i
 *      |---> Ancestral state s
 * </pre> 
 * We call this the <i>likelihood array</i> for each node.
 * In the same way, we store first and second order derivatives.
 *
 * @see RecursiveTreeLikelihoodData
 */
  
  class RecursiveLikelihoodNode :
    public AbstractLikelihoodNode
  {
  private:

    /*
     * @brief Likelihood of the data below the node, and its derivatives
     *
     */
     
    VVdouble nodeLikelihoods_B_;
    VVdouble nodeDLikelihoods_B_;
    VVdouble nodeD2Likelihoods_B_;

    /*
     * @brief sum of the likelihoods below the node multiplied by the
     * transition probabilities from the father, and its derivatives.
     *
     *  \sum_y P_{x,y}(t) . L(D_{node}|y)
     */
    
    VVdouble node_fatherLikelihoods_B_;
    VVdouble node_fatherDLikelihoods_B_;
    VVdouble node_fatherD2Likelihoods_B_;
    
    /*
     * @brief likelihoods above the node.
     *
     *  L(D^{node}|y)
     */
    
    VVdouble nodeLikelihoods_A_;

    /*
     * @brief allocated vector for computation
     *
     */

    VVdouble temp_;
    VVdouble temp2_;
     
    /*
     * @brief Check if likelihood arrays are up to date
     *
     */
    
    bool up2date_B_;
    bool up2dateD_B_;
    bool up2dateD2_B_;
    bool up2date_BF_;
    bool up2dateD_BF_;
    bool up2dateD2_BF_;
    bool up2date_A_;

  public:
    RecursiveLikelihoodNode() :
      AbstractLikelihoodNode(),
      nodeLikelihoods_B_(),
      nodeDLikelihoods_B_(),
      nodeD2Likelihoods_B_(),
      node_fatherLikelihoods_B_(),
      node_fatherDLikelihoods_B_(),
      node_fatherD2Likelihoods_B_(),
      nodeLikelihoods_A_(),
      temp_(),
      temp2_(),
      up2date_B_(false),
      up2dateD_B_(false),
      up2dateD2_B_(false),
      up2date_BF_(false),
      up2dateD_BF_(false),
      up2dateD2_BF_(false),
      up2date_A_(false)
    {}
    
    RecursiveLikelihoodNode(const Node& np):
      AbstractLikelihoodNode(np),
      nodeLikelihoods_B_(),
      nodeDLikelihoods_B_(),
      nodeD2Likelihoods_B_(),
      node_fatherLikelihoods_B_(),
      node_fatherDLikelihoods_B_(),
      node_fatherD2Likelihoods_B_(),
      nodeLikelihoods_A_(),
      temp_(),
      temp2_(),
      up2date_B_(false),
      up2dateD_B_(false),
      up2dateD2_B_(false),
      up2date_BF_(false),
      up2dateD_BF_(false),
      up2dateD2_BF_(false),
      up2date_A_(false)
    {}

    RecursiveLikelihoodNode(int num, std::string st):
      AbstractLikelihoodNode(num, st),
      nodeLikelihoods_B_(),
      nodeDLikelihoods_B_(),
      nodeD2Likelihoods_B_(),
      node_fatherLikelihoods_B_(),
      node_fatherDLikelihoods_B_(),
      node_fatherD2Likelihoods_B_(),
      nodeLikelihoods_A_(),
      temp_(),
      temp2_(),
      up2date_B_(false),
      up2dateD_B_(false),
      up2dateD2_B_(false),
      up2date_BF_(false),
      up2dateD_BF_(false),
      up2dateD2_BF_(false),
      up2date_A_(false)
    {}
    
    RecursiveLikelihoodNode(const RecursiveLikelihoodNode& data) :
      AbstractLikelihoodNode(data),
      nodeLikelihoods_B_(data.nodeLikelihoods_B_),
      nodeDLikelihoods_B_(data.nodeDLikelihoods_B_),
      nodeD2Likelihoods_B_(data.nodeD2Likelihoods_B_),
      node_fatherLikelihoods_B_(data.node_fatherLikelihoods_B_),
      node_fatherDLikelihoods_B_(data.node_fatherDLikelihoods_B_),
      node_fatherD2Likelihoods_B_(data.node_fatherD2Likelihoods_B_),
      nodeLikelihoods_A_(data.nodeLikelihoods_A_),
      temp_(data.temp_),
      temp2_(data.temp2_),
      up2date_B_(data.up2date_B_),
      up2dateD_B_(data.up2dateD_B_),
      up2dateD2_B_(data.up2dateD2_B_),
      up2date_BF_(data.up2date_BF_),
      up2dateD_BF_(data.up2dateD_BF_),
      up2dateD2_BF_(data.up2dateD2_BF_),
      up2date_A_(data.up2date_A_)
    {}
    
    RecursiveLikelihoodNode& operator=(const RecursiveLikelihoodNode& data)
    {
      AbstractLikelihoodNode::operator=(data);
      
      nodeLikelihoods_B_   = data.nodeLikelihoods_B_;
      nodeDLikelihoods_B_  = data.nodeDLikelihoods_B_;
      nodeD2Likelihoods_B_ = data.nodeD2Likelihoods_B_;
      node_fatherLikelihoods_B_   = data.node_fatherLikelihoods_B_;
      node_fatherDLikelihoods_B_  = data.node_fatherDLikelihoods_B_;
      node_fatherD2Likelihoods_B_ = data.node_fatherD2Likelihoods_B_;
      nodeLikelihoods_A_   = data.nodeLikelihoods_A_;

      up2date_B_ = data.up2date_B_;
      up2dateD_B_ = data.up2dateD_B_;
      up2dateD2_B_ = data.up2dateD2_B_;
      up2date_BF_ = data.up2date_BF_;
      up2dateD_BF_ = data.up2dateD_BF_;
      up2dateD2_BF_ = data.up2dateD2_BF_;
      up2date_A_ = data.up2date_A_;

      temp_ = data.temp_;
      temp2_ = data.temp2_;

      return *this;
    }


    RecursiveLikelihoodNode* clone() const
    {
      return new RecursiveLikelihoodNode(*this);
    }


    void setUseLog(bool useLog)
    {
      if (useLog==usesLog())
        return;
      
      AbstractLikelihoodNode::setUseLog(useLog);

      if (isUp2dateAbove())
      {
        size_t nSites=nodeLikelihoods_A_.size();
        size_t nStates=nodeLikelihoods_A_[0].size();
        
        if (useLog)
        {
          for (size_t i = 0; i < nSites; i++)
          {
          Vdouble* nodeLikelihoods_A_i_ = &(nodeLikelihoods_A_[i]);
          
          for(size_t s = 0; s < nStates; s++)
            (*nodeLikelihoods_A_i_)[s]=log((*nodeLikelihoods_A_i_)[s]);
          }
        }
        else
        {
          for (size_t i = 0; i < nSites; i++)
          {
            Vdouble* nodeLikelihoods_A_i_ = &(nodeLikelihoods_A_[i]);
          
            for(size_t s = 0; s < nStates; s++)
              (*nodeLikelihoods_A_i_)[s]=exp((*nodeLikelihoods_A_i_)[s]);
          }
        }
      }

      if (isUp2dateFatherBelow_(ComputingNode::D0))
      {
        size_t nSites=node_fatherLikelihoods_B_.size();
        size_t nStates=node_fatherLikelihoods_B_[0].size();

        if (useLog)
        {
          for (size_t i = 0; i < nSites; i++)
          {
            Vdouble* node_fatherLikelihoods_B_i_ = &(node_fatherLikelihoods_B_[i]);
          
            for(size_t s = 0; s < nStates; s++)
              (*node_fatherLikelihoods_B_i_)[s]=log((*node_fatherLikelihoods_B_i_)[s]);
          }
        }
        else
        {
          for (size_t i = 0; i < nSites; i++)
          {
            Vdouble* node_fatherLikelihoods_B_i_ = &(node_fatherLikelihoods_B_[i]);
            
            for(size_t s = 0; s < nStates; s++)
            (*node_fatherLikelihoods_B_i_)[s]=exp((*node_fatherLikelihoods_B_i_)[s]);
          }
        }
      }

      if (isUp2dateBelow_(ComputingNode::D0))
      {
        size_t nSites=nodeLikelihoods_B_.size();
        size_t nStates=nodeLikelihoods_B_[0].size();

        if (useLog)
        {
          for (size_t i = 0; i < nSites; i++)
          {
            Vdouble* nodeLikelihoods_B_i_ = &(nodeLikelihoods_B_[i]);
          
            for(size_t s = 0; s < nStates; s++)
              (*nodeLikelihoods_B_i_)[s]=log((*nodeLikelihoods_B_i_)[s]);
          }
        }
        else
        {
          for (size_t i = 0; i < nSites; i++)
          {
            Vdouble* nodeLikelihoods_B_i_ = &(nodeLikelihoods_B_[i]);
            
            for(size_t s = 0; s < nStates; s++)
              (*nodeLikelihoods_B_i_)[s]=exp((*nodeLikelihoods_B_i_)[s]);
          }
        }
      }

    }

    /*
     * @brief Updates Above Likelihood flags, recursively down to the
     * leaves if false and needed.
     *
     */

    bool isUp2dateAbove() const
    {
      return up2date_A_;
    }

    void updateAbove(bool check)
    {
      if (!check)
      {
        if (up2date_A_)
        {
          size_t nS=getNumberOfSons();
          for (size_t i=0; i<nS; i++)
            static_cast<RecursiveLikelihoodNode*>(getSon(i))->updateAbove(false);
        }
        update(false, ComputingNode::D0);
      }
      
      up2date_A_=check;
    } 

    /*
     * @brief reset and reallocate likelihood arrays and Below
     * likelihood arrays.
     *
     */
     
    void resetBelowLikelihoods(size_t nbSites, size_t nbStates, unsigned char DX)
    {
      VVdouble& array=getBelowLikelihoodArray_(DX);
      VVdouble& array2=getToFatherBelowLikelihoodArray_(DX);

      array.resize(nbSites);
      if (hasFather())
        array2.resize(nbSites);

      for (size_t i = 0; i < nbSites; i++)
      {
        (array)[i].resize(nbStates);
        if (hasFather())
          (array2)[i].resize(nbStates);
      }

      updateBelow_(false, DX);
      if (hasFather())
        updateFatherBelow_(false, DX);

      if (DX!=ComputingNode::D0){        
        temp_.resize(nbSites);
        for (size_t i = 0; i < nbSites; i++)
          temp_[i].resize(nbStates);
        if (usesLog())
        {
          temp2_.resize(nbSites);
          for (size_t i = 0; i < nbSites; i++)
            temp2_[i].resize(nbStates);
        }
      }
    }
 
    /*
     * @brief reset the Likelihood for downward recursion
     *
     */
    
    void resetAboveLikelihoods(size_t nbSites, size_t nbStates)
    {
      VVdouble& array=getAboveLikelihoodArray_();

      array.resize(nbSites);
      for (size_t i = 0; i < nbSites; i++)
      {
        (array)[i].resize(nbStates);
      }

      if (temp_.size()!=0)
      {
        temp_.resize(nbSites);
        for (size_t i = 0; i < nbSites; i++)
          temp_[i].resize(nbStates);
      }

      if (temp2_.size()!=0)
      {
        temp2_.resize(nbSites);
        for (size_t i = 0; i < nbSites; i++)
          temp2_[i].resize(nbStates);
      }
      
      updateAbove(false);
    }


    /*
     * @brief COMPUTATIONS
     *
     * @{
     */
    
    /*
     * @brief Compute likelihoods from the product of below & above
     * likelihoods.
     *
     */

    
    void computeLikelihoods(const ComputingNode& cNode, unsigned char DX, Vint* brId = NULL)
    {
      if (!isUp2dateBelow_(DX))
      {
        update(false, DX);
        computeUpwardBelowLikelihoods(cNode, DX, brId);
      }

      if (DX==ComputingNode::D0)
      {
        if (!isUp2dateAbove())
        {
          update(false, DX);
          computeDownwardAboveLikelihoods(cNode);
        }
      }

      if (!isUp2date(DX))
      {
        VVdouble& res=getLikelihoodArray(DX);
        
        res=getBelowLikelihoodArray_(DX);

        if (usesLog())
          res+=getAboveLikelihoodArray_();
        else
          res*=getAboveLikelihoodArray_();
        
        update(true, DX);
      }
    }
    

    /*
     * @brief compute the DXToFatherBelowLikelihoods, where
     * derivated branches ids are given in a pointer of vector<int>.
     *
     */
     
    void computeUpwardToFatherBelowLikelihoods(const ComputingNode& cNode, unsigned char DX, Vint* vBrid= NULL)
    {

      // First check below dependencies are up to date
      
      if (!isUp2dateBelow_(ComputingNode::D0))
        computeUpwardBelowLikelihoods(cNode, ComputingNode::D0, vBrid);

      if (DX!= ComputingNode::D0)
      {
        if (!isUp2dateBelow_(ComputingNode::D1))
          computeUpwardBelowLikelihoods(cNode, ComputingNode::D1, vBrid);

        if (!isUp2dateFatherBelow_(ComputingNode::D0))
          computeUpwardToFatherBelowLikelihoods(cNode, ComputingNode::D0, vBrid);
        
        if (DX== ComputingNode::D2){       
          if (!isUp2dateBelow_(ComputingNode::D2))
            computeUpwardBelowLikelihoods(cNode, ComputingNode::D2, vBrid);
          if (!isUp2dateFatherBelow_(ComputingNode::D1))
            computeUpwardToFatherBelowLikelihoods(cNode, ComputingNode::D1, vBrid);
        }
      }


      // now compute
      
      VVdouble* res=&(getToFatherBelowLikelihoodArray_(DX));

      switch(DX){
      case ComputingNode::D0:
        cNode.setUpwardPartialLikelihoods(res,
                                          &(getBelowLikelihoodArray_(ComputingNode::D0)), ComputingNode::D0, usesLog());
        break;
        
      case ComputingNode::D1:

        cNode.setUpwardPartialLikelihoods(res,
                                          &(getBelowLikelihoodArray_(ComputingNode::D1)),
                                          ComputingNode::D0,
                                          false);

        if (vBrid && VectorTools::contains(*vBrid,getId()))
        {
          if (!usesLog())
            cNode.addUpwardPartialLikelihoods(res,
                                              &(getBelowLikelihoodArray_(ComputingNode::D0)),
                                              ComputingNode::D1,
                                              false);
          else
          {
            temp_=VectorTools::exp(getBelowLikelihoodArray_(ComputingNode::D0));
            cNode.addUpwardPartialLikelihoods(res,
                                              &temp_,
                                              ComputingNode::D1,
                                              false);
          }
        }
        
        break;
        
      case ComputingNode::D2:
        if (vBrid && VectorTools::contains(*vBrid,getId()))
        {
          cNode.setUpwardPartialLikelihoods(res,
                                            &(getBelowLikelihoodArray_(ComputingNode::D1)),
                                            ComputingNode::D1,
                                            false);

          (*res)*=2;

          
          if (!usesLog())
            cNode.addUpwardPartialLikelihoods(res,
                                              &(getBelowLikelihoodArray_(ComputingNode::D0)),
                                              ComputingNode::D2,
                                              false);
          else
          {
            temp_=VectorTools::exp(getBelowLikelihoodArray_(ComputingNode::D0));
            cNode.addUpwardPartialLikelihoods(res,
                                              &temp_,
                                              ComputingNode::D2,
                                              false);
          }

          cNode.addUpwardPartialLikelihoods(res,
                                            &(getBelowLikelihoodArray_(ComputingNode::D2)),
                                            ComputingNode::D0,
                                            false);
        }
        else
          cNode.setUpwardPartialLikelihoods(res,
                                            &(getBelowLikelihoodArray_(ComputingNode::D2)),
                                            ComputingNode::D0,
                                            false);
      }

      updateFatherBelow_(true, DX);
    }

    /*
     * @brief compute the DXBelowLikelihoods, where derivated branches
     * ids are given in a pointer of vector<int>.
     *
     */
     
    void computeUpwardBelowLikelihoods(const ComputingNode& cNode, unsigned char DX, Vint* vBrid= NULL)
    {
      // First check below dependencies are up to date
      size_t nbSons=getNumberOfSons();

      if (nbSons==0)
        return;

      for (size_t l = 0; l < nbSons; l++){
        RecursiveLikelihoodNode* son=dynamic_cast<RecursiveLikelihoodNode*>(getSon(l));
        if (!son->isUp2dateFatherBelow_(ComputingNode::D0))
          son->computeUpwardToFatherBelowLikelihoods(*cNode.getSon(l), ComputingNode::D0, vBrid);
        if (DX!=ComputingNode::D0)
        {
          if (!son->isUp2dateFatherBelow_(ComputingNode::D1))
            son->computeUpwardToFatherBelowLikelihoods(*cNode.getSon(l), ComputingNode::D1, vBrid);
          
          if (DX==ComputingNode::D2){
            if (!son->isUp2dateFatherBelow_(ComputingNode::D2))
              son->computeUpwardToFatherBelowLikelihoods(*cNode.getSon(l), ComputingNode::D2, vBrid);
          }
        }
      }
      
      // now compute
      VVdouble* res=&(getBelowLikelihoodArray_(DX));
      RecursiveLikelihoodNode* son;
      
      switch(DX){
      case ComputingNode::D0:
        
        son=dynamic_cast<RecursiveLikelihoodNode*>(getSon(0));
        if (usesLog()==son->usesLog())
          setLikelihoodsFromSon_(res, &son->getToFatherBelowLikelihoodArray_(ComputingNode::D0), 0);
        else
        {
          if (usesLog())
            temp_=VectorTools::log(son->getToFatherBelowLikelihoodArray_(ComputingNode::D0));
          else
            temp_=VectorTools::exp(son->getToFatherBelowLikelihoodArray_(ComputingNode::D0));
          
          setLikelihoodsFromSon_(res, &temp_, 0);
        }
        
        for (size_t l = 1; l < nbSons; l++)
        {
          son=dynamic_cast<RecursiveLikelihoodNode*>(getSon(l));
          if (usesLog()==son->usesLog())
            if (usesLog())
              addLikelihoodsFromSon_(res, &son->getToFatherBelowLikelihoodArray_(ComputingNode::D0), l);
            else
              multiplyLikelihoodsFromSon_(res, &son->getToFatherBelowLikelihoodArray_(ComputingNode::D0), l);
          else
          {
            if (usesLog())
              temp_=VectorTools::log(son->getToFatherBelowLikelihoodArray_(ComputingNode::D0));
            else
              temp_=VectorTools::exp(son->getToFatherBelowLikelihoodArray_(ComputingNode::D0));
            
            if (usesLog())
              addLikelihoodsFromSon_(res, &temp_, l);
            else
              multiplyLikelihoodsFromSon_(res, &temp_, l);
          }
        }
        
        break;
          
      case ComputingNode::D1:

        for (size_t l = 0; l < nbSons; l++){
          setLikelihoodsFromSon_(&temp_, &dynamic_cast<RecursiveLikelihoodNode*>(getSon(l))->getToFatherBelowLikelihoodArray_(ComputingNode::D1), l);
          
          for (size_t k = 0; k < nbSons; k++) 
            if (k!=l)
            {
              son=dynamic_cast<RecursiveLikelihoodNode*>(getSon(k));
              if (son->usesLog())
              {
                temp2_=VectorTools::exp(son->getToFatherBelowLikelihoodArray_(ComputingNode::D0));
                multiplyLikelihoodsFromSon_(&temp_ , &temp2_, k);
              }
              else
                multiplyLikelihoodsFromSon_(&temp_ , &son->getToFatherBelowLikelihoodArray_(ComputingNode::D0), k);
            }
          
          if (l==0)
            (*res)=temp_;
          else
            (*res)+=temp_;
        }

        break;

      case ComputingNode::D2:

        for (size_t l = 0; l < nbSons; l++){
          setLikelihoodsFromSon_(&temp_, &dynamic_cast<RecursiveLikelihoodNode*>(getSon(l))->getToFatherBelowLikelihoodArray_(ComputingNode::D2), l);
          
          for (size_t k = 0; k < nbSons; k++) 
            if (k!=l)
            {
              son=dynamic_cast<RecursiveLikelihoodNode*>(getSon(k));
              if (son->usesLog())
              {
                temp2_=VectorTools::exp(son->getToFatherBelowLikelihoodArray_(ComputingNode::D0));
                multiplyLikelihoodsFromSon_(&temp_ , &temp2_, k);
              }
              else
                multiplyLikelihoodsFromSon_(&temp_ , &son->getToFatherBelowLikelihoodArray_(ComputingNode::D0), k);
            }
          
          if (l==0)
            (*res)=temp_;
          else
            (*res)+=temp_;
        }
        
        for (size_t l = 0; l < nbSons; l++){
          setLikelihoodsFromSon_(&temp_, &dynamic_cast<RecursiveLikelihoodNode*>(getSon(l))->getToFatherBelowLikelihoodArray_(ComputingNode::D1), l);

          for (size_t l2 = l+1; l2 < nbSons; l2++){
            multiplyLikelihoodsFromSon_(&temp_ , &dynamic_cast<RecursiveLikelihoodNode*>(getSon(l2))->getToFatherBelowLikelihoodArray_(ComputingNode::D1), l2);
            
            for (size_t k = 0; k < nbSons; k++) 
              if ((k!=l) && (k!=l2))
              {
                son=dynamic_cast<RecursiveLikelihoodNode*>(getSon(k));
                if (son->usesLog())
                {
                  temp2_=VectorTools::exp(son->getToFatherBelowLikelihoodArray_(ComputingNode::D0));
                  multiplyLikelihoodsFromSon_(&temp_ , &temp2_, k);
                }
                else
                  multiplyLikelihoodsFromSon_(&temp_ , &son->getToFatherBelowLikelihoodArray_(ComputingNode::D0), k);
              }
            
            temp_*=2;
            
            (*res)+=temp_;
          }
        }
      }

      updateBelow_(true, DX);
    }


    /*
     * @brief set the Above Likelihood to a frequencies (or a set of
     * frequencies). Used for root frequencies.
     *
     */

    void setAboveLikelihoods(const Vdouble& rootFreq)
    {
      VVdouble& abArray=getAboveLikelihoodArray_();

      size_t nbSites=abArray.size();

      if (usesLog_)
        for (size_t i=0; i<nbSites; i++)
          abArray[i]=VectorTools::log(rootFreq);
      else
        for (size_t i=0; i<nbSites; i++)
          abArray[i]=rootFreq;
    }

    void setAboveLikelihoods(const VVdouble& initFreq)
    {
      VVdouble& abArray=getAboveLikelihoodArray_();

      size_t nbSites=abArray.size();

      if (usesLog_)
        for (size_t i=0; i<nbSites; i++)
          abArray[i]=VectorTools::log(initFreq[i]);
      else
        for (size_t i=0; i<nbSites; i++)
          abArray[i]=initFreq[i];
    }

    /*
     * @brief compute the Above Father Likelihoods from the sons and the father.
     *
     */
     
    void computeDownwardAboveLikelihoods(const ComputingNode& cNode)
    {
      // First check dependencies are up to date

      if (hasFather())
      {
        RecursiveLikelihoodNode* father = dynamic_cast<RecursiveLikelihoodNode*>(getFather());
        const ComputingNode* cFather=dynamic_cast<const ComputingNode*>(cNode.getFather());

        if (!father->isUp2dateAbove())
          father->computeDownwardAboveLikelihoods(*cFather);
        
        size_t nbBr = father->getNumberOfSons();
        for (size_t i=0; i<nbBr; i++)
        {
          RecursiveLikelihoodNode* bro= dynamic_cast<RecursiveLikelihoodNode*>(father->getSon(i));
          if (bro!=this)
            if (!bro->isUp2dateFatherBelow_(ComputingNode::D0))
              bro->computeUpwardToFatherBelowLikelihoods(*cFather->getSon(i), ComputingNode::D0);
        }

        if (usesLog()==father->usesLog())
          temp_ = father->getAboveLikelihoodArray_();
        else
          if (usesLog())
            temp_ = VectorTools::log(father->getAboveLikelihoodArray_());
          else
            temp_ = VectorTools::exp(father->getAboveLikelihoodArray_());

        for (size_t i=0; i<nbBr; i++)
        {
          const RecursiveLikelihoodNode* bro= dynamic_cast<const RecursiveLikelihoodNode*>(father->getSon(i));
          
          if (bro!=this)
          {
            if (usesLog())
            {
              if (bro->usesLog())
                temp_+=bro->getToFatherBelowLikelihoodArray_(ComputingNode::D0);
              else
              {
                temp2_=VectorTools::log(bro->getToFatherBelowLikelihoodArray_(ComputingNode::D0));
                temp_+=temp2_;
              }
            }
            else
            {
              if (bro->usesLog())
              {
                temp2_=VectorTools::exp(bro->getToFatherBelowLikelihoodArray_(ComputingNode::D0));
                temp_*=temp2_;
              }
              else
                temp_*=bro->getToFatherBelowLikelihoodArray_(ComputingNode::D0);
            }
          }
        }
              
        cNode.setDownwardPartialLikelihoods(&getAboveLikelihoodArray_(), &temp_, ComputingNode::D0, usesLog());
      }
      
      updateAbove(true);

    }

      private:

    /*
     * @brief retrieve the Below DXLikelihood Arrays
     *
     */
    
    VVdouble& getBelowLikelihoodArray_(unsigned char DX)
    {
      switch(DX){
      case ComputingNode::D0:
        return nodeLikelihoods_B_;
      case ComputingNode::D1:
        return nodeDLikelihoods_B_;
      case ComputingNode::D2:
        return nodeD2Likelihoods_B_;
      default:
        throw Exception("Unknown derivative " + TextTools::toString(DX));
      }
    }

    const VVdouble& getBelowLikelihoodArray_(unsigned char DX) const
    {
      switch(DX){
      case ComputingNode::D0:
        return nodeLikelihoods_B_;
      case ComputingNode::D1:
        return nodeDLikelihoods_B_;
      case ComputingNode::D2:
        return nodeD2Likelihoods_B_;
      default:
        throw Exception("Unknown derivative " + TextTools::toString(DX));
      }
    }

    /*
     * @brief retrieve the Below to Father DXLikelihood Arrays
     *
     */

    VVdouble& getToFatherBelowLikelihoodArray_(unsigned char DX)
    {
      switch(DX){
      case ComputingNode::D0:
        return node_fatherLikelihoods_B_;
      case ComputingNode::D1:
        return node_fatherDLikelihoods_B_;
      case ComputingNode::D2:
        return node_fatherD2Likelihoods_B_;
      default:
        throw Exception("Unknown derivative " + TextTools::toString(DX));
      }
    }

    const VVdouble& getToFatherBelowLikelihoodArray_(unsigned char DX) const
    {
      switch(DX){
      case ComputingNode::D0:
        return node_fatherLikelihoods_B_;
      case ComputingNode::D1:
        return node_fatherDLikelihoods_B_;
      case ComputingNode::D2:
        return node_fatherD2Likelihoods_B_;
      default:
        throw Exception("Unknown derivative " + TextTools::toString(DX));
      }
    }

    /*
     * @brief retrieve the Above to Father Likelihood Arrays
     *
     */

    VVdouble& getAboveLikelihoodArray_() { return nodeLikelihoods_A_; }

    const VVdouble& getAboveLikelihoodArray_() const { return nodeLikelihoods_A_; }


    /*
     * @brief  Use patterns or not for computing likelihood arrays from sons
     *
     */
    
    void multiplyLikelihoodsFromSon_(VVdouble* likelihoods_out, const VVdouble* likelihoods_in, size_t sonNb)
    {
      size_t nbSites=likelihoods_out->size();

      if (vPatt_.size()!=0)
      {
        const std::vector<size_t>& patterns=*vPatt_[sonNb];
        for (size_t i = 0; i < nbSites; i++)
          (*likelihoods_out)[i]*=(*likelihoods_in)[patterns[i]];
      }
      else
        for (size_t i = 0; i < nbSites; i++)
          (*likelihoods_out)[i]*=(*likelihoods_in)[i];
    }

    
    void addLikelihoodsFromSon_(VVdouble* likelihoods_out, const VVdouble* likelihoods_in, size_t sonNb)
    {
      size_t nbSites=likelihoods_out->size();

      if (vPatt_.size()!=0)
      {
        const std::vector<size_t>& patterns=*vPatt_[sonNb];
        for (size_t i = 0; i < nbSites; i++)
          (*likelihoods_out)[i]+=(*likelihoods_in)[patterns[i]];
      }
      else
        for (size_t i = 0; i < nbSites; i++)
          (*likelihoods_out)[i]+=(*likelihoods_in)[i];
    }


    void setLikelihoodsFromSon_(VVdouble* likelihoods_out, const VVdouble* likelihoods_in, size_t sonNb)
    {
      size_t nbSites=likelihoods_out->size();

      if (vPatt_.size()!=0)
      {
        const std::vector<size_t>& patterns=*vPatt_[sonNb];
        for (size_t i = 0; i < nbSites; i++)
          (*likelihoods_out)[i]=(*likelihoods_in)[patterns[i]];
      }
      else
        for (size_t i = 0; i < nbSites; i++)
          (*likelihoods_out)[i]=(*likelihoods_in)[i];
    }
    

    /*
     * @brief Check the Update of Below Likelihood Arrays.
     *
     */
    
    bool isUp2dateBelow_(unsigned char DX) const
    {
      switch(DX){
      case ComputingNode::D0:
        return up2date_B_;
      case ComputingNode::D1:
        return up2dateD_B_;
      case ComputingNode::D2:
        return up2dateD2_B_;
      default:
        return false;
      }
    }
    
    /*
     * @brief Updates Below Likelihood flags, recursively up to the
     * root if false and needed.
     *
     */
      
    void updateBelow_(bool check, unsigned char DX)
    {
      if (!check)
      {
        switch(DX){
        case ComputingNode::D0:
          if (up2date_B_)
          {
            updateFatherBelow_(false, ComputingNode::D0);
          }
        case ComputingNode::D1:
          if (up2dateD_B_)
          {
            updateFatherBelow_(false, ComputingNode::D1);
          }
        case ComputingNode::D2:
          if (up2dateD2_B_)
          {
            updateFatherBelow_(false, ComputingNode::D2);
          }
        }
        switch(DX){
        case ComputingNode::D0:
          up2date_B_=false;
        case ComputingNode::D1:
          up2dateD_B_=false;
        case ComputingNode::D2:
          up2dateD2_B_=false;
        }
        
        update(false, DX);
      }
      else
      {
        switch(DX){
        case ComputingNode::D0:
          up2date_B_=true;
          break;
        case ComputingNode::D1:
          up2dateD_B_=true;
          break;
        case ComputingNode::D2:
          up2dateD2_B_=true;
          break;
        }
      }

    }
   

    bool isUp2dateFatherBelow_(unsigned char DX) const
    {
      switch(DX){
      case ComputingNode::D0:
        return up2date_BF_;
      case ComputingNode::D1:
        return up2dateD_BF_;
      case ComputingNode::D2:
        return up2dateD2_BF_;
      default:
        return false;
      }
    }
    
    /*
     * @brief Updates Below to Father DXLikelihood flags, recursively up
     * to the root if false and needed.
     *
     * Recursivity is done through the BelowLikelihoods only, to make
     * it simple.
     */
      
    void updateFatherBelow_(bool check, unsigned char DX)
    {
      if (!check)
      {
        switch(DX){
        case ComputingNode::D0:
          if (up2date_BF_)
          {
            if (hasFather())
              static_cast<RecursiveLikelihoodNode*>(getFather())->updateBelow_(false, ComputingNode::D0);
          }
        case ComputingNode::D1:
          if (up2dateD_BF_)
          {
            if (hasFather())
              static_cast<RecursiveLikelihoodNode*>(getFather())->updateBelow_(false, ComputingNode::D1);
          }
        case ComputingNode::D2:
          if (up2dateD2_BF_)
          {
            if (hasFather())
              static_cast<RecursiveLikelihoodNode*>(getFather())->updateBelow_(false, ComputingNode::D2);
          }
        }
        switch(DX){
        case ComputingNode::D0:
          up2date_BF_=false;
        case ComputingNode::D1:
          up2dateD_BF_=false;
        case ComputingNode::D2:
          up2dateD2_BF_=false;
        }
        update(false, DX);
      }
      else
      {
        switch(DX){
        case ComputingNode::D0:
          up2date_BF_=true;
          break;
        case ComputingNode::D1:
          up2dateD_BF_=true;
          break;
        case ComputingNode::D2:
          up2dateD2_BF_=true;
          break;
        }
      }
    }

    friend class RecursiveLikelihoodTree;
    friend class RecursiveLikelihoodTreeCalculation;
    
  };

  
} //end of namespace bpp.

#endif //_RECURSIVE_LIKELIHOOD_NODE_H_

