//
// File: ComputingNode.h
// Created by: Laurent Guéguen
// Created on: vendredi 29 avril 2016, à 15h 44
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

#ifndef _COMPUTINGNODE_H_
#define _COMPUTINGNODE_H_

//From bpp-core:
#include <Bpp/Numeric/VectorTools.h>

#include <Bpp/Clonable.h>

#include "../Tree/PhyloNode.h"
#include "../Tree/AwareNode.h"

namespace bpp
{

  /**
   * @brief An abstract class for specific ComputingNodes, ie Nodes
   * that perform the likelihood computations
   *
   */
  
  class ComputingNode :
    public AwareNode,
    virtual public Clonable
  {
  private:

    /*
     *For compatibilty during development, should be removed.
     */
    
    const RowMatrix<double> voidprob_;
    
  public:
    static const unsigned char D0=0;
    static const unsigned char D1=1;
    static const unsigned char D2=2;

  public:
    ComputingNode() :
      AwareNode(),
      voidprob_()
    {
    }


    ComputingNode(const ComputingNode& nc) :
      AwareNode(nc),
      voidprob_()
    {
    }

    ComputingNode(const PhyloNode& np) :
      AwareNode(np),
      voidprob_()
    {
    }

    ComputingNode& operator=(const ComputingNode& nc) 
    {
      return *this;
    }
    
    virtual ~ComputingNode()
    {
    }

    virtual ComputingNode* clone() const {
      return (new ComputingNode(*this));
    }
    
    /**
     * @brief updates this node.
     *
     * If flag = true (default), node has to be updated (false otherwise).
     */
    
    virtual void update(bool flag = true) {};

    /**
     * @brief checks if this node is up to date.
     */
    
    virtual bool isUp2date() const { return(false);};

    /**
     * @brief Updates this node and all subtree below.
     *
     */

    void updateAllBelow()
    {
      update();
      for (size_t i=0; i<getNumberOfSons();i++)
        getSon(i)->updateAllBelow();
    }

    /*
     * @brief append the vector of the nodes to be updated in the
     * subtree below this node, this one included if needed.
     *
     */
    
    void toBeUpdatedBelow(Vuint& lId) const
    {
      if (!isUp2date())
        lId.push_back(getId());

      if (!hasNoSon()){
        size_t nS=getNumberOfSons();
        for (size_t i=0; i<nS; i++)
          getSon(i)->toBeUpdatedBelow(lId);
      }
    }

    /*
     *@brief from Node
     */

    ComputingNode* getSon(size_t pos)
    {
      return dynamic_cast<ComputingNode*>(AwareNode::getSon(pos));
    }

    const ComputingNode* getSon(size_t pos) const
    {
      return dynamic_cast<const ComputingNode*>(AwareNode::getSon(pos));
    }


    
    /*
     *@brief compute partial likelihoods from likelihoods, upward or
     * downward, using its specific computation.
     *
     */

    /**
     *@brief adds or sets the partial target (D)likelihood from this
     * partial (D)likelihood.
     *
     * @param likelihoods_target a pointer to the partial target
     *  (D)likelihood [in, out].
     * @param likelihoods_node a pointer to the partial (D)likelihood
     *   of this node [in].
     *
     * @param DX tells which likelihood arrays should be computed,
     *   either D0 for likelihoods, D1 for their first derivate, D2
     *   for their second.
     * @param usesLog says if log-likelihoods are used.
     *
     **/
    

    virtual void addUpwardLikelihoodsAtASite(Vdouble* likelihoods_target, const Vdouble* likelihoods_node, unsigned char DX, bool usesLog) const {};
    
    virtual void setUpwardLikelihoodsAtASite(Vdouble* likelihoods_target, const Vdouble* likelihoods_node, unsigned char DX, bool usesLog) const {};
    
    virtual void setDownwardLikelihoodsAtASite(Vdouble* likelihoods_node, const Vdouble* likelihoods_father, unsigned char DX, bool usesLog) const {};
    
    /**
     *@brief adds or sets the  likelihood using its own
     * likelihood.
     *
     * @param likelihoods a pointer to the partial likelihood updated
     *   [in, out].
     *
     * @param likelihoods_self  the partial likelihood of
     *   the gNode.
     *
     * @param patterns the corresponding positions.
     *
     * @param DX tells which likelihood arrays should be computed,
     *   either D0 for likelihoods, D1 for their first derivate, D2
     *   for their second.
     * @param usesLog says if log-likelihoods are used.
     *
     **/

    
    void addUpwardPartialLikelihoods(VVdouble* likelihoods, const VVdouble* likelihoods_self, const std::vector<size_t>& patterns, unsigned char DX, bool usesLog) const
    {
      size_t nbSites=likelihoods->size();
        
      for (size_t i = 0; i < nbSites; i++)
        addUpwardLikelihoodsAtASite(&(*likelihoods)[i], &(*likelihoods_self)[patterns[i]],DX, usesLog);
    }

    
    void setUpwardPartialLikelihoods(VVdouble* likelihoods, const VVdouble* likelihoods_self, const std::vector<size_t>& patterns, unsigned char DX, bool usesLog) const
    {
      size_t nbSites=likelihoods->size();
        
      for (size_t i = 0; i < nbSites; i++)
        setUpwardLikelihoodsAtASite(&(*likelihoods)[i], &(*likelihoods_self)[patterns[i]],DX, usesLog);
    }

    /***
     * Here, these method should not be called.
     *
     **/
    
    virtual const Matrix<double>& getTransitionProbabilities() const
    {
      throw Exception("ComputingNode::getTransitionProbabilities should not be called.");
      
      return voidprob_;
    }

    
    virtual const Matrix<double>& getTransitionProbabilitiesD1() const
    {
      throw Exception("ComputingNode::getTransitionProbabilities should not be called.");
      
      return voidprob_;
    }

    
    virtual const Matrix<double>& getTransitionProbabilitiesD2() const
    {
      throw Exception("ComputingNode::getTransitionProbabilities should not be called.");
      
      return voidprob_;
    }


  };
} // end namespace bpp

#endif // _COMPUTINGNODE_H_

