//
// File: ComputingNode.h
// Created by: Laurent Guéguen
// Created on: lundi 1 juillet 2013, à 15h 17
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
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/VectorTools.h>

#include "../Tree/Node.h"
#include "../Model/SubstitutionModel.h"


namespace bpp
{

  /**
   * @brief An implementation of the Transition Probabilities and
   * method of computing the partial likelihoods at a node.
   *
   * These classes handle convenient arrays for storing previously
   * computed probabilities.
   *
   * The only parameter they have is "scale", and stands for the
   * scaling to be done of the branch length.
   *
   */
  
  class ComputingNode:
    public Node,
    public AbstractParametrizable
  {
  public:
    static const unsigned char D0=0;
    static const unsigned char D1=1;
    static const unsigned char D2=2;
  private:
  
    const SubstitutionModel* model_;

    size_t nbStates_;
    
    /**
     * @brief the scale on the branch
     */

    double scale_;
    
    /**
     * @brief All transition probabilities on the node.
     */

    mutable  RowMatrix<double> probabilities_;
    mutable  RowMatrix<double> probabilitiesD1_;
    mutable  RowMatrix<double> probabilitiesD2_;

    mutable  bool computeProbabilities_;
    mutable  bool computeProbabilitiesD1_;
    mutable  bool computeProbabilitiesD2_;

  private:

    /*
     * @brief for computational purpose
     *
     */
    
    mutable Vdouble vLogStates_;
    
  public:
    ComputingNode(const SubstitutionModel* model);

    ComputingNode(const ComputingNode& np);

    ComputingNode(const Node& np);

    ComputingNode(int num, std::string st);

    ComputingNode();

    ComputingNode& operator=(const ComputingNode& np);

    ComputingNode* clone() const { return new ComputingNode(*this);}
    
    ~ComputingNode() {}

    void setSubstitutionModel(const SubstitutionModel* pSM);

    void setDistanceToFather(double x)
    {
      Node::setDistanceToFather(x);
      update();
    }

    const SubstitutionModel* getSubstitutionModel() const
    {
      return model_;      
    }

    /*
     * @brief return if transition probabilities need to be
     * recomputed.
     *
     */

    bool isUp2dateTransitionProbabilities() const
    {
      return !computeProbabilities_;
    }
    

    bool isUp2dateTransitionProbabilitiesD1() const
    {
      return !computeProbabilitiesD1_;
    }

    bool isUp2dateTransitionProbabilitiesD2() const
    {
      return !computeProbabilitiesD2_;
    }

    /*
     * @brief append the vector of the nodes to be updated in the
     * subtree below this node.
     *
     */
    
    void updatedSubTreeNodes(Vint& lId) const;

    /*
     * @brief compute transition probabilities.
     *
     */
    
    void computeTransitionProbabilities() const;

    void computeTransitionProbabilitiesD1() const;

    void computeTransitionProbabilitiesD2() const;

    /**
     * @brief Return the transition probability between two states.
     *
     */
    
    double getTransitionProbability(size_t i, size_t j) const
    {
      if (computeProbabilities_) 
        computeTransitionProbabilities();

      return probabilities_(i,j);
    }
  
    /**
     * @brief Return the first order derivate of the transition
     * probability between two states.
     *
     */
    
    double getTransitionProbabilityD1(size_t i, size_t j) const
    {
      if (computeProbabilitiesD1_) 
        computeTransitionProbabilitiesD1();

      return probabilitiesD1_(i,j);
    }
  
    /**
     * @brief Return the second order derivate of the transition
     * probability between two states.
     *
     */

    double getTransitionProbabilityD2(size_t i, size_t j) const
    {
      if (computeProbabilitiesD2_) 
        computeTransitionProbabilitiesD2();

      return probabilitiesD2_(i,j);
    }

    /**
     * @brief Return the transition matrix.
     *
     */

    const Matrix<double>& getTransitionProbabilities() const
    {
      if (computeProbabilities_) 
        computeTransitionProbabilities();

      return probabilities_;
    }

    /**
     * @brief Return the first order derivate of the transition
     * matrix.
     *
     */
    
    const Matrix<double>& getTransitionProbabilitiesD1() const
    {
      if (computeProbabilitiesD1_) 
        computeTransitionProbabilitiesD1();

      return probabilitiesD1_;
    }

    /**
     * @brief Return the second order derivate of the transition
     * matrix.
     *
     */
    
    const Matrix<double>& getTransitionProbabilitiesD2() const
    {
      if (computeProbabilitiesD2_) 
        computeTransitionProbabilitiesD2();

      return probabilitiesD2_;
    }

    void fireParameterChanged(const ParameterList& pl);

    /**
     * @brief Sets the computeProbabilities to true on this node.
     *
     * If flag = true (default), node has to be updated (false otherwise).
     */
    
    void update(bool flag = true);

    /**
     * @brief Sets the computeProbabilities to true on this node and
     * all subtree below.
     *
     */

    void updateAll()
    {
      update();
      for (size_t i=0; i<getNumberOfSons();i++)
        getSon(i)->updateAll();
    }
    
    void clearAllModels()
    {
      model_=0;
      for (size_t i=0; i<getNumberOfSons();i++)
        getSon(i)->clearAllModels();
    }

    bool hasModelOnEachNode() const
    {
      for (size_t i=0; i<getNumberOfSons();i++)
        if (! getSon(i)->hasModelOnEachNode())
          return false;
      
      return true;
    }

    /*
     *@brief from Node
     */

    ComputingNode* getSon(size_t pos) throw (IndexOutOfBoundsException)
    {
      return dynamic_cast<ComputingNode*>(Node::getSon(pos));
    }

    const ComputingNode* getSon(size_t pos) const throw (IndexOutOfBoundsException)
    {
      return dynamic_cast<const ComputingNode*>(Node::getSon(pos));
    }

    
    /*
     *@brief compute partial likelihoods from likelihoods
     */

    /**
     *@brief adds or sets the partial target (D)likelihood with this
     * partial (D)likelihood multiplied with the (D)transition
     * probabilities.
     *
     * @param likelihoods_target a pointer to the partial target
     *  (D)likelihood [in, out].
     * @param likelihoods_node a pointer to the partial (D)likelihood
     *   of this node [in].
     * @param DX tells which matrix should be used as used for
     * transition factors, either D0 for transition probabilities, D1
     * for their first derivate, D2 for their second.
     * @param usesLog says if log-likelihoods are used.
     **/
    

    void addUpwardLikelihoodsAtASite(Vdouble* likelihoods_target, const Vdouble* likelihoods_node, unsigned char DX, bool usesLog) const
    {
      double (ComputingNode::*gtP)(size_t,size_t) const = NULL;

      if (DX==D0)
        gtP=&ComputingNode::getTransitionProbability;
      else
        if (DX==D1)
          gtP=&ComputingNode::getTransitionProbabilityD1;
        else
          if (DX==D2)
            gtP=&ComputingNode::getTransitionProbabilityD2;
          else
            throw Exception("ComputingNode::addUpwardLikelihoodsAtASite: unknown function modifier " + TextTools::toString(DX));

      if (usesLog)
        for (size_t x = 0; x < nbStates_; x++)
        {
          for (size_t y = 0; y < nbStates_; y++)
            vLogStates_[y]=log((*this.*gtP)(x, y)) + (*likelihoods_node)[y];

          double v=VectorTools::logSumExp(vLogStates_);
          double w=(*likelihoods_target)[x];
          (*likelihoods_target)[x] = (v<w)?w+log(1+exp(v-w)):v+log(1+exp(w-v));
        }
      else
        for (size_t x = 0; x < nbStates_; x++)
        {
          // For each initial state,
          double likelihood = 0;
          for (size_t y = 0; y < nbStates_; y++)
            likelihood += (*this.*gtP)(x, y) * (*likelihoods_node)[y];
          
          (*likelihoods_target)[x] += likelihood;
        }
    }

    void setUpwardLikelihoodsAtASite(Vdouble* likelihoods_target, const Vdouble* likelihoods_node, unsigned char DX, bool usesLog) const
    {
      double (ComputingNode::*gtP)(size_t,size_t) const = NULL;

      if (DX==D0)
        gtP=&ComputingNode::getTransitionProbability;
      else
        if (DX==D1)
          gtP=&ComputingNode::getTransitionProbabilityD1;
        else
          if (DX==D2)
            gtP=&ComputingNode::getTransitionProbabilityD2;
          else
            throw Exception("ComputingNode::setUpwardLikelihoodsAtASite: unknown function modifier " + TextTools::toString(DX));

      if (usesLog)
        for (size_t x = 0; x < nbStates_; x++)
        {
          // For each initial state,
          for (size_t y = 0; y < nbStates_; y++)
            vLogStates_[y]=log((*this.*gtP)(x, y)) + (*likelihoods_node)[y];
          
          (*likelihoods_target)[x] = VectorTools::logSumExp(vLogStates_);
        }
      else 
        for (size_t x = 0; x < nbStates_; x++)
        {
          // For each initial state,
          double likelihood = 0;
          for (size_t y = 0; y < nbStates_; y++)
            likelihood += (*this.*gtP)(x, y) * (*likelihoods_node)[y];

          (*likelihoods_target)[x] = likelihood;
        }
    }

    
    /**
     *@brief sets the partial (D)likelihood of this node with the
     * partial (D)likelihood of the father multiplied with the
     * (D)transition probabilities.
     *
     * @param likelihoods_node a pointer to the partial (D)likelihood
     * of this node [in, out].
     * @param likelihoods_father a pointer to the partial (D)likelihood
     * of the father node [in].
     * @param DX tells which matrix should be used as used for
     * transition factors, either D0 for transition probabilities, D1
     * for their first derivate, D2 for their second.
     **/
    
    void setDownwardLikelihoodsAtASite(Vdouble* likelihoods_node, const Vdouble* likelihoods_father, unsigned char DX, bool usesLog) const
    {
      double (ComputingNode::*gtP)(size_t,size_t) const = NULL;

      if (DX==D0)
        gtP=&ComputingNode::getTransitionProbability;
      else
        if (DX==D1)
          gtP=&ComputingNode::getTransitionProbabilityD1;
        else
          if (DX==D2)
            gtP=&ComputingNode::getTransitionProbabilityD2;
          else
            throw Exception("ComputingNode::setDownwardLikelihoodsAtASite: unknown function modifier " + TextTools::toString(DX));

      if (usesLog)
        for (size_t x = 0; x < nbStates_; x++)
        {
          // For each initial state,
          for (size_t y = 0; y < nbStates_; y++)
            vLogStates_[y]=log((*this.*gtP)(y, x)) + (*likelihoods_father)[y];
          
          (*likelihoods_node)[x] = VectorTools::logSumExp(vLogStates_);
        }
      else
        for (size_t x = 0; x < nbStates_; x++)
        {
          // For each initial state,
          double likelihood = 0;
          for (size_t y = 0; y < nbStates_; y++)
            likelihood += (*this.*gtP)(y, x) * (*likelihoods_father)[y];
          
          (*likelihoods_node)[x] = likelihood;
        }
    }

    
    /**
     *@brief adds or sets the partial likelihood using its own partial
     * likelihood.
     *
     * @param likelihoods a pointer to the partial likelihood updated
     * [in, out].
     * @param likelihoods_self  the partial likelihood of
     * the ComputingNode.
     * @param DX tells which matrix should be used as used for
     * transition factors, either D0 for transition probabilities, D1
     * for their first derivate, D2 for their second.
     * @param usesLog says if log-likelihoods are used.
     **/


    void addUpwardPartialLikelihoods(VVdouble* likelihoods, const VVdouble* likelihoods_self, unsigned char DX, bool usesLog) const
    {
      size_t nbSites=likelihoods->size();
        
      for (size_t i = 0; i < nbSites; i++)
        addUpwardLikelihoodsAtASite(&(*likelihoods)[i], &(*likelihoods_self)[i],DX, usesLog);
    }

    void setUpwardPartialLikelihoods(VVdouble* likelihoods, const VVdouble* likelihoods_self, unsigned char DX, bool usesLog) const
    {
      size_t nbSites=likelihoods->size();
        
      for (size_t i = 0; i < nbSites; i++)
        setUpwardLikelihoodsAtASite(&(*likelihoods)[i], &(*likelihoods_self)[i],DX, usesLog);
    }

    /**
     *@brief adds or sets the partial likelihood using its own partial
     * likelihoods and a pattern of corresponding positions.
     *
     * @param likelihoods a pointer to the partial likelihood
     * updated [in, out].
     * @param likelihoods_self  the partial likelihood of
     * the ComputingNode.
     * @param patterns the corresponding positions.
     * @param DX tells which matrix should be used as used for
     * transition factors, either D0 for transition probabilities, D1
     * for their first derivate, D2 for their second.
     * @param usesLog says if log-likelihoods are used.
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


    /**
     *@brief sets its own partial likelihood using the father likelihood.
     *
     * @param likelihoods_self a pointer to the partial likelihood of
     * the ComputingNode [in, out].
     * @param likelihoods a pointer to the partial likelihood of the father.
     * @param DX tells which matrix should be used as used for
     * transition factors, either D0 for transition probabilities, D1
     * for their first derivate, D2 for their second.
     * @param usesLog says if log-likelihoods are used.
     **/


    void setDownwardPartialLikelihoods(VVdouble* likelihoods_self, const VVdouble* likelihoods, unsigned char DX, bool usesLog) const
    {
      size_t nbSites=likelihoods->size();
        
      for (size_t i = 0; i < nbSites; i++)
        setDownwardLikelihoodsAtASite(&(*likelihoods_self)[i], &(*likelihoods)[i], DX, usesLog);
    }


  };
} // end namespace bpp

#endif // _COMPUTINGNODE_H_

