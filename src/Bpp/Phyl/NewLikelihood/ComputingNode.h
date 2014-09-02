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

//From the STL:
//#include <memory>

//From bpp-core:
#include <Bpp/Numeric/Matrix/Matrix.h>

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
    virtual public Node,
    public AbstractParametrizable
  {
  public:
    static unsigned char D0;
    static unsigned char D1;
    static unsigned char D2;
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
     */
    
    void update();

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
     *@brief compute partial likelihood
     */

    /**
     *@brief multiplies the partial (D)likelihood of the father node
     * with this partial (D)likelihood multiplied with the
     * (D)transition probabilities.
     *
     * @param likelihoods_father a pointer to the partial (D)likelihood
     * of the father node [in, out].
     * @param likelihoods_node a pointer to the partial (D)likelihood
     * of this node [in].
     * @param DX tells which matrix should be used as used for
     * transition factors, either D0 for transition probabilities, D1
     * for their first derivate, D2 for their second.
     **/
    
    void multiplyLikelihoodsAtASite(Vdouble* likelihoods_father, const Vdouble* likelihoods_node, unsigned char DX) const
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
            throw Exception("ComputingNode::multiplyDXLikelihoodsAtASite: unknown function modifier " + TextTools::toString(DX));
        
      for (size_t x = 0; x < nbStates_; x++)
      {
        // For each initial state,
        double likelihood = 0;
        for (size_t y = 0; y < nbStates_; y++)
        {
          likelihood += (*this.*gtP)(x, y) * (*likelihoods_node)[y];
        }
        (*likelihoods_father)[x] *= likelihood;
      }
    }

    /**
     *@brief multiplies its partial likelihood using the partial
     * likelihoods of some sons.
     *
     * @param likelihoods_node a pointer to the partial likelihood
     * of this node [in, out].
     * @param vLikelihoods_sons a vector of the partial likelihoods of
     * the sons. For sons that are not used, these pointers are null.
     * @param DX tells which matrix should be used as used for
     * transition factors, either D0 for transition probabilities, D1
     * for their first derivate, D2 for their second.
     **/

    void multiplyPartialLikelihoods(VVdouble* likelihoods_node, const std::vector<const VVdouble*>& vLikelihoods_sons, unsigned char DX)
    {
      size_t nbSites=likelihoods_node->size();
        
      for (size_t l = 0; l < getNumberOfSons(); l++)
      {
        if (vLikelihoods_sons[l]==0)
          continue;

        const ComputingNode* son_l=static_cast<const ComputingNode*>(getSon(l));
        const VVdouble* likelihoods_son_l = vLikelihoods_sons[l];

        for (size_t i = 0; i < nbSites; i++)
          son_l->multiplyLikelihoodsAtASite(&(*likelihoods_node)[i],&(*likelihoods_son_l)[i],DX);
        
      }
    }

    
    /**
     *@brief multiplies its partial likelihood using the partial
     * likelihoods of the sons.
     *
     * @param likelihoods_node a pointer to the partial likelihood
     * of this node [in, out].
     * @param vLikelihoods_sons a vector of the partial likelihoods of
     * some sons. For sons that are not used, these pointers are null.
     * @param vPatterns a vector of the corresponding positions
     * from this node to the sons.
     * @param DX tells which matrix should be used as used for
     * transition factors, either D0 for transition probabilities, D1
     * for their first derivate, D2 for their second.
     **/

    void multiplyPartialLikelihoods(VVdouble* likelihoods_node, const std::vector<const VVdouble*>& vLikelihoods_sons, const std::vector<const std::vector<size_t>* >& vPatterns, unsigned char DX)
    {
      size_t nbSites=likelihoods_node->size();
        
      for (size_t i = 0; i < nbSites; i++)
      {
        Vdouble* likelihoods_node_i = &(*likelihoods_node)[i];
      
        for (size_t l = 0; l < getNumberOfSons(); l++)
        {
          if (vLikelihoods_sons[l]==0)
            continue;

          const ComputingNode* son_l=dynamic_cast<const ComputingNode*>(getSon(l));
                    
          const Vdouble* likelihoods_son_l_i = &(*vLikelihoods_sons[l])[(*vPatterns[l])[i]];
          son_l->multiplyLikelihoodsAtASite(likelihoods_node_i,likelihoods_son_l_i,DX);
        }
      }
    }


    /**
     *@brief multiplies the partial likelihood using the partial
     * likelihoods of a son.
     *
     * @param likelihoods_node a pointer to the partial likelihood
     * of this node [in, out].
     * @param likelihoods_son  the partial likelihood of
     * the used son.
     * @param son the used Computing son.
     * @param DX tells which matrix should be used as used for
     * transition factors, either D0 for transition probabilities, D1
     * for their first derivate, D2 for their second.
     **/

    void multiplyPartialLikelihoods(VVdouble* likelihoods_node, const VVdouble* likelihoods_son, const ComputingNode* son, unsigned char DX)
    {
      size_t nbSites=likelihoods_node->size();
        
      for (size_t i = 0; i < nbSites; i++)
        son->multiplyLikelihoodsAtASite(&(*likelihoods_node)[i], &(*likelihoods_son)[i],DX);
    }

    /**
     *@brief multiplies the partial likelihood using the partial
     * likelihoods of a son.
     *
     * @param likelihoods_node a pointer to the partial likelihood
     * of this node [in, out].
     * @param likelihoods_son  the partial likelihood of
     * the used son.
     * @param son the used Computing son.
     * @param patterns the corresponding positions from this node to
     * the son.
     * @param DX tells which matrix should be used as used for
     * transition factors, either D0 for transition probabilities, D1
     * for their first derivate, D2 for their second.
     **/

    void multiplyPartialLikelihoods(VVdouble* likelihoods_node, const VVdouble* likelihoods_son, const ComputingNode* son, const std::vector<size_t>& patterns, unsigned char DX)
    {
      size_t nbSites=likelihoods_node->size();
        
      for (size_t i = 0; i < nbSites; i++)
        son->multiplyLikelihoodsAtASite(&(*likelihoods_node)[i], &(*likelihoods_son)[patterns[i]],DX);
    }
    
  };
} // end namespace bpp

#endif // _COMPUTINGNODE_H_

