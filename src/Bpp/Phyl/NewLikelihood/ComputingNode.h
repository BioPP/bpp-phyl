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

#ifndef _COMPUTINGUNIT_H_
#define _COMPUTINGUNIT_H_

//From the STL:
#include <memory>

//From bpp-core:
#include <Bpp/Numeric/Matrix/Matrix.h>

#include "../Node.h"
#include "../Model/SubstitutionModel.h"

namespace bpp
{

  /**
   * @brief An implementation of the Transition Probabilities at a node.
   *
   * These classes handle convenient arrays for storing previously
   * computed probabilities.
   *
   */
  
  class ComputingNode:
    virtual public Node,
    public AbstractParametrizable
  {
  private:
  
    const SubstitutionModel* model_;

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

    void computeTransitionProbabilities() const;

    void computeTransitionProbabilitiesD1() const;

    void computeTransitionProbabilitiesD2() const;

    double getTransitionProbability(size_t i, size_t j) const
    {
      if (computeProbabilities_) 
        computeTransitionProbabilities();

      return probabilities_(i,j);
    }
  
    double getTransitionProbabilityD1(size_t i, size_t j) const
    {
      if (computeProbabilitiesD1_) 
        computeTransitionProbabilitiesD1();

      return probabilitiesD1_(i,j);
    }
  
    double getTransitionProbabilityD2(size_t i, size_t j) const
    {
      if (computeProbabilitiesD2_) 
        computeTransitionProbabilitiesD2();

      return probabilitiesD2_(i,j);
    }

    const Matrix<double>& getTransitionProbabilities() const
    {
      if (computeProbabilities_) 
        computeTransitionProbabilities();

      return probabilities_;
    }

    const Matrix<double>& getTransitionProbabilitiesD1() const
    {
      if (computeProbabilitiesD1_) 
        computeTransitionProbabilitiesD1();

      return probabilitiesD1_;
    }

    const Matrix<double>& getTransitionProbabilitiesD2() const
    {
      if (computeProbabilitiesD2_) 
        computeTransitionProbabilitiesD2();

      return probabilitiesD2_;
    }

    void fireParameterChanged(const ParameterList& pl);

    /*
     *@brief from Node
     */

    // ComputingNode* getSon(size_t pos) throw (IndexOutOfBoundsException)
    // {
    //   return dynamic_cast<ComputingNode*>(Node::getSon(pos));
    // }
    
  };

 } // end namespace bpp

#endif // _COMPUTINGUNIT_H_

