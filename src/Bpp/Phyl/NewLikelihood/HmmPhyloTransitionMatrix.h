//
// File: HmmPhyloTransitionMatrix.h
// Created by: Laurent Guéguen
// Created on: samedi 21 septembre 2013, à 00h 41
//

/*
Copyright or © or Copr. Bio++Development Team, (November 16, 2004)

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

#ifndef _HMMPHYLOTRANSITIONMATRIX_H_
#define _HMMPHYLOTRANSITIONMATRIX_H_

#include "HmmProcessAlphabet.h"

#include <Bpp/Numeric/Prob/Simplex.h>
#include <Bpp/Numeric/Hmm/HmmTransitionMatrix.h>
#include <Bpp/Numeric/VectorTools.h>

namespace bpp
{

/**
 * @brief Describe the transition probabilities between hidden states of a Hidden Markov Model.
 *  The states are denoted the process of a collection. 
 *
 */
  
class HmmPhyloTransitionMatrix:
  public virtual HmmTransitionMatrix,
  public AbstractParametrizable
{
private:
  std::vector<Simplex> vSimplex_;

  const HmmProcessAlphabet* procAlph_;
  
  mutable RowMatrix<double> pij_, tmpmat_;

  mutable Vdouble eqFreq_;

  mutable bool upToDate_;
  
public:

  // HmmPhyloTransitionMatrix(size_t size, const std::string& prefix = "");

  HmmPhyloTransitionMatrix(const HmmProcessAlphabet* procalph, const std::string& prefix = "");

  HmmPhyloTransitionMatrix(const HmmPhyloTransitionMatrix& hptm);

  HmmPhyloTransitionMatrix& operator=(const HmmPhyloTransitionMatrix& hptm);

  HmmPhyloTransitionMatrix* clone() const { return new HmmPhyloTransitionMatrix(*this);}

  /**
   * @return The hidden alphabet associated to this model.
   */

  const HmmStateAlphabet* getHmmStateAlphabet() const
  {
    return procAlph_;
  }

  /**
   * @brief Set the new hidden state alphabet
   * @param stateAlphabet The new state alphabet
   * @throw UnvalidStateAlphabetException if the new alphabet is uncorrect (for instance is NULL pointer).
   */
  
  void setHmmStateAlphabet(const HmmStateAlphabet* stateAlphabet) throw (HmmUnvalidAlphabetException);
   
  /**
   * @return The number of states in the model.
   */

  unsigned int getNumberOfStates() const
  {
    return (unsigned int)vSimplex_.size();
  }

  /**
   * @brief Set the matrix of the transition probabilities.
   *
   */
  
  void setTransitionProbabilities(const Matrix<double>& mat);
  
  /**
   * @brief Get the transition probability between two states.
   *
   * @param i initial state.
   * @param j final state.
   * @return the transition probability between the two states.
   */

  double Pij(unsigned int i, unsigned int j) const
  {
    return vSimplex_[i].prob(j);
  }

  /**
   * @brief Get all transition probabilities as a matrix.
   *
   * @return A n*n matrix will all transition probabilities (n being the number of hidden states).
   */

  const Matrix<double>& getPij() const;

  /**
   * @return The vector of equilibrium frequencies of the Markov chain described by the matrix.
   */

  const std::vector<double>& getEquilibriumFrequencies() const;


  /*
   * @brief From AbstractParametrizable interface
   *
   */

  void fireParameterChanged(const ParameterList& parameters);
  
};

} //end of namespace bpp

#endif //_HMMTRANSITIONMATRIX_H_

