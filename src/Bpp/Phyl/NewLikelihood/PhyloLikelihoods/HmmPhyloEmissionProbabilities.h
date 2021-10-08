//
// File: HmmPhyloEmissionProbabilities.h
// Authors:
//   Laurent GuÃ©guen
// Created: mercredi 28 octobre 2015, Ã  17h 58
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef BPP_PHYL_NEWLIKELIHOOD_PHYLOLIKELIHOODS_HMMPHYLOEMISSIONPROBABILITIES_H
#define BPP_PHYL_NEWLIKELIHOOD_PHYLOLIKELIHOODS_HMMPHYLOEMISSIONPROBABILITIES_H

#include <Bpp/Numeric/AbstractParametrizable.h>

#include "../DataFlow/DataFlowCWise.h"
#include "../HmmEmissionProbabilities_Eigen.h"
#include "HmmPhyloAlphabet.h"

namespace bpp
{
using EmissionLogk = CWiseCompound<MatrixLik, ReductionOf<RowLik> >;

/**
 * @brief Emission probabilities in the context of DF phylolikeihoods.
 *
 */

class HmmPhyloEmissionProbabilities :
  public virtual HmmEmissionProbabilities_Eigen,
  public AbstractParametrizable
{
private:
  Context& context_;

  std::shared_ptr<HmmPhyloAlphabet> phylAlph_;

  /*
   *@brief Emission likelihoods are stored in a Matrix from a set of
   * RowVectors.
   *
   */

  ValueRef<MatrixLik> emProb_;

  size_t nbSites_;

public:
  HmmPhyloEmissionProbabilities(std::shared_ptr<HmmPhyloAlphabet> alphabet);

  HmmPhyloEmissionProbabilities(const HmmPhyloEmissionProbabilities& hEP) :
    AbstractParametrizable(hEP),
    context_(hEP.context_),
    phylAlph_(hEP.phylAlph_),
    emProb_(hEP.emProb_),
    nbSites_(hEP.nbSites_)
  {}

  HmmPhyloEmissionProbabilities* clone() const { return new HmmPhyloEmissionProbabilities(*this);}

  const HmmStateAlphabet* getHmmStateAlphabet() const
  {
    return phylAlph_.get();
  }

  size_t getNumberOfStates() const
  {
    return phylAlph_->getNumberOfStates();
  }

  size_t getNumberOfSites() const
  {
    return nbSites_;
  }

  /**
   * @brief Set the new hidden state alphabet.
   * @param stateAlphabet The new state alphabet.
   * @throw UnvalidStateAlphabetException if the new alphabet is uncorrect (for instance is NULL pointer).
   */

  void setHmmStateAlphabet(std::shared_ptr<HmmStateAlphabet> stateAlphabet);

  /**
   * @brief Operator access to the emission probabilities.
   *
   * This is the fastest way to get the values, but no checking is performed on the indices.
   * For debugging purpose, the getPhyloEmissionProbability would be a safer use.
   *
   * @param pos The position of the sequential data to consider.
   * @param state The index of the hidden state to consider, as defined by the HmmStateAlphabet object associated to this class
   *
   */
  DataLik operator()(size_t pos, size_t state) const
  {
    return (emProb_->getTargetValue())(Eigen::Index(state), Eigen::Index(pos));
  }

  ValueRef<MatrixLik> getEmissionProbabilities()
  {
    return emProb_;
  }

  /**
   * @brief Operator access to the emission probabilities.
   *
   * This is the fastest way to get the values, but no checking is performed on the indices.
   * For debugging purpose, the getPhyloEmissionProbability would be a safer use.
   *
   * @param pos The position of the sequential data to consider.
   * @return A vector of probabilities, whose size is the number of hidden states.
   */
  VectorLik operator()(size_t pos) const
  {
    return emProb_->getTargetValue().col(Eigen::Index(pos));
  }

  /**
   * @return The number of positions in the data.
   */
  size_t getNumberOfPositions() const
  {
    return nbSites_;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_NEWLIKELIHOOD_PHYLOLIKELIHOODS_HMMPHYLOEMISSIONPROBABILITIES_H
