// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMPHYLOEMISSIONPROBABILITIES_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMPHYLOEMISSIONPROBABILITIES_H

#include <Bpp/Numeric/AbstractParametrizable.h>

#include "../DataFlow/DataFlowCWise.h"
#include "../HmmEmissionProbabilities_Eigen.h"
#include "HmmPhyloAlphabet.h"

namespace bpp
{
using EmissionLogk = CWiseCompound<MatrixLik, ReductionOf<RowLik>>;

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
   * @throw UnvalidStateAlphabetException if the new alphabet is incorrect (for instance is NULL pointer).
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
    return (emProb_->targetValue())(Eigen::Index(state), Eigen::Index(pos));
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
    return emProb_->targetValue().col(Eigen::Index(pos));
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
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_HMMPHYLOEMISSIONPROBABILITIES_H
