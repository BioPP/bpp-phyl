// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_DECOMPOSITIONREWARD_H
#define BPP_PHYL_MAPPING_DECOMPOSITIONREWARD_H


#include "DecompositionMethods.h"
#include "Reward.h"

namespace bpp
{
/**
 * @brief Analytical reward using the eigen decomposition method.
 *
 * The codes is adapted from the original R code by Paula Tataru and
 * Asger Hobolth to the formula in the article of Minin & Suchard.
 *
 * Minin, V.N. and Suchard, M.A.,
 * Fast, accurate and simulation-free stochastic mapping
 * Philosophical Transactions of the Royal Society B 2008 363:3985-95.
 *
 * Only reversible models are supported for now.
 *
 * @author Laurent Guéguen
 */
class DecompositionReward :
  public AbstractReward,
  public DecompositionMethods
{
private:
  mutable RowMatrix<double> rewards_;
  mutable double currentLength_;

public:
  DecompositionReward(
      std::shared_ptr<const SubstitutionModelInterface> model,
      std::shared_ptr<const AlphabetIndex1> alphIndex);

  DecompositionReward(
      const StateMapInterface& stateMap,
      std::shared_ptr<const AlphabetIndex1> alphIndex);

  DecompositionReward(const DecompositionReward& dr) :
    AbstractReward(dr),
    DecompositionMethods(dr),
    rewards_(dr.rewards_),
    currentLength_(dr.currentLength_)
  {}

  DecompositionReward& operator=(const DecompositionReward& dr)
  {
    AbstractReward::operator=(dr);
    DecompositionMethods::operator=(dr);

    rewards_        = dr.rewards_;
    currentLength_  = dr.currentLength_;
    return *this;
  }

  virtual ~DecompositionReward() {}

  DecompositionReward* clone() const override { return new DecompositionReward(*this); }

public:
  double getReward(size_t initialState, size_t finalState, double length) const override;

  Matrix<double>* getAllRewards(double length) const override;

  void storeAllRewards(double length, Eigen::MatrixXd& mat) const override;

  /**
   * @brief Set the substitution model.
   *
   * @param model A pointer toward the substitution model to use. Only
   * reversible models are currently supported. Setting a
   * non-reversible model will throw an exception.
   */
  void setSubstitutionModel(std::shared_ptr<const SubstitutionModelInterface> model) override;

protected:
  void initRewards_();

  void computeRewards_(double length) const;

  void alphabetIndexHasChanged() override;

private:
  void fillBMatrice_();
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_DECOMPOSITIONREWARD_H
