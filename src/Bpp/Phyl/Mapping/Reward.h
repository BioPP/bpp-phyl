// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_REWARD_H
#define BPP_PHYL_MAPPING_REWARD_H

#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex1.h>

#include "../Model/SubstitutionModel.h"

// From the STL:
#include <vector>

namespace bpp
{
/**
 * @brief The Reward interface.
 *
 * Provide a method to compute the reward of a real-valued function
 * @f$f@f$ of the states on a branch. Namely, on a branch of length
 * @f$t@f$, with initial state @f$x@f$ and final state @f$y@f$, if
 * @f$t_s@f$ is the time spent in each state @f$s@f$, the reward of
 * @f$f@f$ is @f$\sum_{s} f(s).E(t_s) @f$.
 *
 * @author Laurent GuÃÂ©guen
 *
 * See:
 * Minin, V.N. and Suchard, M.A.,
 * Fast, accurate and simulation-free stochastic mapping
 * Philosophical Transactions of the Royal Society B 2008 363:3985-95.
 */
class Reward :
  public virtual Clonable
{
public:
  Reward() {}
  virtual ~Reward() {}
  virtual Reward* clone() const = 0;

public:
  /**
   * @return Tell if an alphabet index  has been attached to this class.
   */
  virtual bool hasAlphabetIndex() const = 0;

  /**
   * @return The AlphabetIndex1 object associated to this instance.
   * The alphabet index contains the value associated to each state.
   */
  virtual std::shared_ptr<const AlphabetIndex1> getAlphabetIndex() const = 0;

  /**
   * @param alphind The new AlphabetIndex1 object to be associated to this instance.
   */

  virtual void setAlphabetIndex(std::shared_ptr<const AlphabetIndex1> alphind) = 0;

  /**
   * @brief Short cut function, equivalent to getSubstitutionRegister()->getAlphabet().
   *
   * @return The alphabet associated to this substitution count.
   */
  virtual std::shared_ptr<const Alphabet> getAlphabet() const
  {
    return getAlphabetIndex()->getAlphabet();
  }

  /**
   * @brief Short cut function, equivalent to getSubstitutionRegister()->getAlphabet()->getSize().
   *
   * @return The number of states in the model/alphabet.
   */
  virtual size_t getNumberOfStates() const { return getAlphabet()->getSize(); }


  /**
   * @brief Get the reward of susbstitutions on a branch, given the initial and final states, and the branch length.
   *
   * @param initialState The initial state.
   * @param finalState   The final state.
   * @param length       The length of the branch.
   * @return The reward of the function on a branch of specified length and
   * according to initial and final states.
   */
  virtual double getReward(size_t initialState, size_t finalState, double length) const = 0;

  /**
   * @brief Get the rewards on a branch, for each initial and final
   * states, and given the branch length.
   *
   * @param length       The length of the branch.
   * @return A matrix with all rewards for each initial and final states.
   */
  virtual Matrix<double>* getAllRewards(double length) const = 0;

  /**
   * @brief Store the rewards on a branch, for each initial and final
   * states, and given the branch length.
   *
   * @param length       The length of the branch.
   * @param mat A matrix to store  all rewards for each initial and final states.
   */
  virtual void storeAllRewards(double length, Eigen::MatrixXd& mat) const = 0;

  /**
   * @brief Set the substitution model associated with this reward, if relevant.
   *
   * @param model The substitution model to use with this reward.
   */
  virtual void setSubstitutionModel(std::shared_ptr<const SubstitutionModelInterface> model) = 0;
};


/**
 * @brief Basic implementation of the the Reward interface.
 *
 * This partial implementation deals with the AlphabetIndex1
 * gestion, by owning a pointer.
 */
class AbstractReward :
  public virtual Reward
{
protected:
  std::shared_ptr<const AlphabetIndex1> alphIndex_;

public:
  AbstractReward(std::shared_ptr<const AlphabetIndex1> alphIndex) :
    alphIndex_(alphIndex)
  {}

  AbstractReward(const AbstractReward& ar) :
    alphIndex_(ar.alphIndex_)
  {}

  AbstractReward& operator=(const AbstractReward& ar)
  {
    alphIndex_ = ar.alphIndex_;
    return *this;
  }

  virtual ~AbstractReward() {}

public:
  bool hasAlphabetIndex() const { return alphIndex_ != 0; }

  /**
   * @brief attribution of an AlphabetIndex1
   *
   * @param alphIndex pointer to a AlphabetIndex1
   */
  void setAlphabetIndex(std::shared_ptr<const AlphabetIndex1> alphIndex)
  {
    alphIndex_ = alphIndex;
    alphabetIndexHasChanged();
  }

  std::shared_ptr<const AlphabetIndex1> getAlphabetIndex() const { return alphIndex_; }

  const AlphabetIndex1& alphabetIndex() const { return *alphIndex_; }

protected:
  virtual void alphabetIndexHasChanged() = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_REWARD_H
