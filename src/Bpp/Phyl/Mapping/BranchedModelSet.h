// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_BRANCHEDMODELSET_H
#define BPP_PHYL_MAPPING_BRANCHEDMODELSET_H

#include <cstddef>
#include <vector>


namespace bpp
{
/**
 * @brief An interface for classes where models are assigned to branch
 * ids.
 */

class TransitionModelInterface;

class BranchedModelSet
{
public:
  BranchedModelSet() {}
  virtual ~BranchedModelSet() {}

  /**
   * @return The current number of distinct substitution models in this set.
   */
  virtual size_t getNumberOfModels() const = 0;

  /**
   * @return The vector of model indexes.
   *
   */

  virtual std::vector<size_t> getModelNumbers() const = 0;

  /**
   * @brief Get the model with a ginev index.
   *
   * @param index The index of the query model.
   * @return A pointer toward the corresponding model.
   */
  virtual std::shared_ptr<const TransitionModelInterface> getModel(size_t index) const = 0;

  /**
   * @brief Get the model associated to a particular branch id.
   *
   * @param branchId The id of the query branch.
   * @return A pointer toward the corresponding model.
   * @throw Exception If no model is found for this branch.
   */
  virtual std::shared_ptr<const TransitionModelInterface> getModelForBranch(unsigned int branchId) const = 0;

  virtual std::shared_ptr<TransitionModelInterface> getModelForBranch(unsigned int branchId) = 0;

  /**
   * @brief Get a list of branches id for which the given model is associated.
   *
   * @param index The index of the model.
   * @return A vector with the ids of the branch associated to this model.
   */

  virtual std::vector<unsigned int> getBranchesWithModel(size_t index) const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_BRANCHEDMODELSET_H
