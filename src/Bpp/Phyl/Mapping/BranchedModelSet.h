//
// File: BranchedModelSet.h
// Created by: Laurent Guéguen
// Created on: jeudi 7 décembre 2017, à 17h 41
//

/*
   Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _BRANCHED_MODELS_H_
#define _BRANCHED_MODELS_H_

#include <cstddef>
#include <vector>

namespace bpp
{
/**
 *
 * @brief An interface for classes where models are assigned to branch
 * ids.
 *
 */

class TransitionModel;

class BranchedModelSet
{
public:
  BranchedModelSet() {}
  virtual ~BranchedModelSet() {}

  /**
   * @return The current number of distinct substitution models in this set.
   *
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

  virtual const TransitionModel* getModel(size_t index) const = 0;

  /**
   * @brief Get the model associated to a particular branch id.
   *
   * @param branchId The id of the query branch.
   * @return A pointer toward the corresponding model.
   * @throw Exception If no model is found for this branch.
   */

  virtual const TransitionModel* getModelForBranch(unsigned int branchId) const = 0;

  virtual TransitionModel* getModelForBranch(unsigned int branchId) = 0;

  /**
   * @brief Get a list of branches id for which the given model is associated.
   *
   * @param index The index of the model.
   * @return A vector with the ids of the branch associated to this model.
   */

  virtual std::vector<unsigned int> getBranchesWithModel(size_t index) const = 0;
};
} // end of namespace bpp.

#endif// _BRANCHED_MODELS_H_
