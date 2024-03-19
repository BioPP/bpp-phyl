// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_MODEL_SUBSTITUTIONMODELSETTOOLS_H
#define BPP_PHYL_LEGACY_MODEL_SUBSTITUTIONMODELSETTOOLS_H


#include "../../Tree/Tree.h"
#include "SubstitutionModelSet.h"

// From the STL:
#include <vector>
#include <map>

namespace bpp
{
/**
 * @brief Tools for automatically creating SubstitutionModelSet objects.
 */
class SubstitutionModelSetTools
{
public:
  /**
   * @brief Create a SubstitutionModelSet object, corresponding to the homogeneous case.
   *
   * This class is mainly for testing purpose.
   *
   * @param model     The model to use.
   * @param rootFreqs A FrequencySet object to parametrize root frequencies.
   * @param tree      The tree to use for the construction of the set.
   */
  static std::unique_ptr<SubstitutionModelSet> createHomogeneousModelSet(
      std::shared_ptr<TransitionModelInterface> model,
      std::shared_ptr<FrequencySetInterface> rootFreqs,
      const Tree& tree
      );

  /**
   * @brief Create a SubstitutionModelSet object, with one model per branch.
   *
   * All branches share the same type of model, but allow one set of parameters per branch.
   * This is also possible to specify some parameters to be common to all branches.
   *
   * @param model                The model to use.
   * @param rootFreqs            A FrequencySet object to parametrize root frequencies.
   * @param tree                 The tree to use for the construction of the set.
   * @param aliasFreqNames       Aliases for the frequencies names
   * @param globalParameterNames map for shared parameters.
   *                             Associated value is a vector of
   *                             vectors of branch ids that share
   *                             the parameter. If the vector is
   *                             empty, the parameter is shared
   *                             among all branches.
   *
   * All other parameters will be considered distinct for all branches.
   */
  static std::unique_ptr<SubstitutionModelSet> createNonHomogeneousModelSet(
      std::shared_ptr<TransitionModelInterface> model,
      std::shared_ptr<FrequencySetInterface> rootFreqs,
      const Tree& tree,
      const std::map<std::string, std::string>& aliasFreqNames,
      const std::map<std::string, std::vector<Vint>>& globalParameterNames
      );
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_MODEL_SUBSTITUTIONMODELSETTOOLS_H
