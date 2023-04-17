//
// File: SubstitutionModelSetTools.h
// Authors:
//   Julien Dutheil
// Created: 2007-09-17 16:57:00
//

/*
  Copyright or ÃÂ© or Copr. CNRS, (November 16, 2004)
  
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
    const std::map<std::string, std::vector<Vint> >& globalParameterNames
    );
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_MODEL_SUBSTITUTIONMODELSETTOOLS_H
