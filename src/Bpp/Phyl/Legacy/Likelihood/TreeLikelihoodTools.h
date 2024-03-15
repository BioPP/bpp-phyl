// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_TREELIKELIHOODTOOLS_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_TREELIKELIHOODTOOLS_H


#include "TreeLikelihood.h"

// From the STL:
#include <vector>
#include <map>

namespace bpp
{
/**
 * @brief Utilitary methods that work with TreeLikelihood objects.
 */
class TreeLikelihoodTools
{
public:
  TreeLikelihoodTools() {}
  virtual ~TreeLikelihoodTools() {}

public:
  /**
   * @brief Compute the expected ancestral frequencies of all states at all (inner) nodes
   * according to a Markov process defined by a given substitution model.
   *
   * The computation is performed for a given site. If the likelihood object has no
   * site partition, then the method will return the same result for all positions.
   *
   * @param tl          [in] A tree likelihood object.
   * @param site        [in] The site for which the frequencies should be computed.
   * @param frequencies [out] A map where to store the results, as a vector of double (the
   * size of which being equal to the number of states in the model), and with nodes id as keys.
   * @param alsoForLeaves [opt] Tell if frequencies should also be estimated for terminal nodes.
   * @throw Exception In case something bad happens, like an unvalid model set.
   */
  static void getAncestralFrequencies(
    const TreeLikelihoodInterface& tl,
    size_t site,
    std::map<int, std::vector<double> >& frequencies,
    bool alsoForLeaves = false);

  /**
   * @brief Compute the expected ancestral frequencies of all states at all (inner) nodes
   * according to a Markov process defined by a given substitution model.
   *
   * The computation is averaged over all sites. If the likelihood object has no
   * site partition, then the method will return the same result as all single site numbers.
   *
   * @param tl          [in] A tree likelihood object.
   * @param frequencies [out] A map where to store the results, as a vector of double (the
   * size of which being equal to the number of states in the model), and with nodes id as keys.
   * @param alsoForLeaves [opt] Tell if frequencies should also be estimated for terminal nodes.
   * @throw Exception In case something bad happens, like an unvalid model set.
   */
  static void getAncestralFrequencies(
    const TreeLikelihoodInterface& tl,
    std::map<int, std::vector<double> >& frequencies,
    bool alsoForLeaves = false);

private:
  /**
   * @brief Recursive method, for internal use only.
   *
   * @see getAncestralFrequencies()
   */
  static void getAncestralFrequencies_(
    const TreeLikelihoodInterface& tl,
    size_t siteIndex,
    int parentId,
    const std::vector<double>& ancestralFrequencies,
    std::map<int, std::vector<double> >& frequencies,
    bool alsoForLeaves);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_TREELIKELIHOODTOOLS_H
