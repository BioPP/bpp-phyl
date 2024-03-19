// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _LEGACY_MARGINAL_ANCESTRAL_STATES_RECONSTRUCTION_H_
#define _LEGACY_MARGINAL_ANCESTRAL_STATES_RECONSTRUCTION_H_

#include "../AncestralStateReconstruction.h"
#include "DRTreeLikelihood.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Sequence.h>

// From the STL:
#include <vector>

namespace bpp
{
/**
 * @brief Likelihood ancestral states reconstruction: marginal method.
 *
 * Reference:
 * Z Yang, S Kumar and M Nei (1995), _Genetics_ 141(4) 1641-50.
 */
class LegacyMarginalAncestralStateReconstruction :
  public virtual LegacyAncestralStateReconstruction
{
private:
  std::shared_ptr<const DRTreeLikelihoodInterface> likelihood_;
  TreeTemplate<Node> tree_;
  std::shared_ptr<const Alphabet> alphabet_;
  size_t nbSites_;
  size_t nbDistinctSites_;
  size_t nbClasses_;
  size_t nbStates_;
  std::vector<size_t> rootPatternLinks_;
  std::vector<double> r_;
  std::vector<double> l_;

public:
  LegacyMarginalAncestralStateReconstruction(std::shared_ptr<const DRTreeLikelihoodInterface> drl) :
    likelihood_      (drl),
    tree_            (drl->tree()),
    alphabet_        (drl->getAlphabet()),
    nbSites_         (drl->likelihoodData().getNumberOfSites()),
    nbDistinctSites_ (drl->likelihoodData().getNumberOfDistinctSites()),
    nbClasses_       (drl->likelihoodData().getNumberOfClasses()),
    nbStates_        (drl->likelihoodData().getNumberOfStates()),
    rootPatternLinks_(drl->likelihoodData().getRootArrayPositions()),
    r_               (drl->rateDistribution().getProbabilities()),
    l_               (drl->likelihoodData().getRootRateSiteLikelihoodArray())
  {}

  LegacyMarginalAncestralStateReconstruction(const LegacyMarginalAncestralStateReconstruction& masr) = default;

  LegacyMarginalAncestralStateReconstruction& operator=(const LegacyMarginalAncestralStateReconstruction& masr) = default;

  LegacyMarginalAncestralStateReconstruction* clone() const { return new LegacyMarginalAncestralStateReconstruction(*this); }

  virtual ~LegacyMarginalAncestralStateReconstruction() {}

public:
  /**
   * @brief Get ancestral states for a given node as a vector of int.
   *
   * The size of the vector is the number of distinct sites in the container
   * associated to the likelihood object.
   * This method is mainly for efficient internal use in other classes.
   * Consider using the getAncestralSequenceForNode() method for a more
   * general output.
   *
   * @param nodeId The id of the node at which the states must be reconstructed.
   * @param probs  A vector to be filled with the probability for each state at each position (will be the same size as the returned vector for states).
   * @param sample Tell if the sequence should be sample from the posterior distribution instead of taking the one with maximum probability.
   * @return A vector of states indices.
   * @see getAncestralSequenceForNode
   */
  std::vector<size_t> getAncestralStatesForNode(int nodeId, VVdouble& probs, bool sample) const;

  std::vector<size_t> getAncestralStatesForNode(int nodeId) const override
  {
    VVdouble probs(nbSites_);
    return getAncestralStatesForNode(nodeId, probs, false);
  }

  std::map<int, std::vector<size_t>> getAllAncestralStates() const override;

  /**
   * @brief Get the ancestral sequence for a given node.
   *
   * The name of the sequence will be the name of the node if there is one, its id otherwise.
   * A new sequence object is created, whose destruction is up to the user.
   *
   * @param nodeId The id of the node at which the sequence must be reconstructed.
   * @param probs  A pointer toward a vector to be filled with the probability for each state at each site (set to NULL if you don't want these probabilities).
   * @param sample Tell if the sequence should be sample from the posterior distribution instead of taking the one with maximum probability.
   * @return A sequence object.
   */
  std::unique_ptr<Sequence> getAncestralSequenceForNode(int nodeId, VVdouble* probs, bool sample) const;

  std::unique_ptr<Sequence> getAncestralSequenceForNode(int nodeId) const override
  {
    return getAncestralSequenceForNode(nodeId, 0, false);
  }

  std::unique_ptr<SiteContainerInterface> getAncestralSequences() const override
  {
    return getAncestralSequences(false);
  }

  std::unique_ptr<SiteContainerInterface> getAncestralSequences(bool sample) const;

private:
  void recursiveMarginalAncestralStates(
    const Node* node,
    std::map<int, std::vector<size_t>>& ancestors,
    AlignedSequenceContainer& data) const;
};
} // end of namespace bpp.

#endif // _LEGACY_MARGINAL_ANCESTRAL_STATES_RECONSTRUCTION_H_
