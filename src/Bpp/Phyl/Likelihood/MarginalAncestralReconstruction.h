// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_MARGINALANCESTRALRECONSTRUCTION_H
#define BPP_PHYL_LIKELIHOOD_MARGINALANCESTRALRECONSTRUCTION_H


#include "../AncestralStateReconstruction.h"
#include "DataFlow/DataFlowCWise.h"
#include "DataFlow/LikelihoodCalculationSingleProcess.h"

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
class MarginalAncestralReconstruction :
  public virtual AncestralStateReconstruction
{
private:
  std::shared_ptr<LikelihoodCalculationSingleProcess> likelihood_;
  std::shared_ptr<const ParametrizablePhyloTree> tree_;
  mutable std::shared_ptr<const Alphabet> alphabet_;
  size_t nbSites_;
  size_t nbDistinctSites_;
  // size_t nbClasses_;
  size_t nbStates_;
  PatternType rootPatternLinks_;

public:
  MarginalAncestralReconstruction(std::shared_ptr<LikelihoodCalculationSingleProcess> drl) :
    likelihood_      (drl),
    tree_            (drl->substitutionProcess().getParametrizablePhyloTree()),
    alphabet_        (drl->stateMap().getAlphabet()),
    nbSites_         (drl->getNumberOfSites()),
    nbDistinctSites_ (drl->getNumberOfDistinctSites()),
    nbStates_        (drl->stateMap().getNumberOfModelStates()),
    rootPatternLinks_(drl->getRootArrayPositions())
  {
    if (!tree_)
      throw Exception("MarginalAncestralReconstruction::MarginalAncestralReconstruction: missing ParametrizablePhyloTree.");
  }

  MarginalAncestralReconstruction(const MarginalAncestralReconstruction& masr) :
    likelihood_      (masr.likelihood_),
    tree_            (masr.tree_),
    alphabet_        (masr.alphabet_),
    nbSites_         (masr.nbSites_),
    nbDistinctSites_ (masr.nbDistinctSites_),
    nbStates_        (masr.nbStates_),
    rootPatternLinks_(masr.rootPatternLinks_)
  {}

  MarginalAncestralReconstruction& operator=(const MarginalAncestralReconstruction& masr)
  {
    likelihood_       = masr.likelihood_;
    tree_             = masr.tree_;
    alphabet_         = masr.alphabet_;
    nbSites_          = masr.nbSites_;
    nbDistinctSites_  = masr.nbDistinctSites_;
    nbStates_         = masr.nbStates_;
    rootPatternLinks_ = masr.rootPatternLinks_;
    return *this;
  }


  MarginalAncestralReconstruction* clone() const { return new MarginalAncestralReconstruction(*this); }

  virtual ~MarginalAncestralReconstruction() {}

public:
  std::shared_ptr<const Alphabet> getAlphabet() const
  {
    return alphabet_;
  }

  /**
   * @brief Get ancestral states  for a given node as a vector of int.
   *
   * The size of the vector is the total number of sites in the container
   * associated to the likelihood object.
   * This method is mainly for efficient internal use in other classes.
   * Consider using the getAncestralSequenceForNode() method for a more
   * general output.
   *
   * @param nodeId The id of the node at which the states must be
   * reconstructed [in].
   * @param probs  A vector to be filled with the probability for
   * each state at each position (will be the same size as the
   * returned vector for states) [out].
   * @param sample Tell if the sequence should be sampled from the
   * posterior distribution instead of taking the one with maximum
   * probability.
   * @return A vector of states indices.
   * @see getAncestralSequenceForNode
   */
  std::vector<size_t> getAncestralStatesForNode(unsigned int nodeId, VVdouble& probs, bool sample) const;

  std::vector<size_t> getAncestralStatesForNode(unsigned int nodeId) const override
  {
    VVdouble probs(nbSites_);
    return getAncestralStatesForNode(nodeId, probs, false);
  }

  std::map<unsigned int, std::vector<size_t>> getAllAncestralStates() const override;

  /**
   * @brief Get an ancestral sequence for a given node.
   *
   * The name of the sequence will be the name of the node if
   * there is one, its id otherwise. A new sequence object is
   * created, whose destruction is up to the user.
   *
   * @param nodeId The id of the node at which the sequence must
   * be reconstructed.
   * @param probs A pointer toward a vector to be filled with the
   * probability for each state at each site (set to NULL if you
   * don't want these probabilities).
   * @param sample Tell if the sequence should be sample from the
   * posterior distribution instead of taking the one with maximum
   * probability.
   * @return A sequence object.
   */
  std::unique_ptr<Sequence> getAncestralSequenceForNode(unsigned int nodeId, VVdouble* probs, bool sample) const;

  std::unique_ptr<Sequence> getAncestralSequenceForNode(unsigned int nodeId) const override
  {
    return getAncestralSequenceForNode(nodeId, 0, false);
  }

  std::unique_ptr<AlignedSequenceContainer> getAncestralSequences() const override
  {
    return getAncestralSequences(false);
  }

  std::unique_ptr<AlignedSequenceContainer> getAncestralSequences(bool sample) const;

private:
  void recursiveMarginalAncestralStates(
      const std::shared_ptr<PhyloNode> node,
      std::map<unsigned int, std::vector<size_t>>& ancestors,
      AlignmentDataInterface& data) const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_MARGINALANCESTRALRECONSTRUCTION_H
