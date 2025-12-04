// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_SIMULATION_DETAILEDSITESIMULATOR_H
#define BPP_PHYL_SIMULATION_DETAILEDSITESIMULATOR_H


#include "../Likelihood/ParametrizablePhyloTree.h"
#include "MutationProcess.h"
#include "SiteSimulator.h"

// From the STL:
#include <map>
#include <vector>

namespace bpp
{
/**
 * @brief Data structure to store the result of a DetailedSiteSimulator.
 *
 * This data structure stores each transitional state, and the time when it occurred.
 */
class SiteSimulationResult
{
private:
  mutable std::map<unsigned int, size_t> indexes_;
  size_t currentIndex_;
  std::vector<MutationPath> paths_;
  std::vector<size_t> ancestralStates_;
  std::shared_ptr<const ParametrizablePhyloTree> tree_;
  std::vector<unsigned int> leavesId_;
  std::shared_ptr<const StateMapInterface> statemap_;

public:
  SiteSimulationResult(
      std::shared_ptr<const ParametrizablePhyloTree> tree,
      std::shared_ptr<const StateMapInterface> statemap,
      size_t ancestralState) :
    indexes_        (),
    currentIndex_   (0),
    paths_          (),
    ancestralStates_(),
    tree_           (tree),
    leavesId_       (tree->getNodeIndexes(tree->getLeavesUnderNode(tree->getRoot()))),
    statemap_       (statemap)
  {
    indexes_[tree->getRootIndex()] = 0;
    // Warning, watch out the indices there!
    paths_.push_back(MutationPath(getAlphabet(), ancestralState, 0));
    ancestralStates_.push_back(ancestralState);
  }

  SiteSimulationResult(const SiteSimulationResult& ssr) :
    indexes_        (ssr.indexes_),
    currentIndex_   (ssr.currentIndex_),
    paths_          (ssr.paths_),
    ancestralStates_(ssr.ancestralStates_),
    tree_           (ssr.tree_),
    leavesId_       (ssr.leavesId_),
    statemap_       (ssr.statemap_)
  {}

  SiteSimulationResult& operator=(const SiteSimulationResult& ssr)
  {
    indexes_         = ssr.indexes_;
    currentIndex_    = ssr.currentIndex_;
    paths_           = ssr.paths_;
    ancestralStates_ = ssr.ancestralStates_;
    tree_            = ssr.tree_;
    leavesId_        = ssr.leavesId_;
    statemap_        = ssr.statemap_;
    return *this;
  }

  virtual ~SiteSimulationResult() {}

public:
  /**
   * @return The alphabet associated to this simulation.
   */
  std::shared_ptr<const Alphabet> getAlphabet() const { return statemap_->getAlphabet(); }

  const Alphabet& alphabet() const { return statemap_->alphabet(); }

  virtual void addNode(unsigned int nodeId, MutationPath path)
  {
    currentIndex_++;
    indexes_[nodeId] = currentIndex_;
    paths_.push_back(path);
    ancestralStates_.push_back(path.getFinalState());
  }

  virtual size_t getAncestralState(size_t i) const { return ancestralStates_[i]; }

  virtual size_t getAncestralState(unsigned int nodeId) const { return ancestralStates_[indexes_[nodeId]]; }

  virtual const MutationPath& getMutationPath(size_t i) const { return paths_[i]; }

  virtual const MutationPath& getMutationPath(unsigned int nodeId) const { return paths_[indexes_[nodeId]]; }

  virtual size_t getSubstitutionCount(size_t i) const { return paths_[i].getNumberOfEvents(); }

  virtual void getSubstitutionCount(
      size_t i,
      const SubstitutionRegisterInterface& reg,
      std::vector<double>& counts) const
  {
    paths_[i].getEventCounts(counts, reg);
  }

  virtual size_t getSubstitutionCount(unsigned int nodeId) const { return paths_[indexes_[nodeId]].getNumberOfEvents(); }

  virtual void getSubstitutionCount(
      unsigned int nodeId,
      const SubstitutionRegisterInterface& reg,
      std::vector<double>& counts) const
  {
    paths_[indexes_[nodeId]].getEventCounts(counts, reg);
  }

  virtual VVdouble getSubstitutionVector(const SubstitutionRegisterInterface& reg) const
  {
    size_t n = paths_.size();
    VVdouble counts(n);
    for (size_t i = 0; i < n; ++i)
    {
      counts[i].resize(reg.getNumberOfSubstitutionTypes());
      paths_[i].getEventCounts(counts[i], reg);
    }
    return counts;
  }

  /**
   * @return The states at the leaves.
   */
  virtual std::vector<size_t> getFinalStates() const
  {
    size_t n = leavesId_.size();
    std::vector<size_t> states(n);
    for (size_t i = 0; i < n; i++)
    {
      states[i] = ancestralStates_[indexes_[leavesId_[i]]];
    }
    return states;
  }

  /**
   * @return The site corresponding to this simulation.
   */
  virtual std::unique_ptr<SiteInterface> getSite() const
  {
    std::vector<size_t> mstates = getFinalStates();
    std::vector<int> astates(mstates.size());
    for (size_t i = 0; i < mstates.size(); ++i)
    {
      astates[i] = statemap_->getAlphabetStateAsInt(mstates[i]);
    }

    auto alphabet = statemap_->getAlphabet();
    return std::make_unique<Site>(astates, alphabet);
  }

  /**
   * @return A vector with the leaves names.
   */
  virtual std::vector<std::string> getLeaveNames() const
  {
    size_t n = leavesId_.size();
    std::vector<std::string> names(n);
    for (size_t i = 0; i < n; i++)
    {
      names[i] = tree_->getNode(leavesId_[i])->getName();
    }
    return names;
  }
};

// ---------------------------------------------------------------------------

/**
 * @brief Data structure to store the result of a DetailedSiteSimulator.
 *
 * This structure inherits from the SequenceSimulationResult class, and add support for
 * rate variation across sites.
 */
class RASiteSimulationResult :
  public SiteSimulationResult
{
protected:
  double rate_;

public:
  RASiteSimulationResult(
      std::shared_ptr<const ParametrizablePhyloTree> tree,
      std::shared_ptr<const StateMapInterface> stateMap,
      size_t ancestralStateIndex,
      double rate) :
    SiteSimulationResult(tree, stateMap, ancestralStateIndex),
    rate_(rate) {}

  virtual ~RASiteSimulationResult() {}

public:
  /**
   * @return The rate of this simulation.
   */
  virtual double getRate() const { return rate_; }
};

// ---------------------------------------------------------------------------

/**
 * @brief This interface adds the dSimulate method to the SiteSimulator interface.
 *
 * Instances of this class should be used when a detailed output of the simulation is needed.
 */
class DetailedSiteSimulatorInterface :
  public virtual SiteSimulatorInterface
{
public:
  DetailedSiteSimulatorInterface() {}
  virtual ~DetailedSiteSimulatorInterface() {}

  DetailedSiteSimulatorInterface* clone() const override = 0;

public:
  /**
   * @brief Get a detailed simulation result for one site.
   *
   * @return A SiteSimulationResult object with all ancestral
   * states for all nodes and branches.
   */
  virtual std::unique_ptr<SiteSimulationResult> dSimulateSite() const = 0;
  virtual std::unique_ptr<SiteSimulationResult> dSimulateSite(size_t rateClass) const = 0;
  virtual std::unique_ptr<SiteSimulationResult> dSimulateSite(double rate) const = 0;
  virtual std::unique_ptr<SiteSimulationResult> dSimulateSite(size_t ancestralStateIndex, double rate) const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_SIMULATION_DETAILEDSITESIMULATOR_H
