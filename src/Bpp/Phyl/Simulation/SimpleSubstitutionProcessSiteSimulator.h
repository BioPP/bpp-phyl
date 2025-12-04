// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_SIMULATION_SIMPLESUBSTITUTIONPROCESSSITESIMULATOR_H
#define BPP_PHYL_SIMULATION_SIMPLESUBSTITUTIONPROCESSSITESIMULATOR_H

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Phyl/Likelihood/ProcessComputationTree.h>

#include "../Likelihood/ParametrizablePhyloTree.h"
#include "../Model/SubstitutionModel.h"
#include "DetailedSiteSimulator.h"

// From SeqLib:
#include <Bpp/Seq/Site.h>

// From the STL:
#include <map>
#include <vector>

#include "../Likelihood/SubstitutionProcess.h"

namespace bpp
{
class SimProcessNode :
  public ProcessComputationNode
{
private:
  // states during simulation
  size_t state_;

  // probabilities to choose, in case of mixture node or root, for all rates
  VVdouble cumProb_;

  // Sons in case of mixture node
  std::vector<std::shared_ptr<SimProcessNode>> sons_;

public:
  SimProcessNode(const ProcessComputationNode& pcn) :
    ProcessComputationNode(pcn), state_(), cumProb_(), sons_() {}

  friend class SimpleSubstitutionProcessSiteSimulator;
  friend class GivenDataSubstitutionProcessSiteSimulator;
};

class SimProcessEdge :
  public ProcessComputationEdge
{
private:
  // Cumulative pxy for all rates
  VVVdouble cumpxy_;

public:
  SimProcessEdge(const ProcessComputationEdge& pce) :
    ProcessComputationEdge(pce), cumpxy_() {}

  friend class SimpleSubstitutionProcessSiteSimulator;
  friend class GivenDataSubstitutionProcessSiteSimulator;
};

typedef AssociationTreeGlobalGraphObserver<SimProcessNode, SimProcessEdge>  SPTree;


/**
 * @brief Site simulation under a unique substitution process.
 */
class SimpleSubstitutionProcessSiteSimulator :
  public virtual DetailedSiteSimulatorInterface
{
protected:
  std::shared_ptr<const SubstitutionProcessInterface> process_;
  std::shared_ptr<const ParametrizablePhyloTree> phyloTree_;


  /**
   * @brief To store states & transition probabilities of the simulator
   */
  SPTree tree_;

  /**
   * @brief cumsum probas of the substitution rates
   */
  Vdouble qRates_;

  /*
   * @brief cumsum probas of the root frequencies, one per rate class
   * (this "per class" is useful for GivenDataSubstitutionProcessSiteSimulator)
   */
  VVdouble qRoots_;

  /**
   * @brief Vector of indexes of sequenced output species
   */
  std::vector<size_t> seqIndexes_;

  /**
   * @brief Vector of names of sequenced output species
   */
  std::vector<std::string> seqNames_;

  /**
   * @brief Map between species Indexes & used nodes, may change at
   * each simulation.
   */
  mutable std::map<size_t, std::shared_ptr<SimProcessNode>> speciesNodes_;

  size_t nbNodes_;
  size_t nbClasses_;
  size_t nbStates_;

  bool continuousRates_;

  // Should we output internal sequences as well?
  bool outputInternalSites_;

  /**
   * @name Stores intermediate results.
   *
   * @{
   */

public:
  SimpleSubstitutionProcessSiteSimulator(
      std::shared_ptr<const SubstitutionProcessInterface> process);

  virtual ~SimpleSubstitutionProcessSiteSimulator() {}

  SimpleSubstitutionProcessSiteSimulator(const SimpleSubstitutionProcessSiteSimulator& nhss) :
    process_        (nhss.process_),
    phyloTree_      (nhss.phyloTree_),
    tree_           (nhss.tree_),
    qRates_         (nhss.qRates_),
    qRoots_         (nhss.qRoots_),
    seqIndexes_     (nhss.seqIndexes_),
    seqNames_       (nhss.seqNames_),
    speciesNodes_  (nhss.speciesNodes_),
    nbNodes_        (nhss.nbNodes_),
    nbClasses_      (nhss.nbClasses_),
    nbStates_       (nhss.nbStates_),
    continuousRates_(nhss.continuousRates_),
    outputInternalSites_(nhss.outputInternalSites_)
  {}

  SimpleSubstitutionProcessSiteSimulator& operator=(const SimpleSubstitutionProcessSiteSimulator& nhss)
  {
    process_        = nhss.process_;
    phyloTree_       = nhss.phyloTree_;
    tree_            = nhss.tree_;
    qRates_          = nhss.qRates_;
    qRoots_          = nhss.qRoots_;
    seqIndexes_      = nhss.seqIndexes_;
    seqNames_        = nhss.seqNames_;
    speciesNodes_   = nhss.speciesNodes_;
    nbNodes_         = nhss.nbNodes_;
    nbClasses_       = nhss.nbClasses_;
    nbStates_        = nhss.nbStates_;
    continuousRates_ = nhss.continuousRates_;
    outputInternalSites_ = nhss.outputInternalSites_;

    return *this;
  }

  SimpleSubstitutionProcessSiteSimulator* clone() const override
  {
    return new SimpleSubstitutionProcessSiteSimulator(*this);
  }

private:
  /**
   * @brief Init all probabilities.
   *
   * Method called by constructors.
   */
  virtual void init();

public:
  /**
   * @name The SiteSimulator interface
   *
   * @{
   */

  std::unique_ptr<Site> simulateSite() const override;

  std::unique_ptr<Site> simulateSite(size_t rateClass) const override;

  std::unique_ptr<Site> simulateSite(double rate) const override;

  std::unique_ptr<Site> simulateSite(size_t ancestralStateIndex, double rate) const override;

  std::vector<std::string> getSequenceNames() const override
  {
    return seqNames_;
  }

  /** @} */

  /**
   * @name SiteSimulator interface
   *
   * @{
   */
  std::shared_ptr<const Alphabet> getAlphabet() const override { return process_->stateMap().getAlphabet(); }

  const Alphabet& alphabet() const override { return process_->stateMap().alphabet(); }
  /** @} */

  /**
   * @name The DetailedSiteSimulator interface.
   *
   * @{
   */

  std::unique_ptr<SiteSimulationResult> dSimulateSite() const override;

  std::unique_ptr<SiteSimulationResult> dSimulateSite(size_t rateClass) const override;

  std::unique_ptr<SiteSimulationResult> dSimulateSite(double rate) const override;

  std::unique_ptr<SiteSimulationResult> dSimulateSite(size_t ancestralStateIndex, double rate) const override;

  /** @} */

  /**
   * @brief Get the substitution process associated to this instance.
   *
   * @return The substitution process associated to this instance.
   */
  std::shared_ptr<const SubstitutionProcessInterface> getSubstitutionProcess() const
  {
    return process_;
  }

  /**
   * @brief Get the tree associated to this instance.
   *
   * @return The Tree object associated to this instance.
   */
  std::shared_ptr<const ParametrizablePhyloTree> getTree() const
  {
    return phyloTree_;
  }

  /**
   * @brief Enable the use of continuous rates instead of discrete rates.
   *
   * To work, the DiscreteDistribution object used should implement
   * the randC method.
   *
   * In this case, sampling is done jointly on the process classes
   * and on the continuous rate.
   *
   * @param yn Tell if we should use continuous rates.
   */
  void enableContinuousRates(bool yn) { continuousRates_ = yn; }

  /**
   * @brief Sets whether we will output the internal sequences or not.
   *
   *
   * @param yn Tell if we should output internal sequences.
   */
  void outputInternalSites(bool yn) override;

protected:
  /**
   * @name The 'Internal' methods.
   *
   * @{
   */

  /**
   * This method uses the states_ variable for saving ancestral states.
   */
  void evolveInternal(
      std::shared_ptr<SimProcessNode> node,
      size_t rateClass, SiteSimulationResult* ssr = nullptr) const;

  /**
   * This method uses the states_ variable for saving ancestral states.
   */
  void evolveInternal(
      std::shared_ptr<SimProcessNode> node,
      double rate,
      SiteSimulationResult* ssr = nullptr) const;

  /** @} */
};
} // end of namespace bpp.
#endif // BPP_PHYL_SIMULATION_SIMPLESUBSTITUTIONPROCESSSITESIMULATOR_H
