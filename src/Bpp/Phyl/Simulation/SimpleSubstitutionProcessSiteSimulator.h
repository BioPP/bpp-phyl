//
// File: SimpleSubstitutionProcessSiteSimulator.h
// Authors:
//   Laurent GuÃ©guen
// Created: dimanche 24 mai 2020, Ã  07h 30
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_SIMULATION_SIMPLESUBSTITUTIONPROCESSSITESIMULATOR_H
#define BPP_PHYL_SIMULATION_SIMPLESUBSTITUTIONPROCESSSITESIMULATOR_H

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Phyl/NewLikelihood/ProcessComputationTree.h>

#include "../Model/SubstitutionModel.h"
#include "../NewLikelihood/ParametrizablePhyloTree.h"
#include "DetailedSiteSimulator.h"

// From SeqLib:
#include <Bpp/Seq/Site.h>

// From the STL:
#include <map>
#include <vector>

#include "../NewLikelihood/SubstitutionProcess.h"

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
  std::vector<std::shared_ptr<SimProcessNode> > sons_;

public:
  SimProcessNode(const ProcessComputationNode& pcn) :
    ProcessComputationNode(pcn), state_() {}

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
 *
 */

class SimpleSubstitutionProcessSiteSimulator :
  public DetailedSiteSimulator
{
protected:
  const SubstitutionProcess*     process_;
  const ParametrizablePhyloTree* phyloTree_;


  /*
   * @brief To store states & transition probabilities of the simulator
   *
   */

  SPTree tree_;

  /*
   *@brief cumsum probas of the substitution rates
   *
   */

  Vdouble qRates_;

  /*
   *@brief cumsum probas of the root frequencies, one per rate class
   *(this "per class" is useful for
   * GivenDataSubstitutionProcessSiteSimulator)
   *
   */

  VVdouble qRoots_;

  /**
   * @brief Vector of indexes of sequenced output species
   *
   */

  std::vector<size_t> seqIndexes_;

  /**
   * @brief Vector of names of sequenced output species
   *
   */

  std::vector<std::string> seqNames_;

  /**
   * @brief Map between species Indexes & used nodes, may change at
   * each simulation.
   *
   */

  mutable std::map<size_t, std::shared_ptr<SimProcessNode> > speciesNodes_;

  size_t nbNodes_;
  size_t nbClasses_;
  size_t nbStates_;

  bool continuousRates_;

  // Should we ouptut internal sequences as well?
  bool outputInternalSites_;

  /**
   * @name Stores intermediate results.
   *
   * @{
   */

public:
  SimpleSubstitutionProcessSiteSimulator(
    const SubstitutionProcess& process);

  virtual ~SimpleSubstitutionProcessSiteSimulator()
  {}

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

  SimpleSubstitutionProcessSiteSimulator* clone() const { return new SimpleSubstitutionProcessSiteSimulator(*this); }

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
  Site* simulateSite() const;

  Site* simulateSite(size_t rateClass) const;

  Site* simulateSite(double rate) const;

  Site* simulateSite(size_t ancestralStateIndex, double rate) const;

  std::vector<std::string> getSequencesNames() const { return seqNames_; }
  /** @} */

  /**
   * @name SiteSimulator interface
   *
   * @{
   */
  const Alphabet* getAlphabet() const { return process_->getStateMap().getAlphabet(); }
  /** @} */

  /**
   * @name The DetailedSiteSimulator interface.
   *
   * @{
   */

  SiteSimulationResult* dSimulateSite() const;

  SiteSimulationResult* dSimulateSite(size_t rateClass) const;

  SiteSimulationResult* dSimulateSite(double rate) const;

  SiteSimulationResult* dSimulateSite(size_t ancestralStateIndex, double rate) const;

  /** @} */

  /**
   * @brief Get the substitution process associated to this instance.
   *
   * @return The substitution process associated to this instance.
   */
  const SubstitutionProcess* getSubstitutionProcess() const { return process_; }

  /**
   * @brief Get the tree associated to this instance.
   *
   * @return The Tree object associated to this instance.
   */
  const ParametrizablePhyloTree* getTree() const { return phyloTree_; }

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
  void outputInternalSites(bool yn);

protected:
  /**
   * @name The 'Internal' methods.
   *
   * @{
   */

  /**
   * This method uses the states_ variable for saving ancestral states.
   */
  void evolveInternal(std::shared_ptr<SimProcessNode> node, size_t rateClass, SiteSimulationResult* ssr = 0) const;

  /**
   * This method uses the states_ variable for saving ancestral states.
   */
  void evolveInternal(std::shared_ptr<SimProcessNode> node, double rate, SiteSimulationResult* ssr = 0) const;

  /** @} */
};
} // end of namespace bpp.
#endif // BPP_PHYL_SIMULATION_SIMPLESUBSTITUTIONPROCESSSITESIMULATOR_H
