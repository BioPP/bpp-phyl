//
// File: SubstitutionProcessCollectionMember.h
// Authors:
//   Laurent Guéguen
// Created: mercredi 13 mai 2015, ÃÂ  22h 32
//

/*
  Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_LIKELIHOOD_SUBSTITUTIONPROCESSCOLLECTIONMEMBER_H
#define BPP_PHYL_LIKELIHOOD_SUBSTITUTIONPROCESSCOLLECTIONMEMBER_H


#include "AbstractSubstitutionProcess.h"

namespace bpp
{
/**
 * @brief A substitution process which objects belong to a SubstitutionProcessCollection.
 *
 * The parameters are the INDEPENDENT parameters of the objects of
 * the Collection.
 */
class SubstitutionProcessCollection;

class SubstitutionProcessCollectionMember :
  public AbstractSubstitutionProcess
{
private:
  /**
   * @brief A pointer towards the collection the SubstitutionProcessCollectionMember
   * belongs to.
   */
  SubstitutionProcessCollection* pSubProColl_;

  /**
   * @brief The number of the process in the collection
   */
  size_t nProc_;

private:
  /**
   * @brief sets the parameters as the independent parameters on the
   * objects
   */
  void updateParameters();

private:
  /**
   * @brief Contains for each node in a tree the index of the corresponding model in modelSet_
   */

  std::map<unsigned int, size_t> nodeToModel_;
  std::map<size_t, std::vector<unsigned int> > modelToNodes_;

  /**
   * @brief The number of the tree: 0 means no assigned tree
   */
  size_t nTree_;

  /**
   * @brief The number of the rate distribution
   */
  size_t nDist_;

  /**
   * @brief The number of the root frequencies (0 if the process is stationary).
   */
  size_t nRoot_;

  /**
   * @brief the number of the set of model path, if needed.
   */
  size_t nPath_;

private:
  /*
   * @brief Constructors are only accessible through a SubstitutionProcessCollection.
   */

  /**
   * @brief Create a model set belonging to the specified SubstitutionProcessCollection.
   * Stationarity is assumed.
   *
   * @param pSubProColl the SubstitutionProcessCollection.
   * @param nProc Number of the process in the collection
   * @param nTree Number of the tree
   * @param nDist Number of the Discrete Distribution
   */
  SubstitutionProcessCollectionMember(
      SubstitutionProcessCollection* pSubProColl,
      size_t nProc,
      size_t nTree,
      size_t nDist) :
    AbstractParameterAliasable(""),
    pSubProColl_(pSubProColl),
    nProc_(nProc),
    nodeToModel_(),
    modelToNodes_(),
    nTree_(nTree),
    nDist_(nDist),
    nRoot_(0),
    nPath_(0)
  {
    updateParameters();
  }

  /**
   * @brief Resets all the information contained in this object.
   */
  SubstitutionProcessCollectionMember(const SubstitutionProcessCollectionMember& set) :
    AbstractParameterAliasable(set),
    pSubProColl_(set.pSubProColl_),
    nProc_(set.nProc_),
    nodeToModel_(set.nodeToModel_),
    modelToNodes_(set.modelToNodes_),
    nTree_(set.nTree_),
    nDist_(set.nDist_),
    nRoot_(set.nRoot_),
    nPath_(set.nPath_)
  {}

  SubstitutionProcessCollectionMember& operator=(const SubstitutionProcessCollectionMember& set)
  {
    AbstractParameterAliasable::operator=(set);
    pSubProColl_ = set.pSubProColl_;
    nProc_ = set.nProc_;

    nodeToModel_ = set.nodeToModel_;
    modelToNodes_ = set.modelToNodes_;
    nTree_ = set.nTree_;
    nDist_ = set.nDist_;
    nRoot_ = set.nRoot_;
    nPath_ = set.nPath_;

    return *this;
  }

  virtual ~SubstitutionProcessCollectionMember() {}
  
  struct Deleter
  {
    void operator()(SubstitutionProcessCollectionMember* sm) const
    {
      delete sm;
    }
  };

  SubstitutionProcessCollectionMember* clone() const override { return new SubstitutionProcessCollectionMember(*this); }

public:

  const SubstitutionProcessCollection& collection() const
  {
    if (! pSubProColl_)
      throw NullPointerException("SubstitutionProcessCollectionMember::collection(). Member is not associated to any collection.");
    return *pSubProColl_;
  }

  SubstitutionProcessCollection& collection()
  {
    if (! pSubProColl_)
      throw NullPointerException("SubstitutionProcessCollectionMember::collection(). Member is not associated to any collection.");
    return *pSubProColl_;
  }
  
  void clear()
  {
    nodeToModel_.clear();
    modelToNodes_.clear();
    nRoot_ = 0;
  }
  
  const SubstitutionProcessCollection* getCollection() const
  {
    return pSubProColl_;
  }

  SubstitutionProcessCollection* getCollection()
  {
    return pSubProColl_;
  }

  const StateMapInterface& stateMap() const override
  {
    return model(modelToNodes_.begin()->first).stateMap();
  }

  std::shared_ptr<const StateMapInterface> getStateMap() const override
  {
    return model(modelToNodes_.begin()->first).getStateMap();
  }

  /**
   * @return the number of the process in the collection
   */
  size_t getNumberOfProcesses() const
  {
    return nProc_;
  }

  /**
   * @return The current number of distinct substitution models in this set.
   */
  size_t getNumberOfModels() const override { return modelToNodes_.size(); }

  /**
   * @return True iff there is a MixedTransitionModel in the SubstitutionProcessCollectionMember
   */
  bool hasMixedTransitionModel() const;

  /**
   * @return True iff is stationary.
   */
  bool isStationary() const
  {
    return nRoot_ == 0;
  }

  const BranchModelInterface& model(size_t n) const override;
  
  std::shared_ptr<const BranchModelInterface> getModel(size_t n) const override;

  std::shared_ptr<BranchModelInterface> getModel(size_t n);

  std::vector<size_t> getModelNumbers() const override;

  /**
   * @brief Get the number of the model associated to a particular node id.
   *
   * @param nodeId The id of the query node.
   * @return The number of the model associated to the given node.
   * @throw Exception If no model is found for this node.
   */
  size_t getModelNumberForNode(unsigned int nodeId) const override
  {
    auto i = nodeToModel_.find(nodeId);
    if (i == nodeToModel_.end())
      throw Exception("SubstitutionProcessCollectionMember::getModelNumberForNode(). No model associated to node with id " + TextTools::toString(nodeId));
    return i->second;
  }

  /**
   * @brief Get the model associated to a particular node id.
   *
   * @param nodeId The id of the query node.
   * @return A pointer toward the corresponding model.
   * @throw Exception If no model is found for this node.
   */
  std::shared_ptr<const BranchModelInterface> getModelForNode(unsigned int nodeId) const override
  {
    std::map<unsigned int, size_t>::const_iterator i = nodeToModel_.find(nodeId);
    if (i == nodeToModel_.end())
      throw Exception("SubstitutionProcessCollectionMember::getModelForNode(). No model associated to node with id " + TextTools::toString(nodeId));
    return getModel(i->second);
  }

  /**
   * @brief Get a list of nodes id for which the given model is associated.
   *
   * @param i The index of the model in the set.
   * @return A vector with the ids of the node associated to this model.
   * @throw IndexOutOfBoundsException If the index is not valid.
   */
  const std::vector<unsigned int> getNodesWithModel(size_t i) const override
  {
    const auto it = modelToNodes_.find(i);
    if (it == modelToNodes_.end())
      throw Exception("SubstitutionProcessCollectionMember::getNodesWithModel(). No nodes associated with model " + TextTools::toString(i));

    return it->second;
  }

  /**
   * @brief Add a new model to the set, and set relationships with nodes.
   *
   * @param numModel The number of a model in the SubstitutionProcessCollection.
   * @param nodesId the set of nodes in the tree that points toward this model.
   * This will override any previous affectation.
   */
  void addModel(size_t numModel, const std::vector<unsigned int>& nodesId);

  /**
   * @brief Get the rate distribution
   */
  const DiscreteDistribution& rateDistribution() const override;

  DiscreteDistribution& rateDistribution() override;

  std::shared_ptr<const DiscreteDistribution> getRateDistribution() const override;

  std::shared_ptr<DiscreteDistribution> getRateDistribution() override;

  const size_t getRateDistributionNumber() const { return nDist_; }

  /**
   * @brief Set the root Frequencies Set
   * @param numFreq the index of the frequencies in the collection.
   *
   */
  void setRootFrequencies(size_t numFreq);

  const size_t getRootFrequenciesNumber() const { return nRoot_; }

  bool hasRootFrequencySet() const override { return !isStationary(); }

  const FrequencySetInterface& rootFrequencySet() const override;

  FrequencySetInterface& rootFrequencySet() override;

  std::shared_ptr<const FrequencySetInterface> getRootFrequencySet() const override;

  std::shared_ptr<FrequencySetInterface> getRootFrequencySet() override;

  /*
   * @brief Set the Set of Model Path
   * @param numPath the index of the frequencies in the collection.
   *
   * This method checks that the models defined in the process match
   * those of the model paths.
   */
  void setModelScenario(size_t numPath);

  const size_t getModelScenarioNumber() const { return nPath_; }

  std::shared_ptr<const ModelScenario> getModelScenario() const override;

  std::shared_ptr<ModelScenario> getModelScenario();

  /**
   * @brief AbsractParametrizable interface
   */
  bool matchParametersValues(const ParameterList& parameters) override;

  /**
   * @brief Check if the model set is fully specified for a given tree.
   *
   * This include:
   * - that each node as a model set up,
   * - that each model in the set is attributed to a node,
   * - all nodes ids in the set refer to an existing node in the tree.
   *
   * @param throwEx Tell if an exception have to be thrown in case of test not passed.
   */
  bool isFullySetUp(bool throwEx = true) const
  {
    return checkOrphanNodes(throwEx)
           && checkUnknownNodes(throwEx);
  }

protected:
  /**
   * @name Check function.
   *
   * @{
   */
  bool checkOrphanNodes(bool throwEx) const;

  bool checkUnknownNodes(bool throwEx) const;
  /** @} */

public:
  /**
   * @return The values of the root frequencies.
   */
  const std::vector<double>& getRootFrequencies() const override;

  /**
   * @return the Tree
   */
  const ParametrizablePhyloTree& parametrizablePhyloTree() const override;
  
  std::shared_ptr<const ParametrizablePhyloTree> getParametrizablePhyloTree() const override;

  ParametrizablePhyloTree& parametrizablePhyloTree();
  
  std::shared_ptr<ParametrizablePhyloTree> getParametrizablePhyloTree();

public:
  size_t getTreeNumber() const { return nTree_; }

  const BranchModelInterface& model(unsigned int nodeId, size_t classIndex) const override  {
    return model(nodeToModel_.at(nodeId));
  }

  std::shared_ptr<const BranchModelInterface> getModel(unsigned int nodeId, size_t classIndex) const override
  {
    return getModel(nodeToModel_.at(nodeId));
  }


  /**
   * @brief Get the parameters of the substitution models.
   */
  ParameterList getSubstitutionModelParameters(bool independent) const override;

  /**
   * @brief Get the parameters of the rate distribution.
   */
  ParameterList getRateDistributionParameters(bool independent) const override;

  /**
   * @brief Get the parameters of the tree.
   */
  ParameterList getBranchLengthParameters(bool independent) const override;

  /**
   * @brief Get the parameters of the root frequencies set.
   */
  ParameterList getRootFrequenciesParameters(bool independent) const override;

  /**
   * @brief get all NonDerivable parameters.
   */
  ParameterList getNonDerivableParameters() const override;

  double getProbabilityForModel(size_t classIndex) const override;

  Vdouble getClassProbabilities() const override;

  double getRateForModel(size_t classIndex) const override;

  friend class SubstitutionProcessCollection;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_SUBSTITUTIONPROCESSCOLLECTIONMEMBER_H
