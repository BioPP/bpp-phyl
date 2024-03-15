// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_NONHOMOGENEOUSSUBSTITUTIONPROCESS_H
#define BPP_PHYL_LIKELIHOOD_NONHOMOGENEOUSSUBSTITUTIONPROCESS_H


#include "../Model/FrequencySet/FrequencySet.h"
#include "AbstractAutonomousSubstitutionProcess.h"

// From bpp-core:
#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>

// From the STL:
#include <vector>
#include <map>
// #include <algorithm>
#include <memory>
#include <numeric>

namespace bpp
{
/**
 * @brief Substitution process manager for non-homogeneous /
 * non-reversible models of evolution.
 *
 * This class contains a set of substitution models, and their
 * assignment toward the branches of a phylogenetic tree. Each branch
 * in the tree corresponds to a model in the set, but a susbstitution
 * model may correspond to several branches. The particular case where
 * all branches point toward a unique model is the homogeneous case.
 *
 * This class also deals with the parameters associated to the models.
 * The models may have their own parameters or share some of them. To
 * deal with this issue, the NonHomogeneousSubstitutionProcess class
 * contains its own parameter list and an index which tells to which
 * models these parameters apply to. Since parameters in a list must
 * have unique names, the names are suffixed with numbers according to
 * the order of the model in the list.
 *
 * To track the relationships between names in the list and names in
 * each model, the parameter list is duplicated in modelParameters_.
 * The user only act on parameters_, the fireParameterChanged
 * function, automatically called, will update the modelParameters_
 * field.
 *
 * In the non-homogeneous and homogeneous non-reversible cases, the
 * likelihood depends on the position of the root. The states
 * frequencies at the root of the tree are hence distinct parameters.
 * Theses are accounted by a FrequencySet objet, managed by the
 * NonHomogeneousSubstitutionProcess class. The corresponding
 * parameters, if any, are added at the begining of the global
 * parameter list.
 *
 * If the heterogenity of the model does not affect the equilibrium
 * frequencies, the model can be considered as stationary. In such a
 * model, the process is supposed to be at equilibrium all along the
 * trees, including at the root. Whether a model should be considered
 * as stationary or not is left to the user. If the "assume
 * stationarity" option is set when building the set, then no
 * FrequencySet object is used, but the frequencies are taken to be
 * the same as the one at the first model in the set. Nothing hence
 * prevents you to build a "supposingly stationary model which
 * actually is not", so be careful!!
 *
 * This class provides several methods to specify which model and/or
 * which parameter is associated to which branch/clade. Several check
 * points are provided, but some are probably missing due to the large
 * set of possible models that this class allows to build, so be
 * carefull!
 *
 */


class NonHomogeneousSubstitutionProcess :
  public AbstractAutonomousSubstitutionProcess
{
private:
  /**
   * @brief Contains all models used in this tree.
   */

  std::vector<std::shared_ptr<BranchModelInterface>> modelSet_;

  /**
   *  @brief Rate Distribution
   */
  std::shared_ptr<DiscreteDistributionInterface> rDist_;

  /**
   * @brief Contains for each node in a tree the index of the corresponding model in modelSet_
   */
  mutable std::map<unsigned int, size_t> nodeToModel_;
  mutable std::map<size_t, std::vector<unsigned int>> modelToNodes_;

  /**
   * @brief Parameters for each model in the set.
   */
  std::vector<ParameterList> modelParameters_;

public:

  /**
   * @brief Create a model set according to the specified alphabet and root frequencies.
   * Stationarity is not assumed.
   *
   * @param rdist  The DiscreteDistribution for the rates
   * @param tree the phylo tree tree
   * @param rootFreqs The frequencies at root node. The underlying object will be owned by this instance ( = 0 if stationary)
   */
  NonHomogeneousSubstitutionProcess(
      std::shared_ptr<DiscreteDistributionInterface> rdist,
      std::shared_ptr<const PhyloTree> tree = 0,
      std::shared_ptr<FrequencySetInterface> rootFreqs = 0) :
    AbstractParameterAliasable(""),
    AbstractAutonomousSubstitutionProcess(tree, rootFreqs),
    modelSet_(),
    rDist_(rdist),
    nodeToModel_(),
    modelToNodes_(),
    modelParameters_()
  {
    if (rDist_)
      addParameters_(rDist_->getIndependentParameters());
  }
  
  /**
   * @brief Create a model set according to the specified alphabet and root frequencies.
   * Stationarity is not assumed.
   *
   * @param rdist  The DiscreteDistribution for the rates
   * @param tree the parametrizable tree
   * @param rootFreqs The frequencies at root node. The underlying object will be owned by this instance ( = 0 if stationary)
   */
  NonHomogeneousSubstitutionProcess(
      std::shared_ptr<DiscreteDistributionInterface> rdist,
      std::shared_ptr<ParametrizablePhyloTree> tree,
      std::shared_ptr<FrequencySetInterface> rootFreqs = 0) :
    AbstractParameterAliasable(""),
    AbstractAutonomousSubstitutionProcess(tree, rootFreqs),
    modelSet_(),
    rDist_(rdist),
    nodeToModel_(),
    modelToNodes_(),
    modelParameters_()
  {
    if (rDist_)
      addParameters_(rDist_->getIndependentParameters());
  }

  NonHomogeneousSubstitutionProcess(const NonHomogeneousSubstitutionProcess& set);

  NonHomogeneousSubstitutionProcess& operator=(const NonHomogeneousSubstitutionProcess& set);

  virtual ~NonHomogeneousSubstitutionProcess()
  {
    clear();
  }

  NonHomogeneousSubstitutionProcess* clone() const override { return new NonHomogeneousSubstitutionProcess(*this); }

  /**
   * @brief Resets all the information contained in this object.
   *
   */

  void clear();

  /**
   * To be called when a parameter has changed.
   * Depending on parameters, this will actualize the rootFrequencies_
   * vector or the corresponding models in the set.
   *
   * @param parameters The modified parameters.
   */
  void fireParameterChanged(const ParameterList& parameters) override;

  const StateMapInterface& stateMap() const override
  {
    if (modelSet_.size() == 0)
      throw Exception("NonHomogeneousSubstitutionProcess::getStateMap : no model associated");
    else
      return modelSet_[0]->stateMap();
  }

  std::shared_ptr<const StateMapInterface> getStateMap() const override
  {
    if (modelSet_.size() == 0)
      throw Exception("NonHomogeneousSubstitutionProcess::getStateMap : no model associated");
    else
      return modelSet_[0]->getStateMap();
  }

  /**
   * @return The current number of distinct substitution models in this set.
   */
  size_t getNumberOfModels() const override { return modelSet_.size(); }

  /**
   * @return True iff there is a MixedTransitionModel in the NonHomogeneousSubstitutionProcess
   **/

  bool hasMixedTransitionModel() const;

  /**
   * @brief Set the modelPath, after checking  it is valid
   * (ie modelpath has only the model of the process).
   *
   */

  void setModelScenario(std::shared_ptr<ModelScenario> modelscenario) override;

  std::vector<size_t> getModelNumbers() const override
  {
    std::vector<size_t> v(getNumberOfModels());
    std::iota(std::begin(v), std::end(v), 1);
    return v;
  }

  const BranchModelInterface& model(size_t n) const override
  {
    if ((n == 0) || (n > modelSet_.size())) throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::model().", 1, modelSet_.size(), n);
    return *modelSet_[n - 1];
  }

  std::shared_ptr<const BranchModelInterface> getModel(size_t n) const override
  {
    if ((n == 0) || (n > modelSet_.size())) throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::getModel().", 1, modelSet_.size(), n);
    return modelSet_[n - 1];
  }

  std::shared_ptr<const BranchModelInterface> getModel(size_t n)
  {
    if ((n == 0) || (n > modelSet_.size())) throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::getModel().", 1, modelSet_.size(), n);
    return modelSet_[n - 1];
  }

  /**
   * @brief Get the number of the model associated to a
   * particular node id.
   *
   * @param nodeId The id of the query node.
   * @return The index of the model associated to the given node.
   * @throw Exception If no model is found for this node.
   */
  size_t getModelNumberForNode(unsigned int nodeId) const override
  {
    const auto i = nodeToModel_.find(nodeId);
    if (i == nodeToModel_.end())
      throw Exception("NonHomogeneousSubstitutionProcess::getModelNumberForNode(). No model associated to node with id " + TextTools::toString(nodeId));
    return i->second + 1;
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
      throw Exception("NonHomogeneousSubstitutionProcess::getModelForNode(). No model associated to node with id " + TextTools::toString(nodeId));
    return modelSet_[i->second];
  }

  /**
   * @brief Get a list of nodes id for which the given model is associated.
   *
   * @param i The number of the model in the set (ie the index in the set, plus 1).
   * @return A vector with the ids of the node associated to this model.
   * @throw IndexOutOfBoundsException If the index is not valid.
   */
  const std::vector<unsigned int> getNodesWithModel(size_t i) const override
  {
    if (i > modelSet_.size()) throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::getNodesWithModel().", i, 0, modelSet_.size());
    return modelToNodes_[i - 1];
  }

  /**
   * @brief Add a new model to the set, and set relationships with nodes and params.
   *
   * @param model A pointer toward a susbstitution model, that will added to the set.
   * Warning! The set will now be the owner of the pointer, and will destroy it if needed!
   * Copy the model first if you don't want it to be lost!
   * @param nodesId the set of nodes in the tree that points toward this model.
   * This will override any previous affectation.
   * @throw Exception in case of error:
   * <ul>
   * <li>if the new model does not match the alphabet<li>
   * <li>if the new model does not have the same number of states than existing ones<li>
   * <li>etc.</li>
   * </ul>
   */
  void addModel(std::shared_ptr<BranchModelInterface> model, const std::vector<unsigned int>& nodesId);

  /**
   * @brief Change a given model.
   *
   * The new model will be copied and will replace the old one.
   * All previous associations will be kept the same.
   * @param model A pointer toward a susbstitution model, that will added to the set.
   * Warning! The set will now be the owner of the pointer, and will destroy it if needed!
   * Copy the model first if you don't want it to be lost!
   * @param modelIndex The index of the existing model to replace.
   */
  void setModel(std::shared_ptr<BranchModelInterface> model, size_t modelIndex);

  /**
   * @brief Associate an existing model with a given node.
   *
   * If the node was was previously associated to a model, the old association is deleted.
   * If other nodes are associated to this model, the association is conserved.
   *
   * @param modelIndex The position of the model in the set.
   * @param nodeNumber The id of the corresponding node.
   */
  void setModelToNode(size_t modelIndex, unsigned int nodeNumber);

  /**
   * @brief list all model names.
   */
  void listModelNames(std::ostream& out = std::cout) const;

  ParameterList getBranchLengthParameters(bool independent) const override
  {
    if (getParametrizablePhyloTree())
      return getParametrizablePhyloTree()->getParameters();
    else
      return ParameterList();
  }

  /**
   * @brief Get the parameters attached to the rate distribution.
   */
  ParameterList getRateDistributionParameters(bool independent) const override
  {
    return rDist_.get() ? (independent ? rDist_->getIndependentParameters() : rDist_->getParameters()) : ParameterList();
  }

  const DiscreteDistributionInterface& rateDistribution() const override
  {
    if (!rDist_)
      throw NullPointerException("NonHomogeneousSubstitutionProcess::rateDistribution. No associated rate distribution.");
    return *rDist_;
  }

  DiscreteDistributionInterface& rateDistribution() override
  {
    if (!rDist_)
      throw NullPointerException("NonHomogeneousSubstitutionProcess::rateDistribution. No associated rate distribution.");
    return *rDist_;
  }

  std::shared_ptr<const DiscreteDistributionInterface> getRateDistribution() const override
  {
    return rDist_ ? rDist_ : nullptr;
  }

  std::shared_ptr<DiscreteDistributionInterface> getRateDistribution() override
  {
    return rDist_ ? rDist_ : nullptr;
  }

  /**
   * @brief Get the INDEPENDENT parameters corresponding to the models.
   */
  ParameterList getSubstitutionModelParameters(bool independent) const override;

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

public:
  /**
   * @return The values of the root frequencies.
   */
  const std::vector<double>& getRootFrequencies() const override
  {
    if (!hasRootFrequencySet())
    {
      if (std::dynamic_pointer_cast<const TransitionModelInterface>(modelSet_[0]))
        return std::dynamic_pointer_cast<const TransitionModelInterface>(modelSet_[0])->getFrequencies();
      else
        throw Exception("NonHomogeneousSubstitutionProcess::getRootFrequencies not callable.");
    }
    else
      return rootFrequencySet().getFrequencies();
  }

  const BranchModelInterface& model(unsigned int nodeId, size_t classIndex) const override
  {
    return *modelSet_[nodeToModel_[nodeId]];
  }

  std::shared_ptr<const BranchModelInterface> getModel(unsigned int nodeId, size_t classIndex) const override
  {
    return modelSet_[nodeToModel_[nodeId]];
  }

  double getProbabilityForModel(size_t classIndex) const override
  {
    if (classIndex >= (rDist_ ? rDist_->getNumberOfCategories() : 1))
      throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::getProbabilityForModel.", classIndex, 0, rDist_->getNumberOfCategories());
    return rDist_ ? rDist_->getProbability(classIndex) : 1.;
  }

  Vdouble getClassProbabilities() const override
  {
    Vdouble vProb;

    if (!rDist_)
      vProb.push_back(1.);
    else
      for (size_t i = 0; i < rDist_->getNumberOfCategories(); i++)
      {
        vProb.push_back(rDist_->getProbability(i));
      }

    return vProb;
  }

  double getRateForModel(size_t classIndex) const override
  {
    if (classIndex >= (rDist_ ? rDist_->getNumberOfCategories() : 1))
      throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::getRateForModel.", classIndex, 0, (rDist_ ? rDist_->getNumberOfCategories() : 1));
    return rDist_ ? rDist_->getCategory(classIndex) : 1;
  }

  /**
   * Static methods to create "simply" NonHomogeneousSubstitutionProcess.
   *
   */

  /**
   * @brief Create a NonHomogeneousSubstitutionProcess object,
   * corresponding to the homogeneous case.
   *
   * This class is mainly for testing purpose.
   *
   * @param model     The model to use.
   * @param rdist     The rate distribution
   * @param tree      The tree to use for the construction of the set.
   * @param rootFreqs A FrequencySet object to parametrize root frequencies
   *        (0 if stationary).
   * @param scenario (optional) the scenario used (in case of Mixed Models)
   */
  static std::unique_ptr<AutonomousSubstitutionProcessInterface> createHomogeneousSubstitutionProcess(
    std::shared_ptr<BranchModelInterface> model,
    std::shared_ptr<DiscreteDistributionInterface> rdist,
    std::shared_ptr<PhyloTree> tree,
    std::shared_ptr<FrequencySetInterface> rootFreqs = 0,
    std::shared_ptr<ModelScenario> scenario = 0
    );

  /**
   * @brief Create a NonHomogeneousSubstitutionProcess object, with one model per branch.
   *
   * All branches share the same type of model, but allow one set of parameters per branch.
   * This is also possible to specify some parameters to be common to all branches.
   *
   * @param model                The model to use.
   * @param rdist                The rate distribution
   * @param rootFreqs            A FrequencySet object to parametrize root frequencies.
   * @param tree                 The tree to use for the construction of the set.
   * @param globalParameterNames Common parameters for all branches.
   * All other parameters will be considered distinct for all branches.
   *
   * @param scenario (optional) the scenario used (in case of Mixed Models)
   */
  static std::unique_ptr<NonHomogeneousSubstitutionProcess> createNonHomogeneousSubstitutionProcess(
    std::shared_ptr<BranchModelInterface> model,
    std::shared_ptr<DiscreteDistributionInterface> rdist,
    std::shared_ptr<PhyloTree> tree,
    std::shared_ptr<FrequencySetInterface> rootFreqs,
    const std::vector<std::string>& globalParameterNames,
    std::shared_ptr<ModelScenario> scenario = 0
    );


  /** @} */
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_NONHOMOGENEOUSSUBSTITUTIONPROCESS_H
