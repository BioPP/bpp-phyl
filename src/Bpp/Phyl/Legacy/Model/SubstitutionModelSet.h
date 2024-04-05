// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_MODEL_SUBSTITUTIONMODELSET_H
#define BPP_PHYL_LEGACY_MODEL_SUBSTITUTIONMODELSET_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/VectorTools.h>

#include "../../Model/AbstractSubstitutionModel.h"
#include "../../Model/FrequencySet/FrequencySet.h"
#include "../../Tree/Tree.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

// From the STL:
#include <vector>
#include <map>
#include <algorithm>
#include <memory>
#include <typeinfo>

namespace bpp
{
/**
 * @brief Substitution models manager for non-homogeneous / non-reversible models of evolution.
 *
 * This class contains a set of substitution models, and their
 * assignment toward the branches of a phylogenetic tree. Each branch
 * in the tree corresponds to a model in the set, but a susbstitution
 * model may correspond to several branches. The particular case where
 * all branches point toward a unique model is the homogeneous case.
 *
 * This class also deals with the parameters associated to the models.
 * In the homogeneous case, the parameter list is the same as the list
 * in susbstitution model. When two models at least are specified,
 * these models may have their own parameters or share some of them.
 * To deal with this issue, the SubstitutionModelSet class contains
 * its own parameter list, where parameters are numbered according to
 * the model they belong to.
 *
 * The user only act on parameters_, the fireParameterChanged
 * function, automatically called, will update the modelParameters_
 * field.
 *
 * In the non-homogeneous and homogeneous non-reversible cases, the likelihood depends on the position of the root.
 * The states frequencies at the root of the tree are hence distinct parameters.
 * These are accounted by a FrequencySet object, managed by the SubstitutionModelSet class.
 * The corresponding parameters, if any, are added at the beginning of the global parameter list.
 *
 * If the heterogenity of the model does not affect the equilibrium frequencies, the model can be considered as stationary.
 * In such a model, the process is supposed to be at equilibrium all along the trees, including at the root.
 * Whether a model should be considered as stationary or not is left to the user. If the "asumme stationarity" option is set when
 * building the set, then no FrequencySet object is used, but the frequencies are taken to be the same as the one at the first
 * model in the set. Nothing hence prevents you to build a "supposedly stationary model which actually is not", so be careful!!
 *
 * This class provides several methods to specify which model and/or which parameter is associated to which branch/clade.
 * Several check points are provided, but some are probably missing due to the large set of possible models that this class allows to build,
 * so be careful!
 *
 * @see SubstitutionModelSetTools for methods that provide instances of the SubstitutionModelSet for general cases.
 */
class SubstitutionModelSet :
  public AbstractParameterAliasable
{
protected:
  /**
   * @brief A pointer toward the common alphabet to all models in the set.
   */
  std::shared_ptr<const Alphabet> alphabet_;

  size_t nbStates_;

  /**
   * @brief Contains all models used in this tree.
   */
  std::vector< std::shared_ptr<TransitionModelInterface>> modelSet_;

private:
  /**
   * @brief Root frequencies.
   */
  std::shared_ptr<FrequencySetInterface> rootFrequencies_;

  /**
   * @brief Contains for each node in a tree the index of the corresponding model in modelSet_
   */
  mutable std::map<int, size_t> nodeToModel_;
  mutable std::map<size_t, std::vector<int>> modelToNodes_;

  /**
   * @brief Parameters for each model in the set.
   *
   * The parameters_ field, inherited from AbstractSubstitutionModel contains all parameters, with unique names.
   * To make the correspondence with parameters for each model in the set, we duplicate them in this array.
   * In most cases, this is something like 'theta_1 <=> theta', 'theta_2 <=> theta', etc.
   */
  std::vector<ParameterList> modelParameters_;

  bool stationarity_;

public:
  /**
   * @brief Create a model set according to the specified alphabet.
   * Stationarity is assumed.
   *
   * @param alpha The alphabet to use for this set.
   */
  SubstitutionModelSet(std::shared_ptr<const Alphabet> alpha) :
    AbstractParameterAliasable(""),
    alphabet_(alpha),
    nbStates_(0),
    modelSet_(),
    rootFrequencies_(),
    nodeToModel_(),
    modelToNodes_(),
    modelParameters_(),
    stationarity_(true)
  {}

  /**
   * @brief Create a model set according to the specified alphabet and root frequencies.
   * Stationarity is not assumed.
   *
   * @param alpha The alphabet to use for this set.
   * @param rootFreqs The frequencies at root node. The underlying object will be owned by this instance.
   */
  SubstitutionModelSet(
      std::shared_ptr<const Alphabet> alpha,
      std::shared_ptr<FrequencySetInterface> rootFreqs) :
    AbstractParameterAliasable(""),
    alphabet_(alpha),
    nbStates_(0),
    modelSet_(),
    rootFrequencies_(),
    nodeToModel_(),
    modelToNodes_(),
    modelParameters_(),
    stationarity_(true)
  {
    setRootFrequencies(rootFreqs);
  }

  /**
   * @brief Resets all the information contained in this object.
   *
   */
  void clear();

  bool isStationary() const
  {
    return stationarity_;
  }

  /**
   * @brief Sets a given FrequencySet for root frequencies.
   *
   * @param rootFreqs The FrequencySet for root frequencies.
   */
  void setRootFrequencies(std::shared_ptr<FrequencySetInterface> rootFreqs);

  SubstitutionModelSet(const SubstitutionModelSet& set);

  SubstitutionModelSet& operator=(const SubstitutionModelSet& set);

  virtual ~SubstitutionModelSet() {}

  SubstitutionModelSet* clone() const { return new SubstitutionModelSet(*this); }

public:
  /**
   * @brief Get the number of states associated to this model set.
   *
   * @return The number of states.
   * @throw Exception if no model is associated to the set.
   */
  size_t getNumberOfStates() const
  {
    return nbStates_;
  }

  /**
   * To be called when a parameter has changed.
   * Depending on parameters, this will actualize the _initialFrequencies vector or the corresponding models in the set.
   * @param parameters The modified parameters.
   */
  virtual void fireParameterChanged(const ParameterList& parameters);

  /**
   * @return The current number of distinct substitution models in this set.
   */
  size_t getNumberOfModels() const { return modelSet_.size(); }

  /**
   * @return True iff there is a MixedTransitionModel in the SubstitutionModelSet
   **/

  bool hasMixedTransitionModel() const;

  /**
   * @brief Get one model from the set knowing its index.
   *
   * @param i Index of the model in the set.
   * @return A pointer toward the corresponding model.
   */
  std::shared_ptr<const TransitionModelInterface> getModel(size_t i) const
  {
    if (i >= modelSet_.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::getModel(i).", i, 0, modelSet_.size() - 1);
    return modelSet_[i];
  }

  const TransitionModelInterface& model(size_t i) const
  {
    if (i >= modelSet_.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::model(i).", i, 0, modelSet_.size() - 1);
    return *modelSet_[i];
  }

  std::shared_ptr<TransitionModelInterface> getModel(size_t i)
  {
    if (i >= modelSet_.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::getModel(i).", i, 0, modelSet_.size() - 1);
    return modelSet_[i];
  }

  TransitionModelInterface& model(size_t i)
  {
    if (i >= modelSet_.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::model(i).", i, 0, modelSet_.size() - 1);
    return *modelSet_[i];
  }
  /**
   * @brief Return a markovian substitution model (or null)
   */
  std::shared_ptr<const SubstitutionModelInterface> getSubstitutionModel(size_t i) const
  {
    auto m = std::dynamic_pointer_cast<const SubstitutionModelInterface>(getModel(i));
    if (m) return m;
    else
      throw Exception("SubstitutionModelSet::getSubstitutionModel : " + model(i).getName() + " is  not a substitution model.");
  }


  std::shared_ptr<SubstitutionModelInterface> getSubstitutionModel(size_t i)
  {
    auto m = std::dynamic_pointer_cast<SubstitutionModelInterface>(getModel(i));
    if (m) return m;
    else
      throw Exception("SubstitutionModelSet::getSubstitutionModel : " + model(i).getName() + " is  not a substitution model.");
  }

  const SubstitutionModelInterface& substitutionModel(size_t i) const
  {
    try
    {
      auto& m = dynamic_cast<const SubstitutionModelInterface&>(model(i));
      return m;
    }
    catch (Exception& ex)
    {
      throw Exception("SubstitutionModelSet::substitutionModel : " + model(i).getName() + " is not a substitution model.");
    }
  }

  /**
   * @brief check if has only markovian substitution models
   */
  bool hasOnlySubstitutionModels() const
  {
    for (const auto& mod : modelSet_)
    {
      auto model = std::dynamic_pointer_cast<const SubstitutionModelInterface>(mod);
      if (!model) return false;
    }
    return true;
  }

  /**
   * @brief Get the index in the set of the model associated to a particular node id.
   *
   * @param nodeId The id of the query node.
   * @return The index of the model associated to the given node.
   * @throw Exception If no model is found for this node.
   */
  size_t getModelIndexForNode(int nodeId) const
  {
    auto i = nodeToModel_.find(nodeId);
    if (i == nodeToModel_.end())
      throw Exception("SubstitutionModelSet::getModelIndexForNode(). No model associated to node with id " + TextTools::toString(nodeId));
    return i->second;
  }

  /**
   * @brief Get the model associated to a particular node id.
   *
   * @param nodeId The id of the query node.
   * @return A pointer toward the corresponding model.
   * @throw Exception If no model is found for this node.
   */
  std::shared_ptr<const TransitionModelInterface> getModelForNode(int nodeId) const
  {
    auto i = nodeToModel_.find(nodeId);
    if (i == nodeToModel_.end())
      throw Exception("SubstitutionModelSet::getModelForNode(). No model associated to node with id " + TextTools::toString(nodeId));
    return modelSet_[i->second];
  }
  std::shared_ptr<TransitionModelInterface> getModelForNode(int nodeId)
  {
    auto i = nodeToModel_.find(nodeId);
    if (i == nodeToModel_.end())
      throw Exception("SubstitutionModelSet::getModelForNode(). No model associated to node with id " + TextTools::toString(nodeId));
    return modelSet_[i->second];
  }

  std::shared_ptr<const SubstitutionModelInterface> getSubstitutionModelForNode(int nodeId) const
  {
    return std::dynamic_pointer_cast<const SubstitutionModelInterface>(getModelForNode(nodeId));
  }

  std::shared_ptr<SubstitutionModelInterface> getSubstitutionModelForNode(int nodeId)
  {
    return std::dynamic_pointer_cast<SubstitutionModelInterface>(getModelForNode(nodeId));
  }

  /**
   * @brief Get a list of nodes id for which the given model is associated.
   *
   * @param i The index of the model in the set.
   * @return A vector with the ids of the node associated to this model.
   * @throw IndexOutOfBoundsException If the index is not valid.
   */
  const std::vector<int>& getNodesWithModel(size_t i) const
  {
    if (i >= modelSet_.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::getNodesWithModel().", i, 0, modelSet_.size());
    return modelToNodes_[i];
  }

  /**
   * @param name The name of the parameter to look for.
   * @return The list of nodes with a model containing the specified parameter.
   * @throw ParameterNotFoundException If no parameter with the specified name is found.
   */

  std::vector<int> getNodesWithParameter(const std::string& name) const;

  /**
   * @brief Add a new model to the set, and set relationships with nodes and params.
   *
   * @param model A pointer toward a susbstitution model, that will added to the set.
   * Warning! The set will now be the owner of the pointer, and will destroy it if needed!
   * Copy the model first if you don't want it to be lost!
   * @param nodesId the set of nodes in the tree that points toward this model.
   * This will override any previous affectation.
   *
   * @throw Exception in case of error:
   * <ul>
   * <li>if the new model does not match the alphabet<li>
   * <li>if the new model does not have the same number of states than existing ones<li>
   * <li>etc.</li>
   * </ul>
   */
  void addModel(std::shared_ptr<TransitionModelInterface> model, const std::vector<int>& nodesId); // , const std::vector<std::string>& newParams);

  /**
   * @brief Sets an assignment of a given model index to a given onde id
   *
   * @param modelIndex The index of the model in the set.
   * @param nodeId      The node ID
   *
   * @throw Exception if the model index doesn't correspond to an existing model in the modelSet
   */
  void setNodeToModel(size_t modelIndex, int nodeId); // Keren: added on my own to allow alternation of nodes assignments to existing nodes

  /**
   * @brief Reset model indices to node ids assignment
   */
  void resetModelToNodeIds();  // Keren: added on my own to allow alternation of nodes assignments to existing nodes

  /**
   * @brief Replace a model in the set, and all corresponding
   * parameters. The replaced model deleted.
   *
   * @param modelIndex The index of the model in the set.
   * @param model the new model. This model will be owned by the Set.
   *
   * @throw Exception if a parameter becomes orphan because of the removal.
   */
  void replaceModel(size_t modelIndex, std::shared_ptr<TransitionModelInterface> model);

  void listModelNames(std::ostream& out = std::cout) const;

  /**
   * @return The set of root frequencies.
   */
  const std::shared_ptr<FrequencySetInterface> getRootFrequencySet() const { return rootFrequencies_; }

  /**
   * @return The values of the root frequencies.
   */
  std::vector<double> getRootFrequencies() const
  {
    if (stationarity_)
      return modelSet_[0]->getFrequencies();
    else
      return rootFrequencies_->getFrequencies();
  }

  /**
   * @brief Get the parameters corresponding to the root frequencies.
   *
   * @return The parameters corresponding to the root frequencies.
   */
  ParameterList getRootFrequenciesParameters() const
  {
    if (stationarity_)
      return ParameterList();
    else
      return rootFrequencies_->getParameters();
  }

  /**
   * @brief Get the parameters corresponding attached to the nodes of the tree.
   *
   * That is, all the parameters without the root frequencies.
   *
   * @return The parameters attached to the nodes of the tree.
   */
  ParameterList getNodeParameters() const
  {
    ParameterList pl;
    for (size_t i = stationarity_ ? 0 : rootFrequencies_->getNumberOfParameters();
        i < getNumberOfParameters(); i++)
    {
      pl.addParameter(getParameter_(i));
    }
    return pl;
  }

  /**
   * @brief Get the parameters attached to a Model.
   *
   * @param modelIndex the index of the model in the set
   *
   * @return The parameters attached to the model.
   */

  ParameterList getModelParameters(size_t modelIndex) const;

  std::shared_ptr<const Alphabet> getAlphabet() const { return alphabet_; }

  const Alphabet& alphabet() const { return *alphabet_; }

  /**
   * @return The supported states of the model set, as a vector of int codes.
   *
   * @see Alphabet
   */
  const std::vector<int>& getAlphabetStates() const
  {
    return model(0).getAlphabetStates();
  }

  const StateMapInterface& stateMap() const
  {
    return model(0).stateMap();
  }

  std::shared_ptr<const StateMapInterface> getStateMap() const
  {
    return model(0).getStateMap();
  }

  std::vector<size_t> getModelStates(int code) const
  {
    return model(0).getModelStates(code);
  }

  std::vector<size_t> getModelStates(const std::string& code) const
  {
    return model(0).getModelStates(code);
  }

  /**
   * @param index The model state.
   * @return The corresponding alphabet state as character code.
   */
  int getAlphabetStateAsInt(size_t index) const
  {
    return model(0).getAlphabetStateAsInt(index);
  }

  /**
   * @param index The model state.
   * @return The corresponding alphabet state as character code.
   */
  std::string getAlphabetStateAsChar(size_t index) const
  {
    return model(0).getAlphabetStateAsChar(index);
  }

  /**
   * @brief Check if the model set is fully specified for a given tree.
   *
   * This include:
   * - that each node as a model set up,
   * - that each model in the set is attributed to a node,
   * - that each parameter in the set actually correspond to a model.
   * - all nodes ids in the set refer to an existing node in the tree.
   *
   * @param tree The tree to check.
   * @param throwEx Tell if an exception have to be thrown in case of test not passed.
   */
  bool isFullySetUpFor(const Tree& tree, bool throwEx = true) const
  {
    return checkOrphanModels(throwEx)
           //           && checkOrphanParameters(throwEx)
           && checkOrphanNodes(tree, throwEx)
           && checkUnknownNodes(tree, throwEx);
  }

protected:
  /**
   * Set rootFrequencies_ from parameters.
   */
  void updateRootFrequencies()
  {
    if (!stationarity_)
      rootFrequencies_->matchParametersValues(getParameters());
  }

  /**
   * @name Check function.
   *
   * @{
   */
  bool checkOrphanModels(bool throwEx) const;

  bool checkOrphanNodes(const Tree& tree, bool throwEx) const;

  bool checkUnknownNodes(const Tree& tree, bool throwEx) const;
  /** @} */
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_MODEL_SUBSTITUTIONMODELSET_H
