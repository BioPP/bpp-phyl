//
// File: SubstitutionModelSet.h
// Created by: Bastien Boussau
//             Julien Dutheil
// Created on: Tue Aug 21 2007
//

/*
   Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef _SUBSTITUTIONMODELSET_H_
#define _SUBSTITUTIONMODELSET_H_


#include "../Tree.h"
#include "SubstitutionModel.h"
#include "AbstractSubstitutionModel.h"
#include "FrequenciesSet/FrequenciesSet.h"

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/VectorTools.h>

// From Seqlib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

// From the STL:
#include <vector>
#include <map>
#include <algorithm>
#include <memory>

namespace bpp
{
/**
 * @brief Substitution models manager for non-homogeneous / non-reversible models of evolution.
 *
 * This class contains a set of substitution models, and their assigment toward the branches of a phylogenetic tree.
 * Each branch in the tree corresponds to a model in the set, but a susbstitution model may correspond to several branches.
 * The particular case where all branches point toward a unique model is the homogeneous case.
 *
 * This class also deals with the parameters associated to the models.
 * In the homogeneous case, the parameter list is the same as the list in susbstitution model.
 * When two models at least are specified, these models may have their own parameters or share some of them.
 * To deal with this issue, the SubstitutionModelSet class contains its own parameter list and an index which tells to which
 * models these parameters apply to.
 * Since parameters in a list must have unique names, duplicated names are numbered according to the order in the list.
 * To track the relationships between names in the list and names in each model, the parameter list is duplicated in modelParameters_.
 * The user only act on parameters_, the fireParameterChanged function, automatically called, will update the modelParameters_ field.
 *
 * In the non-homogeneous and homogeneous non-reversible cases, the likelihood depends on the position of the root.
 * The states frequencies at the root of the tree are hence distinct parameters.
 * Theses are accounted by a FrequenciesSet objet, managed by the SubstitutionModelSet class.
 * The corresponding parameters, if any, are added at the begining of the global parameter list.
 *
 * If the heterogenity of the model does not affect the equilibrium frequencies, the model can be considered as stationary.
 * In such a model, the process is supposed to be at equilibrium all along the trees, including at the root.
 * Whether a model should be considered as stationary or not is left to the user. If the "asumme stationarity" option is set when
 * building the set, then no FrequenciesSet object is used, but the frequencies are taken to be the same as the one at the first
 * model in the set. Nothing hence prevents you to build a "supposingly stationary model which actually is not", so be careful!!
 *
 * This class provides several methods to specify which model and/or which parameter is associated to which branch/clade.
 * Several check points are provided, but some are probably missing due to the large set of possible models that this class allows to build,
 * so be carefull!
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
  const Alphabet* alphabet_;

  size_t nbStates_;

  /**
   * @brief Contains all models used in this tree.
   */
  std::vector<SubstitutionModel*> modelSet_;

private:
  /**
   * @brief Root frequencies.
   */
  std::auto_ptr<FrequenciesSet> rootFrequencies_;

  /**
   * @brief Contains for each node in a tree the index of the corresponding model in modelSet_
   */
  mutable std::map<int, size_t> nodeToModel_;
  mutable std::map<size_t, std::vector<int> > modelToNodes_;

  /**
   * @brief Contains for each parameter in the list the indexes of the corresponding models in modelSet_ that share this parameter.
   */
  std::vector< std::vector<size_t> > paramToModels_;

  std::map<std::string, size_t> paramNamesCount_;

  /**
   * @brief Contains for each parameter in the list the corresponding name in substitution models.
   */
  std::vector<std::string> modelParameterNames_;

  /**
   * @brief Parameters for each model in the set.
   *
   * The parameters_ field, inherited from AbstractSubstitutionModel contains all parameters, with unique names.
   * To make the correspondance with parameters for each model in the set, we duplicate them in this array.
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
  SubstitutionModelSet(const Alphabet* alpha):
    AbstractParameterAliasable(""),
    alphabet_(alpha),
    nbStates_(0),
    modelSet_(),
    rootFrequencies_(0),
    nodeToModel_(),
    modelToNodes_(),
    paramToModels_(),
    paramNamesCount_(),
    modelParameterNames_(),
    modelParameters_(),
    stationarity_(true)
  {
  }

  /**
   * @brief Create a model set according to the specified alphabet and root frequencies.
   * Stationarity is not assumed.
   *
   * @param alpha The alphabet to use for this set.
   * @param rootFreqs The frequencies at root node. The underlying object will be owned by this instance.
   */
  SubstitutionModelSet(const Alphabet* alpha, FrequenciesSet* rootFreqs):
    AbstractParameterAliasable(""),
    alphabet_(alpha),
    nbStates_(0),
    modelSet_(),
    rootFrequencies_(0),
    nodeToModel_(),
    modelToNodes_(),
    paramToModels_(),
    paramNamesCount_(),
    modelParameterNames_(),
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
  
  /**
   * @brief Sets a given FrequenciesSet for root frequencies.
   *
   * @param rootFreqs The FrequenciesSet for root frequencies.
   */

  void setRootFrequencies(FrequenciesSet* rootFreqs);

  SubstitutionModelSet(const SubstitutionModelSet& set);

  SubstitutionModelSet& operator=(const SubstitutionModelSet& set);

  virtual ~SubstitutionModelSet()
  {
    for (size_t i = 0; i < modelSet_.size(); i++) { delete modelSet_[i]; }
  }

#ifndef NO_VIRTUAL_COV
  SubstitutionModelSet*
#else
  Clonable*
#endif
  clone() const { return new SubstitutionModelSet(*this); }

public:
  /**
   * @brief Get the number of states associated to this model set.
   *
   * @return The number of states.
   * @throw Exception if no model is associated to the set.
   */
  size_t getNumberOfStates() const throw (Exception)
  {
    return nbStates_;
  }

  /**
   * @brief Get the index of a given parameter in the list of all parameters.
   *
   * @param name The name of the parameter to look for.
   * @return The position of the parameter in the global parameter list.
   * @throw ParameterNotFoundException If no parameter with this name is found.
   */
  size_t getParameterIndex(const std::string& name) const throw (ParameterNotFoundException)
  {
    for (size_t i = 0; i < getNumberOfParameters(); i++)
    {
      if (getParameter_(i).getName() == name) return i;
    }
    throw ParameterNotFoundException("SubstitutionModelSet::getParameterIndex().", name);
  }

  /**
   * @brief Get the model name of a given parameter in the list of all parameters.
   *
   * @param name The name of the parameter to look for.
   * @return The model name of the parameter in the global parameter list.
   * @throw ParameterNotFoundException If no parameter with this name is found.
   * @throw Exception If the parameter is not a 'model' parameter (that is, it is a root frequency parameter).
   */
  std::string getParameterModelName(const std::string& name) const throw (ParameterNotFoundException, Exception)
  {
    size_t pos = getParameterIndex(name);
    if (stationarity_)
      return modelParameterNames_[pos];
    else
    {
      size_t rfs = rootFrequencies_->getNumberOfParameters();
      if (pos < rfs) throw Exception("SubstitutionModelSet::getParameterModelName(). This parameter as no model name: " + name);
      return modelParameterNames_[pos - rfs];
    }
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
   * @return True iff there is a MixedSubstitutionModel in the SubstitutionModelSet
   **/

  bool hasMixedSubstitutionModel() const;

  /**
   * @brief Get one model from the set knowing its index.
   *
   * @param i Index of the model in the set.
   * @return A pointer toward the corresponding model.
   */
  const SubstitutionModel* getModel(size_t i) const throw (IndexOutOfBoundsException)
  {
    if (i > modelSet_.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::getNumberOfModels().", 0, modelSet_.size() - 1, i);
    return modelSet_[i];
  }

  SubstitutionModel* getModel(size_t i) throw (IndexOutOfBoundsException)
  {
    if (i > modelSet_.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::getNumberOfModels().", 0, modelSet_.size() - 1, i);
    return modelSet_[i];
  }

  /**
   * @brief Get the index in the set of the model associated to a particular node id.
   *
   * @param nodeId The id of the query node.
   * @return The index of the model associated to the given node.
   * @throw Exception If no model is found for this node.
   */
  size_t getModelIndexForNode(int nodeId) const throw (Exception)
  {
   std::map<int, size_t>::iterator i = nodeToModel_.find(nodeId);
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
  const SubstitutionModel* getModelForNode(int nodeId) const throw (Exception)
  {
   std::map<int, size_t>::const_iterator i = nodeToModel_.find(nodeId);
    if (i == nodeToModel_.end())
      throw Exception("SubstitutionModelSet::getModelForNode(). No model associated to node with id " + TextTools::toString(nodeId));
    return modelSet_[i->second];
  }
  SubstitutionModel* getModelForNode(int nodeId) throw (Exception)
  {
   std::map<int, size_t>::iterator i = nodeToModel_.find(nodeId);
    if (i == nodeToModel_.end())
      throw Exception("SubstitutionModelSet::getModelForNode(). No model associated to node with id " + TextTools::toString(nodeId));
    return modelSet_[i->second];
  }

  /**
   * @brief Get a list of nodes id for which the given model is associated.
   *
   * @param i The index of the model in the set.
   * @return A vector with the ids of the node associated to this model.
   * @throw IndexOutOfBoundsException If the index is not valid.
   */
  const std::vector<int>& getNodesWithModel(size_t i) const throw (IndexOutOfBoundsException)
  {
    if (i >= modelSet_.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::getNodesWithModel().", i, 0, modelSet_.size());
    return modelToNodes_[i];
  }

  /**
   * @param name The name of the parameter to look for.
   * @return The list of nodes with a model containing the specified parameter.
   * @throw ParameterNotFoundException If no parameter with the specified name is found.
   */
  std::vector<int> getNodesWithParameter(const std::string& name) const throw (ParameterNotFoundException);

  /**
   * @param name The name of the parameter to look for.
   * @return The list of model indices containing the specified parameter.
   * @throw ParameterNotFoundException If no parameter with the specified name is found.
   */
  std::vector<size_t> getModelsWithParameter(const std::string& name) const throw (ParameterNotFoundException);

  /**
   * @brief Add a new model to the set, and set relationships with nodes and params.
   *
   * @param model A pointer toward a susbstitution model, that will added to the set.
   * Warning! The set will now be the owner of the pointer, and will destroy it if needed!
   * Copy the model first if you don't want it to be lost!
   * @param nodesId the set of nodes in the tree that points toward this model.
   * This will override any previous affectation.
   * @param newParams The names of the parameters that have to be added to the global list.
   * These parameters will only be affected to this susbstitution model.
   * You can use the setParameterToModel function to assign this parameter to an additional model, and the
   * unsetParameterToModel to remove the relationship with this model for instance.
   * Parameters not specified in newParams will be ignored, unless you manually assign them to another parameter with
   * setParameterToModel.
   * @throw Exception in case of error:
   * <ul>
   * <li>if the new model does not match the alphabet<li>
   * <li>if the new model does not have the same number of states than existing ones<li>
   * <li>etc.</li>
   * </ul>
   */
  void addModel(SubstitutionModel* model, const std::vector<int>& nodesId, const std::vector<std::string>& newParams) throw (Exception);

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
  void setModel(SubstitutionModel* model, size_t modelIndex) throw (Exception, IndexOutOfBoundsException);

  /**
   * @brief Associate an existing model with a given node.
   *
   * If the node was was previously associated to a model, the old association is deleted.
   * If other nodes are associated to this model, the association is conserved.
   *
   * @param modelIndex The position of the model in the set.
   * @param nodeNumber The id of the corresponding node.
   */
  void setModelToNode(size_t modelIndex, int nodeNumber) throw (IndexOutOfBoundsException)
  {
    if (modelIndex >= nodeToModel_.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::setModelToNode.", modelIndex, 0, nodeToModel_.size() - 1);
    nodeToModel_[nodeNumber] = modelIndex;
  }

  /**
   * @brief Assign a parameter to a model.
   *
   * @param parameterIndex The index of the parameter in the list.
   * @param modelIndex     The index of the model in the list.
   * @throw IndexOutOfBoundsException If one of the index is not valid.
   */
  void setParameterToModel(size_t parameterIndex, size_t modelIndex) throw (IndexOutOfBoundsException);

  /**
   * @brief Unset a given parameter to the specified model.
   *
   * @param parameterIndex The index of the parameter in the list.
   * @param modelIndex     The index of the model in the list.
   * @throw IndexOutOfBoundsException If one of the index is not valid.
   * @throw Exception If the pseicified parameter is not currently associated to the specified model.
   */
  void unsetParameterToModel(size_t parameterIndex, size_t modelIndex) throw (IndexOutOfBoundsException, Exception);

  /**
   * @brief Add a parameter to the list, and link it to specified existing nodes.
   *
   * @param parameter The parameter to add. Its name must match model parameter names.
   * @param nodesId The list of ids of the nodes to link with these parameters.
   * Nodes must have a corresponding model in the set.
   * @throw Exception If one of the above requirement is not true.
   */
  void addParameter(const Parameter& parameter, const std::vector<int>& nodesId) throw (Exception);

  /**
   * @brief Add several parameters to the list, and link them to specified existing nodes.
   *
   * @param parameters The list of parameters to add. Its name must match model parameter names.
   * @param nodesId The list of ids of the nodes to link with these parameters.
   * Nodes must have a corresponding model in the set.
   * @throw Exception If one of the above requirement is not true.
   */
  void addParameters(const ParameterList& parameters, const std::vector<int>& nodesId) throw (Exception);

  /**
   * @brief Remove a parameter from the list, and unset it to all linked nodes and models.
   *
   * @param name The name of the parameter to remove.
   * @throw ParameterNotFoundException If no parameter with the given name is found in the list.
   */
  void removeParameter(const std::string& name) throw (ParameterNotFoundException);

  /**
   * @brief Remove a model from the set, and all corresponding parameters.
   *
   * @param modelIndex The index of the model in the set.
   * @throw Exception if a parameter becomes orphan because of the removal.
   */
  void removeModel(size_t modelIndex) throw (Exception);

  void listModelNames(std::ostream& out = std::cout) const;

  /**
   * @return The set of root frequencies.
   */
  const FrequenciesSet* getRootFrequenciesSet() const { return rootFrequencies_.get(); }

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

  const Alphabet* getAlphabet() const { return alphabet_; }

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
  bool isFullySetUpFor(const Tree& tree, bool throwEx = true) const throw (Exception)
  {
    return checkOrphanModels(throwEx)
           && checkOrphanParameters(throwEx)
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
  bool checkOrphanModels(bool throwEx) const throw (Exception);

  bool checkOrphanParameters(bool throwEx) const throw (Exception);

  bool checkOrphanNodes(const Tree& tree, bool throwEx) const throw (Exception);

  bool checkUnknownNodes(const Tree& tree, bool throwEx) const throw (Exception);
  /** @} */
};
} // end of namespace bpp.

#endif // _SUBSTITUTIONMODELSET_H_

