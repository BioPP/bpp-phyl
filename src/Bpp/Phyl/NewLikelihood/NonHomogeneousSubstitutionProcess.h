//
// File: NonHomogeneousSubstitutionProcess.h
// Created by: Julien Dutheil, Bastien Boussau, Laurent Guéguen
// Created on: jeudi 20 juin 2013, à 23h 08
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

#ifndef _NONHOMOGENEOUSSUBSTITUTIONPROCESS_H_
#define _NONHOMOGENEOUSSUBSTITUTIONPROCESS_H_

#include "AbstractSubstitutionProcess.h"
#include "../Model/FrequenciesSet/FrequenciesSet.h"

//From bpp-core:
#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>

// From the STL:
#include <vector>
#include <map>
//#include <algorithm>
#include <memory>

namespace bpp
{
  /**
   * @brief Substitution process manager for non-homogeneous /
   * non-reversible models of evolution.
   *
   * This class contains a set of substitution models, and their
   * assigment toward the branches of a phylogenetic tree. Each branch
   * in the tree corresponds to a model in the set, but a susbstitution
   * model may correspond to several branches. The particular case where
   * all branches point toward a unique model is the homogeneous case.
   *
   * This class also deals with the parameters associated to the models.
   * The models may have their own parameters or share some of them. To
   * deal with this issue, the NonHomogeneousSubstitutionProcess class
   * contains its own parameter list and an index which tells to which
   * models these parameters apply to. Since parameters in a list must
   * have unique names, the names are numbered according to the order of
   * the model in the list.
   * To track the relationships between names in the list and names in
   * each model, the parameter list is duplicated in modelParameters_.
   * The user only act on parameters_, the fireParameterChanged
   * function, automatically called, will update the modelParameters_
   * field.
   *
   * In the non-homogeneous and homogeneous non-reversible cases, the
   * likelihood depends on the position of the root. The states
   * frequencies at the root of the tree are hence distinct parameters.
   * Theses are accounted by a FrequenciesSet objet, managed by the
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
   * FrequenciesSet object is used, but the frequencies are taken to be
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
    public AbstractSubstitutionProcess
  {
  private:
    /**
     * @brief Contains all models used in this tree.
     * Note that auto_ptr objects cannot be stored in a vector.
     */
    
    std::vector<SubstitutionModel* > modelSet_;

    /**
     * @brief Root frequencies.
     */
    std::auto_ptr<FrequenciesSet> rootFrequencies_;

    /**
     *  @brief Rate Distribution
     */

    std::auto_ptr<DiscreteDistribution> rDist_;

    /**
     * @brief Contains for each node in a tree the index of the corresponding model in modelSet_
     */
    mutable std::map<int, size_t> nodeToModel_;
    mutable std::map<size_t, std::vector<int> > modelToNodes_;

    /**
     * @brief Parameters for each model in the set.
     *
     */

    std::vector<ParameterList> modelParameters_;

    bool stationarity_;

    /**
     * @brief The related Computing Tree
     *
     */

    mutable std::auto_ptr<ComputingTree> computingTree_;


  public:
    /**
     * @brief Create an empty model set with stationarity assumed.
     *
     * @param rdist  The DiscreteDistribution for the rates
     * @param tree the parametrizable tree
     */
    NonHomogeneousSubstitutionProcess(DiscreteDistribution*  rdist, ParametrizableTree* tree) :
      AbstractParameterAliasable(""),
      AbstractSubstitutionProcess(tree, rdist ? rdist->getNumberOfCategories() : 0),
      modelSet_(),
      rootFrequencies_(0),
      rDist_(rdist),
      nodeToModel_(),
      modelToNodes_(),
      modelParameters_(),
      stationarity_(true),
      computingTree_()
    {
      // Add parameters:
      addParameters_(tree->getParameters());  //Branch lengths
      addParameters_(rdist->getIndependentParameters());

      computingTree_.reset(new ComputingTree(*pTree_.get(), *rDist_.get()));
    }

    /**
     * @brief Create a model set according to the specified alphabet and root frequencies.
     * Stationarity is not assumed.
     *
     * @param rdist  The DiscreteDistribution for the rates
     * @param tree the parametrizable tree
     * @param rootFreqs The frequencies at root node. The underlying object will be owned by this instance.
     */
    NonHomogeneousSubstitutionProcess(DiscreteDistribution*  rdist, ParametrizableTree* tree, FrequenciesSet* rootFreqs):
      AbstractParameterAliasable(""),
      AbstractSubstitutionProcess(tree, rdist ? rdist->getNumberOfCategories() : 0),
      modelSet_(),
      rootFrequencies_(rootFreqs),
      rDist_(rdist),
      nodeToModel_(),
      modelToNodes_(),
      modelParameters_(),
      stationarity_(false),
      computingTree_()
    {
      addParameters_(tree->getParameters());  //Branch lengths
      addParameters_(rdist->getIndependentParameters());  
      setRootFrequencies(rootFreqs);
      
      computingTree_.reset(new ComputingTree(*pTree_.get(), *rDist_.get()));
    }

    NonHomogeneousSubstitutionProcess(const NonHomogeneousSubstitutionProcess& set);

    NonHomogeneousSubstitutionProcess& operator=(const NonHomogeneousSubstitutionProcess& set);

    virtual ~NonHomogeneousSubstitutionProcess()  {
      clear();
    }

    NonHomogeneousSubstitutionProcess* clone() const { return new NonHomogeneousSubstitutionProcess(*this); }

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
    void fireParameterChanged(const ParameterList& parameters);

    /**
     * @brief Get the alphabet
     * @throw Exception if no model is associated to the set.
     *
     */

    const Alphabet* getAlphabet() const
    {
      if (modelSet_.size()==0)
        throw Exception("NonHomogeneousSubstitutionProcess::getAlphabet : no model associated");
      else
        return modelSet_[0]->getAlphabet();
    }

    /**
     * @return The current number of distinct substitution models in this set.
     */

    size_t getNumberOfModels() const { return modelSet_.size(); }

    /**
     * @return True iff there is a MixedSubstitutionModel in the NonHomogeneousSubstitutionProcess
     **/

    bool hasMixedSubstitutionModel() const;

    /**
     * @brief Get one model from the set knowing its index.
     *
     * @param i Index of the model in the set.
     * @return A pointer toward the corresponding model.
     */

    const SubstitutionModel* getModel(size_t i) const
    {
      if (i > modelSet_.size()) throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::getNumberOfModels().", 0, modelSet_.size() - 1, i);
      return modelSet_[i];
    }

    SubstitutionModel* getModel(size_t i)
    {
      if (i > modelSet_.size()) throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::getNumberOfModels().", 0, modelSet_.size() - 1, i);
      return modelSet_[i];
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
      std::map<int, size_t>::iterator i = nodeToModel_.find(nodeId);
      if (i == nodeToModel_.end())
        throw Exception("NonHomogeneousSubstitutionProcess::getModelIndexForNode(). No model associated to node with id " + TextTools::toString(nodeId));
      return i->second;
    }

    /**
     * @brief Get the model associated to a particular node id.
     *
     * @param nodeId The id of the query node.
     * @return A pointer toward the corresponding model.
     * @throw Exception If no model is found for this node.
     */
  
    const SubstitutionModel* getModelForNode(int nodeId) const
    {
      std::map<int, size_t>::const_iterator i = nodeToModel_.find(nodeId);
      if (i == nodeToModel_.end())
        throw Exception("NonHomogeneousSubstitutionProcess::getModelForNode(). No model associated to node with id " + TextTools::toString(nodeId));
      return modelSet_[i->second];
    }
  
    SubstitutionModel* getModelForNode(int nodeId)
    {
      std::map<int, size_t>::iterator i = nodeToModel_.find(nodeId);
      if (i == nodeToModel_.end())
        throw Exception("NonHomogeneousSubstitutionProcess::getModelForNode(). No model associated to node with id " + TextTools::toString(nodeId));
      return modelSet_[i->second];
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
      if (i >= modelSet_.size()) throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::getNodesWithModel().", i, 0, modelSet_.size());
      return modelToNodes_[i];
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

    void addModel(SubstitutionModel* model, const std::vector<int>& nodesId);

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

    void setModel(SubstitutionModel* model, size_t modelIndex);

    /**
     * @brief Associate an existing model with a given node.
     *
     * If the node was was previously associated to a model, the old association is deleted.
     * If other nodes are associated to this model, the association is conserved.
     *
     * @param modelIndex The position of the model in the set.
     * @param nodeNumber The id of the corresponding node.
     */

    void setModelToNode(size_t modelIndex, int nodeNumber);
 
    /**
     * @brief list all model names.
     *
     */
  
    void listModelNames(std::ostream& out = std::cout) const;

    /**
     * @brief Sets a given FrequenciesSet for root frequencies.
     *
     * @param rootFreqs The FrequenciesSet for root frequencies.
     */

    void setRootFrequencies(FrequenciesSet* rootFreqs);

    /**
     * @return The set of root frequencies.
     *
     */
  
    const FrequenciesSet* getRootFrequenciesSet() const
    {
      return rootFrequencies_.get();
    }

    /**
     * @brief Get the parameters corresponding to the root frequencies.
     *
     * @return The parameters corresponding to the root frequencies.
     */

    ParameterList getRootFrequenciesParameters(bool independent) const
    {
      if (stationarity_)
        return ParameterList();
      else
        return rootFrequencies_->getParameters();
    }

    ParameterList getBranchLengthParameters(bool independent) const
    {
      return getParametrizableTree().getParameters();
    }

    /**
     * @brief Get the parameters attached to the rate distribution.
     *
     */
     
    ParameterList getRateDistributionParameters(bool independent) const
    {
      return (rDist_.get()?(independent?rDist_->getIndependentParameters():rDist_->getParameters()):ParameterList());
    }

    const DiscreteDistribution* getRateDistribution() const
    {
      return rDist_.get();
    }

    /**
     * @brief Get the INDEPENDENT parameters corresponding to the models.
     *
     */

    ParameterList getSubstitutionModelParameters(bool independent) const;
    
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
    bool checkOrphanNodes(bool throwEx) const;
    
    bool checkUnknownNodes(bool throwEx) const;

  public:

    /*
     * Inheriting from SubstitutionProcess
     */
  
    bool isCompatibleWith(const SiteContainer& data) const;

    bool hasDerivableParameter(const std::string& name) const;

    /**
     * @brief Get the number of states associated to this model set.
     *
     * @return The number of states, or 0 if no model is associated to
     * the set.
     */
    
    size_t getNumberOfStates() const
    {
      if (modelSet_.size()==0)
        return 0;
      else
        return modelSet_[0]->getNumberOfStates();
    }

    /**
     * @return The values of the root frequencies.
     */
  
    const std::vector<double>& getRootFrequencies() const
    {
      if (stationarity_)
        return modelSet_[0]->getFrequencies();
      else
        return rootFrequencies_->getFrequencies();
    }

    /**
     * @brief Get the substitution model corresponding to a certain branch, site pattern, and model class.
     *
     * @param nodeId The id of the node.
     * @param classIndex The model class index.
     */

    const SubstitutionModel& getSubstitutionModel(int nodeId, size_t classIndex) const
    {
      return *modelSet_[nodeToModel_[nodeId]];
    }
    
    const Matrix<double>& getGenerator(int nodeId, size_t classIndex) const
    {
      return getSubstitutionModel(nodeId, classIndex).getGenerator();
    }

    /**
     * This method is used to initialize likelihoods in reccursions.
     * It typically sends 1 if i = state, 0 otherwise, where
     * i is one of the possible states of the alphabet allowed in the model
     * and state is the observed state in the considered sequence/site.
     *
     * The model used is the first one in the list of the models. 
     *
     * @param i the index of the state in the model.
     * @param state An observed state in the sequence/site.
     * @return 1 or 0 depending if the two states are compatible.
     * @throw BadIntException if states are not allowed in the associated alphabet.
     * @see getStates();
     * @see SubstitutionModel
     */

    double getInitValue(size_t i, int state) const throw (BadIntException)
    {
      if (modelSet_.size()==0)
        throw Exception("NonHomogeneousSubstitutionProcess::getInitValue : no model associated");
      else
        return modelSet_[0]->getInitValue(i,state);
    }
    
    double getProbabilityForModel(size_t classIndex) const
    {
      if (classIndex >= rDist_->getNumberOfCategories())
        throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::getProbabilityForModel.", classIndex, 0, rDist_->getNumberOfCategories());
      return rDist_->getProbability(classIndex);
    }

    Vdouble getClassProbabilities() const
    {
      Vdouble vProb;

      for (size_t i=0;i<rDist_->getNumberOfCategories(); i++)
        vProb.push_back(rDist_->getProbability(i));

      return vProb;
    }

    double getRateForModel(size_t classIndex) const {
      if (classIndex >= rDist_->getNumberOfCategories())
        throw IndexOutOfBoundsException("NonHomogeneousSubstitutionProcess::getRateForModel.", classIndex, 0, rDist_->getNumberOfCategories());
      return rDist_->getCategory(classIndex);
    }

    const ComputingTree& getComputingTree() const
    {
      return *computingTree_.get();
    }
  
    ComputingTree& getComputingTree()
    {
      return *computingTree_.get();
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
     * @param rootFreqs A FrequenciesSet object to parametrize root frequencies
     *        (0 if stationary).
     * 
     * @param tree      The tree to use for the construction of the set.
     */
    
    static NonHomogeneousSubstitutionProcess* createHomogeneousSubstitutionProcess(
      SubstitutionModel* model,
      DiscreteDistribution* rdist,
      FrequenciesSet* rootFreqs,
      ParametrizableTree* tree
      );

    /**
     * @brief Create a NonHomogeneousSubstitutionProcess object, with one model per branch.
     *
     * All branches share the same type of model, but allow one set of parameters per branch.
     * This is also possible to specify some parameters to be common to all branches.
     *
     * @param model                The model to use.
     * @param rdist                The rate distribution
     * @param rootFreqs            A FrequenciesSet object to parametrize root frequencies.
     * @param tree                 The tree to use for the construction of the set.
     * @param globalParameterNames Common parameters for all branches.
     * All other parameters will be considered distinct for all branches.
     */
    static NonHomogeneousSubstitutionProcess* createNonHomogeneousSubstitutionProcess(
      SubstitutionModel* model,
      DiscreteDistribution* rdist,
      FrequenciesSet* rootFreqs,
      ParametrizableTree* tree,
      const std::vector<std::string>& globalParameterNames
      );
    


    /** @} */
  };
} // end of namespace bpp.

#endif // _NONHOMOGENEOUSSUBSTITUTIONPROCESS_H_

 
