//
// File: SubstitutionProcessCollectionMember.h
// Created by: Laurent Guéguen
// Created on: mercredi 13 mai 2015, à 22h 32
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

#ifndef _SUBSTITUTIONPROCESSCOLLECTIONMEMBER_H_
#define _SUBSTITUTIONPROCESSCOLLECTIONMEMBER_H_


#include "SubstitutionProcess.h"

#include "ComputingTree.h"

namespace bpp
{
  /**
   * @brief A substitution process which objects belong to a SubstitutionProcessCollection.
   *
   * This object is a link between SubstitutionProcessCollection and
   * ComputingTree, and does not have any specific Parameter.
   *
   * The parameters are the INDEPENDENT parameters of the objects of
   * the Collection.
   */

  class SubstitutionProcessCollection;
  
  class SubstitutionProcessCollectionMember :
    public virtual SubstitutionProcess,
    public virtual AbstractParameterAliasable
  {
  private:

    /**
     * @brief A pointer towards the collection the SubstitutionProcessCollectionMember
     * belongs to.
     *
     */
  
    SubstitutionProcessCollection* pSubProColl_;

    /**
     * @brief The number of the process in the collection
     *
     **/

    size_t nProc_;

  private:

    /**
     * @brief sets the parameters as the independent parameters on the
     * objects
     *
     **/
    
    void updateParameters();

  private:

    /**
     * @brief Contains for each node in a tree the index of the corresponding model in modelSet_
     */

    std::map<int, size_t> nodeToModel_;
    std::map<size_t, std::vector<int> > modelToNodes_;

    /**
     *@brief The number of the tree
     *
     */

    size_t nTree_;

    /**
     *@brief The number of the rate distribution
     *
     */

    size_t nDist_;

    /**
     * @brief A boolean if the model is stationary, and the number of
     *  the root frequencies.
     *
     */

    bool stationarity_;

    size_t nRoot_;

    /**
     * @brief The related Computing Tree
     *
     */

    mutable std::auto_ptr<ComputingTree> computingTree_;
    
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
  
    SubstitutionProcessCollectionMember( SubstitutionProcessCollection* pSubProColl, size_t nProc, size_t nTree, size_t nDist);

    /**
     * @brief Resets all the information contained in this object.
     *
     */
   
    void clear();
  
    SubstitutionProcessCollectionMember(const SubstitutionProcessCollectionMember& set);

    SubstitutionProcessCollectionMember& operator=(const SubstitutionProcessCollectionMember& set);

    virtual ~SubstitutionProcessCollectionMember() {}

    SubstitutionProcessCollectionMember* clone() const { return new SubstitutionProcessCollectionMember(*this); }

  private:

    /**
     * @brief Method to inform changes in models.
     *
     */

    void changedModel(size_t m)
    {
      if (modelToNodes_.size()==1)
        computingTree_->updateAll(); // faster
      else
        computingTree_->update(modelToNodes_[m]);
    }

    /**
     * @brief Method to inform changes in root
     *
     * If flag = true (default), node has to be updated (false otherwise).
     */
    
    void changedRoot(bool flag = true)
    {
      computingTree_->update((*computingTree_)[0]->getRootId(), flag);
    }
    
  public:

    const Alphabet* getAlphabet() const;

    /**
     * @return the number of the process in the collection
     */

    size_t getNProcess() const
    {
      return nProc_;
    }
    
    /**
     * @return The current number of distinct substitution models in this set.
     */
    size_t getNumberOfModels() const { return modelToNodes_.size(); }

    /**
     * @return True iff there is a MixedSubstitutionModel in the SubstitutionProcessCollectionMember
     **/

    bool hasMixedSubstitutionModel() const;

    /**
     * @return True iff is stationary.
     **/

    bool isStationary() const
    {
      return stationarity_;
    }

    /**
     * @brief Get one model from the set knowing its index.
     *
     * @param i Index of the model in the set.
     * @return A pointer toward the corresponding model.
     */
  
    const SubstitutionModel* getModel(size_t i) const;

    std::vector<size_t> getModelNumbers() const;

    /**
     * @brief Get the index in the set of the model associated to a particular node id.
     *
     * @param nodeId The id of the query node.
     * @return The index of the model associated to the given node.
     * @throw Exception If no model is found for this node.
     */

    size_t getModelIndexForNode(int nodeId) const throw (Exception)
    {
      std::map<int, size_t>::const_iterator i = nodeToModel_.find(nodeId);
      if (i == nodeToModel_.end())
        throw Exception("SubstitutionProcessCollectionMember::getModelIndexForNode(). No model associated to node with id " + TextTools::toString(nodeId));
      return i->second;
    }

    /**
     * @brief Get the model associated to a particular node id.
     *
     * @param nodeId The id of the query node.
     * @return A pointer toward the corresponding model.
     * @throw Exception If no model is found for this node.
     */

    const SubstitutionModel* getModelForNode(int nodeId) const throw (Exception);
  
    /**
     * @brief Get a list of nodes id for which the given model is associated.
     *
     * @param i The index of the model in the set.
     * @return A vector with the ids of the node associated to this model.
     * @throw IndexOutOfBoundsException If the index is not valid.
     */
  
    const std::vector<int>& getNodesWithModel(size_t i) const 
    {
      std::map<size_t, std::vector<int> >::const_iterator it = modelToNodes_.find(i);
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

    void addModel(size_t numModel, const std::vector<int>& nodesId);

    /**
     * @brief Get the rate distribution
     *
     **/

    const DiscreteDistribution* getRateDistribution() const;

    const size_t getRateDistributionNumber() const { return nDist_;}
    
    /*
     * @brief Set the root Frequencies Set
     * @param freqIndex the index of the frequencies in the collection.
     *
     */

    void setRootFrequencies(size_t numFreq);

    const size_t getRootFrequenciesNumber() const { return nRoot_;}

    /**
     * @return The set of root frequencies.
     *
     */
  
    const FrequenciesSet* getRootFrequenciesSet() const;

    /**
     * @brief AbsractParametrizable interface
     *
     **/

    bool matchParametersValues(const ParameterList& parameters) throw (bpp::ConstraintException);
    
    void fireParameterChanged(const ParameterList& parameters)
    {
      computingTree_->matchParametersValues(parameters);
    }
    
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
    bool checkOrphanNodes(bool throwEx) const throw (Exception);

    bool checkUnknownNodes(bool throwEx) const throw (Exception);
    /** @} */

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
    
    size_t getNumberOfStates() const;

    /**
     * @return The values of the root frequencies.
     */
  
    const std::vector<double>& getRootFrequencies() const;

    /**
     * @return the Tree
     */
    
    const TreeTemplate<Node>& getTree() const;
    
    const ParametrizableTree& getParametrizableTree() const;

    size_t getTreeNumber() const { return nTree_;}
    
    size_t getNumberOfClasses() const;
    
    /**
     * @brief Get the substitution model corresponding to a certain branch, site pattern, and model class.
     *
     * @param nodeId The id of the node.
     * @param classIndex The model class index.
     */

    const SubstitutionModel& getSubstitutionModel(int nodeId, size_t classIndex) const;

    /**
     * @brief Get the parameters of the substitution models.
     *
     **/
     
    ParameterList getSubstitutionModelParameters(bool independent) const;
    
    /**
     * @brief Get the parameters of the rate distribution.
     *
     **/

    ParameterList getRateDistributionParameters(bool independent) const;

    /**
     * @brief Get the parameters of the tree.
     *
     **/

    ParameterList getBranchLengthParameters(bool independent) const;

    bool hasBranchLengthParameter(const std::string& name) const;

    /**
     * @brief Get the parameters of the root frequencies set.
     *
     **/

    ParameterList getRootFrequenciesParameters(bool independent) const;

    /**
     * @brief get (Non)Derivable INDEPENDENT parameters
     *
     **/
    
    ParameterList getDerivableParameters() const;

    ParameterList getNonDerivableParameters() const;


    /**
     * @brief Get the transition probabilities corresponding to a
     * certain branch, and model class.
     *
     * @param nodeId The id of the node.
     * @param classIndex The model class index.
     */
    const Matrix<double>& getTransitionProbabilities(int nodeId, size_t classIndex) const
    {
      return (*computingTree_)[classIndex]->getNode(nodeId)->getTransitionProbabilities();
    }
 
    /**
     * @brief Get the first order derivatives of the transition
     * probabilities according to time, corresponding to a certain
     * branch, and model class.
     *
     * @param nodeId The id of the node.
     * @param classIndex The model class index.
     */
    const Matrix<double>& getTransitionProbabilitiesD1(int nodeId, size_t classIndex) const
    {
      return (*computingTree_)[classIndex]->getNode(nodeId)->getTransitionProbabilitiesD1();
    }
 
    /**
     * @brief Get the second order derivatives of the transition
     * probabilities according to time, corresponding to a certain
     * branch, and model class.
     *
     * @param nodeId The id of the node.
     * @param classIndex The model class index.
     */
    const Matrix<double>& getTransitionProbabilitiesD2(int nodeId, size_t classIndex) const
    {
      return (*computingTree_)[classIndex]->getNode(nodeId)->getTransitionProbabilitiesD2();
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

    double getInitValue(size_t i, int state) const throw (BadIntException);
    
    double getProbabilityForModel(size_t classIndex) const;

    Vdouble getClassProbabilities() const;

    double getRateForModel(size_t classIndex) const;

    /**
     * @brief Methods for computing partial likelihoods. See
     * class ComputingTree for details.
     *
     **/

    /**
     * @brief A virtual method to retrieve the ComputingTree defined in
     * inheriting classes.
     *
     */
  
    const ComputingTree& getComputingTree() const
    {
      return *computingTree_.get();
    }
    
    ComputingTree& getComputingTree()
    {
      return *computingTree_.get();
    }

    friend class SubstitutionProcessCollection;
    
  };
} // end of namespace bpp.

#endif // _SUBSTITUTIONPROCESSCOLLECTIONMEMBER_H_

