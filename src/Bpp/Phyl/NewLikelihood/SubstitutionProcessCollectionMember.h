//
// File: SubstitutionProcessCollectionMember.h
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

#ifndef _SUBSTITUTIONPROCESSCOLLECTIONMEMBER_H_
#define _SUBSTITUTIONPROCESSCOLLECTIONMEMBER_H_


#include "SubstitutionProcessCollection.h"
#include "SubstitutionProcess.h"

#include <Bpp/Exceptions.h>

// From Seqlib:
#include <Bpp/Seq/Alphabet/Alphabet.h>


namespace bpp
{
  /**
   * @brief A substitution process which objects belong to a SubstitutionProcessCollection.
   *
   */
  
  class SubstitutionProcessCollectionMember :
    public SubstitutionProcess
  {
  private:

    /**
     * @brief A pointer towards the collection de SubstitutionProcessCollectionMember
     * belongs to.
     *
     */
  
    const SubstitutionProcessCollection* pSubProColl_;
  

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

    unsigned int nTree_;

    /**
     *@brief The number of the rate distribution
     *
     */

    unsigned int nRate_;

    /**
     * @brief A boolean if the model is stationary, and the number of
     *  the root frequencies.
     *
     */

    bool stationarity_;

    unsigned int nRoot_;

  public:
    /**
     * @brief Create a model set belonging to the specified SubstitutionProcessCollection.
     * Stationarity is assumed.
     *
     * @param pSubProColl the SubstitutionProcessCollection.
     */
  
    SubstitutionProcessCollectionMember(const SubstitutionProcessCollection* pSubProColl) :
    pSubProColl_(pSubProColl),
    nodeToModel_(),
    modelToNodes_(),
    nTree_(0),
    nRate_(0),
    stationarity_(true),
    nRoot_(0)
    {
    }

    /**
     * @brief Resets all the information contained in this object.
     *
     */
   
    void clear();
  
    SubstitutionProcessCollectionMember(const SubstitutionProcessCollectionMember& set);

    SubstitutionProcessCollectionMember& operator=(const SubstitutionProcessCollectionMember& set);

    virtual ~SubstitutionProcessCollectionMember() {}

#ifndef NO_VIRTUAL_COV
    SubstitutionProcessCollectionMember*
#else
    Clonable*
#endif
    clone() const { return new SubstitutionProcessCollectionMember(*this); }

  public:
    /**
     * @return The current number of distinct substitution models in this set.
     */
    size_t getNumberOfClasses() const { return modelToNodes_.size(); }

    /**
     * @return The current number of distinct substitution models in this set.
     */
    size_t getNumberOfStates() const { return modelToNodes_.size(); }
    
    /**
     * @return True iff there is a MixedSubstitutionModel in the SubstitutionProcessCollectionMember
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
      return pSubProColl_->getModel(i);
    }

    SubstitutionModel* getModel(size_t i) 
    {
      return pSubProColl_->getModel(i);
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
    const SubstitutionModel* getModelForNode(int nodeId) const throw (Exception)
    {
      std::map<int, size_t>::const_iterator i = nodeToModel_.find(nodeId);
      if (i == nodeToModel_.end())
        throw Exception("SubstitutionProcessCollectionMember::getModelForNode(). No model associated to node with id " + TextTools::toString(nodeId));
      return getModel(i->second);
    }
  
    SubstitutionModel* getModelForNode(int nodeId) 
    {
      std::map<int, size_t>::iterator i = nodeToModel_.find(nodeId);
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
  
    const std::vector<int>& getNodesWithModel(size_t i) const 
    {
      std::map<size_t, std::vector<int> >::const_iterator it = modelToNodes_.find(i);
      if (it == modelParameters_.end())
        throw Exception("SubstitutionProcessCollectionMember::getNodesWithModel(). No nodes associated with model " + TextTools::toString(i));
    
      return it->second;
    }

    // /**
    //  * @param name The name of the parameter to look for.
    //  * @return The list of nodes with a model containing the specified parameter.
    //  * @throw ParameterNotFoundException If no parameter with the specified name is found.
    //  */

    // std::vector<int> getNodesWithParameter(const std::string& name) const;

    /**
     * @brief Add a new model to the set, and set relationships with nodes.
     *
     * @param numModel The number of a model in the SubstitutionProcessCollection.
     * @param nodesId the set of nodes in the tree that points toward this model.
     * This will override any previous affectation.
     */

    void addModel(unsigned int numModel, const std::vector<int>& nodesId);//, const std::vector<std::string>& newParams) throw (Exception);

    /**
     * @brief Set the tree
     * @param numTree the number of  the tree in the collection.
     *
     **/

    void setTree(unsigned int numTree);

    /**
     * @brief Get the tree
     *
     **/

    Tree* getTree()
    {
      return pSubProColl_->getTree(nTree_);
    }

    const Tree* getTree() const
    {
      return pSubProColl_->getTree(nTree_);
    }

    /**
     * @brief Set the rate distribution
     * @param numRate the number of  the rate in the collection.
     *
     **/

    void setRate(unsigned int numRate);

    DiscreteDistribution* getRate()
    {
      return pSubProColl_->getDistribution(nRate_);
    }

    const DiscreteDistribution* getRate() const
    {
      return pSubProColl_->getDistribution(nRate_);
    }
  

    // /**
    //  * @brief Remove a model from the set, and all corresponding parameters.
    //  *
    //  * @param modelIndex The index of the model in the collection.
    //  * @throw Exception if a parameter becomes orphan because of the removal.
    //  */
    // void removeModel(size_t modelIndex) throw (Exception);

    //  void listModelNames(std::ostream& out = std::cout) const;

    /*
     * @brief Set the root Frequencies Set
     * @param freqIndex the index of the frequencies in the collection.
     *
     */

    void addRootFrequencies(unsigned int numFreq);
  
    /**
     * @return The set of root frequencies.
     *
     */
  
    const FrequenciesSet* getRootFrequenciesSet() const
    {
      if (stationarity_)
        return 0;
      else
        return pSubProColl_->getFrequencies(nRoot_);
    }

    /**
     * @return The values of the root frequencies.
     */
  
    std::vector<double> getRootFrequencies() const
    {
      if (stationarity_)
        return (pSubProColl_->getModel(modelToNodes_.begin()->first))->getFrequencies();
      else
        return (pSubProColl_->getFrequencies(nRoot_))->getFrequencies();
    }

    // /**
    //  * @brief Get the parameters corresponding to the root frequencies.
    //  *
    //  * @return The parameters corresponding to the root frequencies.
    //  */
    // ParameterList getRootFrequenciesParameters() const
    // {
    //   if (stationarity_)
    //     return ParameterList();
    //   else
    //     return rootFrequencies_->getParameters();
    // }

    // /**
    //  * @brief Get the parameters corresponding attached to the nodes of the tree.
    //  *
    //  * That is, all the parameters without the root frequencies.
    //  *
    //  * @return The parameters attached to the nodes of the tree.
    //  */
    // ParameterList getNodeParameters() const
    // {
    //   ParameterList pl;
    //   for (size_t i = stationarity_ ? 0 : rootFrequencies_->getNumberOfParameters();
    //       i < getNumberOfParameters(); i++)
    //   {
    //     pl.addParameter(getParameter_(i));
    //   }
    //   return pl;
    // }

    /**
     * @brief Get the parameters attached to a Model.
     *
     * @param modelIndex the index of the model in the set
     *
     * @return The parameters attached to the model.
     */

    //  ParameterList getModelParameters(size_t modelIndex) const;

    const Alphabet* getAlphabet() const
    {
      return pSubProColl_->getModel(modelToNodes_.begin()->first)->getAlphabet();
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
    bool isFullySetUpFor(const Tree& tree, bool throwEx = true) const throw (Exception)
    {
      return checkOrphanModels(throwEx)
        && checkOrphanNodes(tree, throwEx)
        && checkUnknownNodes(tree, throwEx);
    }

    /**
     *  @brief void function for backward compatibilities. To be removed afterwards.
     *
     */
   
    void matchParametersValues(ParameterList&) {}
  
  protected:
    /**
     //  * Set rootFrequencies_ from parameters.
     //  */
    // void updateRootFrequencies()
    // {
    //   if (!stationarity_)
    //     rootFrequencies_->matchParametersValues(getParameters());
    // }

    /**
     * @name Check function.
     *
     * @{
     */
    bool checkOrphanModels(bool throwEx) const throw (Exception);

    //  bool checkOrphanParameters(bool throwEx) const throw (Exception);

    bool checkOrphanNodes(const Tree& tree, bool throwEx) const throw (Exception);

    bool checkUnknownNodes(const Tree& tree, bool throwEx) const throw (Exception);
    /** @} */
  };
} // end of namespace bpp.

#endif // _SUBSTITUTIONPROCESSCOLLECTIONMEMBER_H_

