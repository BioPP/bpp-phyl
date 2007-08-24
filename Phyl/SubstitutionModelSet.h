//
// File: SubstitutionModelSet.h
// Created by: Bastien Bousseau
//             Julien Duteil
// Created on: Tue Aug 21 2007
//

/*
Copyright or <A9> or Copr. CNRS, (November 16, 2004)

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


#include "Phyl/SubstitutionModel.h"
#include "Phyl/AbstractSubstitutionModel.h"

// From NumCalc:
#include <NumCalc/RandomTools.h>

// From utils:
#include <Utils/Exceptions.h>

//From Seqlib:
#include <Seq/Alphabet.h>

// From the STL:
#include <vector>
#include <map>
#include <algorithm>
using namespace std;

/**
 * @brief Substitution models manager for heterogenous models of evolution.
 *
 * This class contains a set of substitution models, and their assigment towoard the branches of a phylogenetic tree.
 * Each branch in the tree corresponds to a model in the set, but a susbstitution model may correspond to several branches.
 * The particular case where all branche points toward a unique model is the homogeneous case.
 *
 * This class also deals with the parameters associated to the models.
 * In the homogeneous case, the parameter list is the same as the list in susbstitution model.
 * When two models at least are specified, these model may have their own parameters or share some of them.
 * To deal with this issue, the SubstitutionModelSet class contains its own parameter list and an index wich tells to which
 * models this parameter applies to.
 * Since parameters in a list must have unique names, duplicated names are numbered according to the order in the list.
 * To track the relation between names in the list and names in each model, the parameter list is duplicated in _modelParameters.
 * The user only act on _parameters, the fireParameterChanged function, automatically called, will update the _modelParameters field.
 */
class SubstitutionModelSet: 
  public AbstractParametrizable
{
  protected:
    /**
     * @brief contains all models used in this tree.
     */
    vector<SubstitutionModel *> _modelSet;

    /**
     * @brief contains for each node in a tree the index of the corresponding model in _modelSet
     */
    mutable map<int, unsigned int> _nodeToModel;

    /**
     * @brief contains for each parameter in the list the indexes of the corresponding models in _modelSet that share this parameter.
     */
    vector< vector<unsigned int> > _paramToModels;

    /**
     * @brief contains for each parameter in the list the corresponding name in substitution models.
     *
     */
    vector<string> _modelParameterNames;

    /**
     * @brief Parameters for each model in the set.
     *
     * The _parameters field, inherited from AbstractSubstitutionModel contains all parameters, with unique names.
     * To make the correspondance with parameters for each model in the set, we duplicate them in this array.
     * In most cases, this is something like 'theta0 <=> theta', 'theta2 <=> theta', etc.
     */
    vector<ParameterList> _modelParameters;

    /**
     * @brief A pointer toward the comon alphabet to all models in the set.
     */
    const Alphabet * _alphabet;

    IncludingInterval _freqConstraint;

    /**
     * @brief Root frequencies.
     */
    vector <double> _rootFrequencies;

  public:
  
    SubstitutionModelSet(const Alphabet *alpha): _alphabet(alpha), _freqConstraint(0, 1)
    {
      double temp = 1. / (double)(_alphabet->getSize());
      for(unsigned int i = 0; i < _alphabet->getSize(); i++)
      {
        _parameters.addParameter(Parameter("RootFreq" + TextTools::toString(i), temp, & _freqConstraint));
        _modelParameterNames.push_back("RootFreq" + TextTools::toString(i)); //Just to maintain correspondance between indices. These names are actually never used.
      }
      _rootFrequencies.resize(_alphabet->getSize());
      updateRootFrequencies();
    }

    SubstitutionModelSet(const Alphabet *alpha, const vector<double> & rootFreqs): _alphabet(alpha), _freqConstraint(0, 1), _rootFrequencies(rootFreqs) {}

    SubstitutionModelSet(): _alphabet(NULL), _freqConstraint(0, 1) {}


//Plutot en classe utilitaire
/*  SubstitutionModelSet(map<string, string> & params)
  {
    _alphabet = SequenceApplicationTools::getAlphabet(params, "", false);
    if (_alphabet->getSize()==4){
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancA", params, 0.25));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancC", params, 0.25));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancG", params, 0.25));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancT", params, 0.25));
    }
    else if (_alphabet->getSize()==20) {
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancA", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancC", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancD", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancE", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancF", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancG", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancH", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancI", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancK", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancL", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancM", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancN", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancP", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancQ", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancR", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancS", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancT", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancV", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancW", params, 0.05));
      _initialFrequencies.push_back(ApplicationTools::getDoubleParameter("ancY", params, 0.05));
    }
    else {
      for (int i=0;i<_alphabet->getSize();i++) {
	double temp =1.0/(double)(_alphabet->getSize());
	_initialFrequencies.push_back(temp);
      }
    }

    double sum=0.0;
    for (int i=0;i<_alphabet->getSize();i++) {
      sum += _initialFrequencies[i];
    }
    if( fabs(1-(sum)) > 0.00000000000001 )
      {
	throw Exception("Equilibrium base frequencies must equal 1. Aborting...");
      }
  }*/

    SubstitutionModelSet(const SubstitutionModelSet & set);
    
    SubstitutionModelSet & operator=(const SubstitutionModelSet & set);
    
    virtual ~SubstitutionModelSet()
    {
      for(unsigned int i = 0; i < _modelSet.size(); i++) delete _modelSet[i];
    }

    SubstitutionModelSet * clone() const { return new SubstitutionModelSet(*this); }

  public:

    /**
     * To be called when a parameter has changed.
     * Depending on parameters, this will actualize the _initialFrequencies vector or the corresponding models in the set.
     * @param parameters The modified parameters.
     */
    void fireParameterChanged(const ParameterList & parameters);

    //Dans classe static!
    //void makeRandomNodeToModelLinks();
  
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
     * @throw Exception in case of error (...).
     */
    void addModel(SubstitutionModel * model, const vector<int> & nodesId, const vector<string> & newParams) throw (Exception);
  
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
    void setModel(SubstitutionModel * model, unsigned int modelIndex) throw (Exception, IndexOutOfBoundsException);
 
    /**
     * @brief Associate an existing model with a given node.
     *
     * If the node was was previously associated to a model, the old association is deleted.
     * If other nodes are associated to this model, the association is conserved.
     *
     * @param modelIndex The position of the model in the set.
     * @param nodeNumber The id of the corresponding node.
     */
    void setModelToNode(unsigned int modelIndex, int nodeNumber) throw (IndexOutOfBoundsException)
    {
      if(modelIndex >= _nodeToModel.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::setModelToNode.", modelIndex, 0, _nodeToModel.size() - 1);
      _nodeToModel[nodeNumber] = modelIndex;
    }
    
    void setParameterToModel(unsigned int parameterIndex, unsigned int modelIndex) throw (IndexOutOfBoundsException)
    {
      if(parameterIndex >= _paramToModels.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::setParameterToModel.", parameterIndex, 0, _paramToModels.size() - 1);
      _paramToModels[parameterIndex].push_back(modelIndex);
    }

    void unsetParameterToModel(unsigned int parameterIndex, unsigned int modelIndex) throw (IndexOutOfBoundsException)
    {
      if(parameterIndex >= _paramToModels.size()) throw IndexOutOfBoundsException("SubstitutionModelSet::unsetParameterToModel.", parameterIndex, 0, _paramToModels.size() - 1);
      remove(_paramToModels[parameterIndex].begin(), _paramToModels[parameterIndex].end(), modelIndex);
      if(!checkOrphanModels())     throw Exception("DEBUG: SubstitutionModelSet::unsetParameterToModel. Orphan model!");
      if(!checkOrphanParameters()) throw Exception("DEBUG: SubstitutionModelSet::unsetParameterToModel. Orphan parameterl!");
    }

    /**
     * @brief Add a parameter to the list, and link it to specified existing nodes.
     *
     * @param parameter The parameter to add. Its name must match model parameter names.
     * @param nodesId The list of ids of the nodes to link with these parameters.
     * Nodes must have a corresponding model in the set.
     * @throw Exception If one of the above requirement is not true.
     */
    void addParameter(const Parameter & parameter, const vector<int> & nodesId) throw (Exception);
 
    /**
     * @brief Add several parameters to the list, and link them to specified existing nodes.
     *
     * @param parameters The list of parameters to add. Its name must match model parameter names.
     * @param nodesId The list of ids of the nodes to link with these parameters.
     * Nodes must have a corresponding model in the set.
     * @throw Exception If one of the above requirement is not true.
     */
    void addParameters(const ParameterList & parameters, const vector<int> & nodesId) throw (Exception);
 
     void removeModel(unsigned int modelIndex)
    {
      _modelSet.erase(_modelSet.begin() + modelIndex);
      if(!checkOrphanParameters()) throw Exception("DEBUG: SubstitutionModelSet::removeModel. Orphan parameter!");
      if(!checkOrphanNodes())      throw Exception("DEBUG: SubstitutionModelSet::removeModel. Orphan node!");
    }

    /**
     * @brief Get a reference toward the model associated to node with the specified id, if any.
     *
     * @param nodeId The id of the node to search.
     * @return A pointer toward the associated model, NULL if no model is associated with this node.
     */
    const SubstitutionModel * getModelForNode(int nodeId) const;
 
    /**
     * @brief Get a reference toward the model associated to node with the specified id, if any.
     *
     * @param nodeId The id of the node to search.
     * @return A pointer toward the associated model, NULL if no model is associated with this node.
     */
    SubstitutionModel * getModelForNode(int nodeId);

    void listModelNames(ostream & out = cout) const;

    //Check sum to 1?
    void setInitialFrequencies(const vector<double> & initFreqs) { _rootFrequencies = initFreqs; }

    vector<double> getRootFrequencies() const { return _rootFrequencies; }

    const Alphabet * getAlphabet() const { return _alphabet; }

  protected:

    /**
     * Set _rootFrequencies from _parameters.
     */
    void updateRootFrequencies();

    //TODO
    bool checkOrphanModels() const { return true; }

    bool checkOrphanParameters() const { return true; }

    bool checkOrphanNodes() const { return true; }

};

#endif // _SUBSTITUTIONMODELSET_H_

