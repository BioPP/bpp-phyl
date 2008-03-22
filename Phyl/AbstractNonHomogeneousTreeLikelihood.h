//
// File: AbstractNonHomogeneousTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Tue Oct 09 16:07 2007
// From file: AbstractHomogeneousTreeLikelihood.h
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _ABSTRACTNONHOMOGENEOUSTREELIKELIHOOD_H_
#define _ABSTRACTNONHOMOGENEOUSTREELIKELIHOOD_H_

#include "NonHomogeneousTreeLikelihood.h"
#include "AbstractDiscreteRatesAcrossSitesTreeLikelihood.h"

namespace bpp
{

/**
 * @brief Partial implementation for non-homogeneous models of the TreeLikelihood interface.
 *
 * This class provides a pointer toward a single substitution model + several utilitary variables.
 */
class AbstractNonHomogeneousTreeLikelihood:
  public virtual NonHomogeneousTreeLikelihood,
  public AbstractDiscreteRatesAcrossSitesTreeLikelihood
{
  protected:
    SubstitutionModelSet * _modelSet;
    ParameterList _brLenParameters;
    
    mutable map<int, VVVdouble> _pxy;

    mutable map<int, VVVdouble> _dpxy;

    mutable map<int, VVVdouble> _d2pxy;
        
    vector<double> _rootFreqs;
        
    /**
     * @brief Pointer toward all nodes in the tree.
     *
     * The position in the array is the number used in the parameter name.
     * This may be different from the node id, unless you used the resetNodeId method on the input tree.
     */
    vector<Node *> _nodes;

    /**
     * @brief An index linking nodes to their id, for faster access than the getNode() method.
     */
    mutable map<int, const Node *> _idToNode;
 
    //some values we'll need:
    unsigned int _nbSites,         //the number of sites in the container
                 _nbDistinctSites, //the number of distinct sites
                 _nbClasses,       //the number of rate classes
                 _nbStates,        //the number of states in the alphabet
                 _nbNodes;         //the number of nodes in the tree

    bool _verbose;

    double _minimumBrLen;
    Constraint * _brLenConstraint;

    int _root1, _root2;

  public:
    AbstractNonHomogeneousTreeLikelihood(
      const Tree & tree,
      SubstitutionModelSet * modelSet,
      DiscreteDistribution * rDist,
      bool verbose = true)
      throw (Exception);

    /**
     * @brief Copy constructor
     *
     * This constructor is to be called by the derived class copy constructor.
     */
    AbstractNonHomogeneousTreeLikelihood(const AbstractNonHomogeneousTreeLikelihood & lik);
    
    /**
     * @brief Assignation operator
     *
     * This operator is to be called by the derived class operator.
     */
    AbstractNonHomogeneousTreeLikelihood & operator=(const AbstractNonHomogeneousTreeLikelihood & lik);
 
    virtual ~AbstractNonHomogeneousTreeLikelihood();
    
  private:

    /**
     * @brief Method called by constructor.
     */
    void _init(const Tree & tree,
      SubstitutionModelSet * modelSet,
      DiscreteDistribution * rDist,
      bool verbose) throw (Exception);

  public:
    
    /**
     * @name The TreeLikelihood interface.
     *
     * Other methods are implemented in the AbstractTreeLikelihood class.
     *
     * @{
     */
    void initialize() throw(Exception);
    
    ParameterList getBranchLengthsParameters() const;
    
    ParameterList getSubstitutionModelParameters() const;
    
    ParameterList getRateDistributionParameters() const
    {
      return AbstractDiscreteRatesAcrossSitesTreeLikelihood::getRateDistributionParameters();
    }

    const SubstitutionModel * getSubstitutionModelForNode(int nodeId) const throw (NodeNotFoundException) 
    {
      return _modelSet->getModelForNode(nodeId);
    }

    SubstitutionModel * getSubstitutionModelForNode(int nodeId) throw (NodeNotFoundException)
    {
      return _modelSet->getModelForNode(nodeId);
    }

    vector<double> getRootFrequencies() const { return _rootFreqs; }
    
    const VVVdouble & getTransitionProbabilitiesForNode(const Node* node) const { return _pxy[node->getId()]; }
    /** @} */

    /**
     * @name The NonHomogeneousTreeLikelihood interface.
     *
     * Other methods are implemented in the AbstractTreeLikelihood class.
     *
     * @{
     */
    const SubstitutionModelSet * getSubstitutionModelSet() const { return _modelSet; }
    
    SubstitutionModelSet * getSubstitutionModelSet() { return _modelSet; }
    
    void setSubstitutionModelSet(SubstitutionModelSet * modelSet) throw (Exception);

    ParameterList getRootFrequenciesParameters() const
    {
      return _modelSet->getRootFrequenciesParameters();
    }
    /** @} */
    
  public: //Specific methods:

    /**
     * @brief This builds the <i>parameters</i> list from all parametrizable objects,
     * <i>i.e.</i> substitution model, rate distribution and tree.
     */
    virtual void initParameters();

    /**
     * @brief All parameters are stored in a parameter list.
     * This function apply these parameters to the substitution model,
     * to the rate distribution and to the branch lengths.
     */
    virtual void applyParameters() throw (Exception);  

    virtual void initBranchLengthsParameters();

    virtual void setMinimumBranchLength(double minimum)
    {
      _minimumBrLen = minimum;
      if(_brLenConstraint != NULL) delete _brLenConstraint;
      _brLenConstraint = new IncludingPositiveReal(_minimumBrLen);
      initBranchLengthsParameters();
    }

    virtual double getMinimumBranchLength() const { return _minimumBrLen; }

  protected:
    /**
     * @brief Fill the _pxy, _dpxy and _d2pxy arrays for all nodes.
     */
    void computeAllTransitionProbabilities();
    /**
     * @brief Fill the _pxy, _dpxy and _d2pxy arrays for one node.
     */
    void computeTransitionProbabilitiesForNode(const Node * node);

};

} //end of namespace bpp.

#endif //_ABSTRACTNONHOMOGENEOUSTREELIKELIHOOD_H_

