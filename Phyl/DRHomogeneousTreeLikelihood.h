//
// File: DRHomogeneousTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Oct 17 18:14:51 2003
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

#ifndef _DRHOMOGENEOUSTREELIKELIHOOD_H_
#define _DRHOMOGENEOUSTREELIKELIHOOD_H_

#include "AbstractHomogeneousTreeLikelihood.h"
#include "DRASDRTreeLikelihoodData.h"
#include "NNISearchable.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/BrentOneDimension.h>

/**
 * @brief Compute likelihood for a 4-tree.
 *
 * This class is used internally by DRHomogeneousTreeLikelihood to test NNI movements.
 * This function needs:
 * - two likelihood arrays corresponding to the conditional likelihoods at top and bottom nodes,
 * - a substitution model and a rate distribution, whose parameters will not be estimated but taken "as is",
 * It takes only one parameter, the branch length.
 */
class BranchLikelihood :
  public virtual Function,
  public virtual AbstractParametrizable
{
  protected:
    const VVVdouble *_array1, *_array2;
    VVVdouble _arrayTmp;
    const SubstitutionModel * _model;
    const DiscreteDistribution * _rDist;
    unsigned int _nbStates, _nbClasses;
    VVVdouble _pxy;
    double _lnL;
    vector<unsigned int> _weights;

  public:
    BranchLikelihood(const vector<unsigned int> & weights):
      _array1(NULL), _array2(NULL), _arrayTmp(),
      _model(NULL), _rDist(NULL),
      _nbStates(0), _nbClasses(0),
      _pxy(), _lnL(log(0.)), _weights(weights)
    {
      _parameters.addParameter(Parameter("BrLen", 1, NULL));
    }
    virtual ~BranchLikelihood() {}

#if defined(VIRTUAL_COV)
    BranchLikelihood * clone() const { return new BranchLikelihood(*this); }
#else
    Clonable * clone() const { return new BranchLikelihood(*this); }
#endif

  public:
    void initModel(const SubstitutionModel *model, const DiscreteDistribution *rDist)
    {
      _model = model;
      _rDist = rDist;
      _nbStates = model->getNumberOfStates();
      _nbClasses  = rDist->getNumberOfCategories();
      _pxy.resize(_nbClasses);
      for(unsigned int i = 0; i < _nbClasses; i++)
      {
        _pxy[i].resize(_nbStates);
        for(unsigned int j = 0; j < _nbStates; j++)
          _pxy[i][j].resize(_nbStates);
      }
    }
    /**
     * @warning No checking on alphabet size or number of rate classes is performed,
     * use with care!
     */
    void initLikelihoods(const VVVdouble *array1, const VVVdouble *array2)
    {
      _array1 = array1;
      _array2 = array2;
    }

    void resetLikelihoods()
    {
      _array1 = NULL;
      _array2 = NULL;
    }

    void setParameters(const ParameterList &parameters)
      throw (ParameterNotFoundException, ConstraintException)
    {
      _parameters.setParametersValues(parameters);
      fireParameterChanged(parameters);
    }

    double getValue() const throw (Exception) { return _lnL; }

    void fireParameterChanged(const ParameterList & parameters)
    {
      computeAllTransitionProbabilities();
      computeLogLikelihood();
    }

  protected:
    void computeAllTransitionProbabilities();
    void computeLogLikelihood();
};



/**
 * @brief This class implements the likelihood computation for a tree using the double-recursive
 * algorithm.
 *
 * The substitution model is the same over the tree (homogeneous model).
 * A non-uniform distribution of rates among the sites is allowed (ASRV models).</p>
 *
 * The Felsenstein recursive algorithm is used for conputation.
 * As in the HomogeneousTreeLikelihood class, a likelihood tensor is defined:
 * 
 * -Site
 * -Rate class
 * -Ancestral state
 * 
 * However, in this class, a node will be attached a set of tensor instead of one single tensor,
 * one tensor for each subtree it defines.
 * These tensors are stored into a map of map, with each node as a primary key and each neighbor
 * of the node (sons + father) a a secondary key.
 *
 * All nodes share the same site patterns.
 */
class DRHomogeneousTreeLikelihood :
  public AbstractHomogeneousTreeLikelihood,
  public NNISearchable
{
  protected:
    mutable DRASDRTreeLikelihoodData *_likelihoodData;
    BranchLikelihood * _brLikFunction;
    /**
     * @brief Optimizer used for testing NNI.
     */
    BrentOneDimension * _brentOptimizer;

    /**
     * @brief Hash used for backing up branch lengths when testing NNIs.
     */
    mutable map<int, double> _brLenNNIValues;

    ParameterList _brLenNNIParams;
    
  public:
    /**
     * @brief Build a new DRHomogeneousTreeLikelihood object.
     *
     * @param tree The tree to use.
     * @param model The substitution model to use.
     * @param rDist The rate across sites distribution to use.
     * @param checkRooted Tell if we have to check for the tree to be unrooted.
     * If true, any rooted tree will be unrooted before likelihood computation.
     * @param verbose Should I display some info?
     * @throw Exception in an error occured.
     */
    DRHomogeneousTreeLikelihood(
      const Tree & tree,
      SubstitutionModel * model,
      DiscreteDistribution * rDist,
      bool checkRooted = true,
      bool verbose = true)
      throw (Exception);
  
    /**
     * @brief Build a new DRHomogeneousTreeLikelihood object.
     *
     * @param tree The tree to use.
     * @param data Sequences to use.
     * @param model The substitution model to use.
     * @param rDist The rate across sites distribution to use.
     * @param checkRooted Tell if we have to check for the tree to be unrooted.
     * If true, any rooted tree will be unrooted before likelihood computation.
     * @param verbose Should I display some info?
     * @throw Exception in an error occured.
     */
    DRHomogeneousTreeLikelihood(
      const Tree & tree,
      const SiteContainer & data,
      SubstitutionModel * model,
      DiscreteDistribution * rDist,
      bool checkRooted = true,
      bool verbose = true)
      throw (Exception);

    /**
     * @brief Copy constructor.
     */ 
    DRHomogeneousTreeLikelihood(const DRHomogeneousTreeLikelihood & lik);
    
    DRHomogeneousTreeLikelihood & operator=(const DRHomogeneousTreeLikelihood & lik);

    virtual ~DRHomogeneousTreeLikelihood();

#if defined(VIRTUAL_COV)
    DRHomogeneousTreeLikelihood * clone() const { return new DRHomogeneousTreeLikelihood(*this); }
#else
    Clonable * clone() const { return new DRHomogeneousTreeLikelihood(*this); }
#endif

  public:

    /**
     * @name The TreeLikelihood interface.
     *
     * Other methods are implemented in the AbstractTreeLikelihood class.
     *
     * @{
     */
    void setData(const SiteContainer & sites) throw (Exception);
    double getLikelihood () const;
    double getLogLikelihood() const;
    double getLikelihoodForASite (unsigned int site) const;
    double getLogLikelihoodForASite(unsigned int site) const;
    /** @} */


    void computeTreeLikelihood();

    
    /**
     * @name The DiscreteRatesAcrossSites interface implementation:
     *
     * @{
     */
    double getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
    double getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
    double getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const;
    double getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const;
    /** @} */
  
    /**
     * @brief Implements the Function interface.
     *
     * Update the parameter list and call the fireParameterChanged() method.
     *
     * If a subset of the whole parameter list is passed to the function,
     * only these parameters are updated and the other remain constant (i.e.
     * equal to their last value).
     *
     * @param parameters The parameter list to pass to the function.
     */
    void setParameters(const ParameterList & parameters) throw (ParameterNotFoundException, ConstraintException);
    
    /**
     * @brief Function and NNISearchable interface.
     */
    double getValue() const throw (Exception);
    
    /**
     * @name DerivableFirstOrder interface.
     *
     * @{
     */
    double getFirstOrderDerivative(const string & variable) const throw (Exception);
    /** @{ */

    /**
     * @name DerivableSecondOrder interface.
     *
     * @{
     */
    double getSecondOrderDerivative(const string & variable) const throw (Exception);
    double getSecondOrderDerivative(const string & variable1, const string & variable2) const throw (Exception) { return 0; } // Not implemented for now.
    /** @} */

    /**
     * @name The NNISearchable interface.
     *
     * Current implementation:
     * When testing a particular NNI, only the branch length of the parent node is optimized (and roughly).
     * All other parameters (substitution model, rate distribution and other branch length are kept at there current value.
     * When performing a NNI, only the topology change is performed.
     * This is up to the user to re-initialize the underlying likelihood data to match the new topology.
     * Usually, this is achieved by calling the topologyChangePerformed() method, which call the reInit() method of the LikelihoodData object.
     * @{
     */
		Tree * getTopology() { return AbstractHomogeneousTreeLikelihood::getTree(); }
		
		const Tree * getTopology() const { return getTree(); }

    double getTopologyValue() const throw (Exception) { return getValue(); }

    double testNNI(int nodeId) const throw (NodeException);
    
    void doNNI(int nodeId) throw (NodeException);

    void topologyChangeTested(const TopologyChangeEvent & event)
    {
      _likelihoodData->reInit();
      //if(_brLenNNIParams.size() > 0)
        fireParameterChanged(_brLenNNIParams);
      _brLenNNIParams.reset();
    }
    void topologyChangeSuccessful(const TopologyChangeEvent & event)
    {
      _brLenNNIValues.clear();
    }
    /** @} */
    
  public:  // Specific methods:

    DRASDRTreeLikelihoodData * getLikelihoodData() { return _likelihoodData; }
    const DRASDRTreeLikelihoodData * getLikelihoodData() const { return _likelihoodData; }
  
    virtual VVVdouble computeLikelihoodAtNode(const Node * node) const;

    /**
     * @brief Retrieves all Pij(t) for a particular node.
     *
     * These intermediate results may be used by other methods.
     */
    virtual const VVVdouble & getTransitionProbabilitiesForNode(const Node * node) const { return _pxy[node->getId()]; }
       
  protected:
  
    /**
     * Initialize the arrays corresponding to each son node for the node passed as argument.
     * The method is called for each son node and the result stored in the corresponding array.
     */
    virtual void computeSubtreeLikelihoodPostfix(const Node * node); //Recursive method.
    /**
     * This method initilize the remaining likelihood arrays, corresponding to father nodes.
     * It must be called after the postfix method because it requires that the arrays for
     * son nodes to be be computed.
     */
    virtual void computeSubtreeLikelihoodPrefix(const Node * node); //Recursive method.

    virtual void computeRootLikelihood();

    virtual void computeTreeDLikelihoodAtNode(const Node * node);
    virtual void computeTreeDLikelihoods();
    
    virtual void computeTreeD2LikelihoodAtNode(const Node * node);
    virtual void computeTreeD2Likelihoods();

    void fireParameterChanged(const ParameterList & params);

    void resetLikelihoodArrays(const Node * node);
  
    /**
     * @brief This method is mainly for debugging purpose.
     *
     * @param node The node at which likelihood values must be displayed.
     */
    virtual void displayLikelihood(const Node * node);

    /**
     * @brief Compute conditional likelihoods.
     *
     * This method is the "core" likelihood computation function, performing all the product uppon all nodes, the summation for each ancestral state and each rate class.
     * It is designed for inner usage, and a maximum efficiency, so no checking is performed on the input parameters.
     * Use with care!
     * 
     * @param iLik A vector of likelihood arrays, one for each conditional node.
     * @param tProb A vector of transition probabilities, one for each node.
     * @param oLik The likelihood array to store the computed likelihoods.
     * @param nbNodes The number of nodes = the size of the input vectors.
     * @param nbDistinctSites The number of distinct sites (the first dimension of the likelihood array).
     * @param nbClasses The number of rate classes (the second dimension of the likelihood array).
     * @param nbStates The number of states (the third dimension of the likelihood array).
     * @param reset Tell if the output likelihood array must be initalized prior to computation.
     * If true, the resetLikelihoodArray method will be called.
     */
    static void computeLikelihoodFromArrays(const vector<const VVVdouble *> & iLik, const vector<const VVVdouble *> & tProb, VVVdouble & oLik, unsigned int nbNodes, unsigned int nbDistinctSites, unsigned int nbClasses, unsigned int nbStates, bool reset = true);

};

#endif  //_DRHOMOGENEOUSTREELIKELIHOOD_H_

