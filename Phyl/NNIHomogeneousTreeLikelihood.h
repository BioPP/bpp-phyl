//
// File: NNIHomogeneousTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Apr 06 14:16 2007
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

#ifndef _NNIHOMOGENEOUSTREELIKELIHOOD_H_
#define _NNIHOMOGENEOUSTREELIKELIHOOD_H_

#include "DRHomogeneousTreeLikelihood.h"
#include "NNISearchable.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/BrentOneDimension.h>
#include <NumCalc/Parametrizable.h>

namespace bpp
{

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
  public Function,
  public AbstractParametrizable
{
  protected:
    const VVVdouble *_array1, *_array2;
    VVVdouble _arrayTmp;
    const SubstitutionModel * _model;
    const DiscreteDistribution * _rDist;
    unsigned int nbStates_, nbClasses_;
    VVVdouble pxy_;
    double _lnL;
    vector<unsigned int> _weights;

  public:
    BranchLikelihood(const vector<unsigned int> & weights) :
      AbstractParametrizable(""),
      _array1(NULL), _array2(NULL), _arrayTmp(),
      _model(NULL), _rDist(NULL),
      nbStates_(0), nbClasses_(0),
      pxy_(), _lnL(log(0.)), _weights(weights)
    {
      Parameter p("BrLen", 1, NULL);
      addParameter_(p);
    }
    virtual ~BranchLikelihood() {}

#ifndef NO_VIRTUAL_COV
    BranchLikelihood*
#else
    Clonable*
#endif
    clone() const { return new BranchLikelihood(*this); }

  public:
    void initModel(const SubstitutionModel *model, const DiscreteDistribution *rDist);

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
      //parameters_.setParametersValues(parameters);
      //fireParameterChanged(parameters);
      setParametersValues(parameters);
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
 * @brief This class adds support for NNI topology estimation to the DRHomogeneousTreeLikelihood class.
 */
class NNIHomogeneousTreeLikelihood:
  public DRHomogeneousTreeLikelihood,
  public NNISearchable
{
  protected:
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
     * @brief Build a new NNIHomogeneousTreeLikelihood object.
     *
     * @param tree The tree to use.
     * @param model The substitution model to use.
     * @param rDist The rate across sites distribution to use.
     * @param checkRooted Tell if we have to check for the tree to be unrooted.
     * If true, any rooted tree will be unrooted before likelihood computation.
     * @param verbose Should I display some info?
     * @throw Exception in an error occured.
     */
    NNIHomogeneousTreeLikelihood(
      const Tree & tree,
      SubstitutionModel * model,
      DiscreteDistribution * rDist,
      bool checkRooted = true,
      bool verbose = true)
      throw (Exception);
  
    /**
     * @brief Build a new NNIHomogeneousTreeLikelihood object.
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
    NNIHomogeneousTreeLikelihood(
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
    NNIHomogeneousTreeLikelihood(const NNIHomogeneousTreeLikelihood & lik);
    
    NNIHomogeneousTreeLikelihood & operator=(const NNIHomogeneousTreeLikelihood & lik);

    virtual ~NNIHomogeneousTreeLikelihood();

#ifndef NO_VIRTUAL_COV
    NNIHomogeneousTreeLikelihood*
#else
    Clonable*
#endif
    clone() const { return new NNIHomogeneousTreeLikelihood(*this); }

  public:

    void setData(const SiteContainer & sites) throw (Exception)
    {
      DRHomogeneousTreeLikelihood::setData(sites);
      if(_brLikFunction) delete _brLikFunction;
      _brLikFunction = new BranchLikelihood(_likelihoodData->getWeights());
    }

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
    
};

} //end of namespace bpp.

#endif  //_NNIHOMOGENEOUSTREELIKELIHOOD_H_

