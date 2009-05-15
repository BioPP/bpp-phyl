//
// File: DistanceEstimation.h
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Wed jun 08 10:39 2005
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

#ifndef _DISTANCEESTIMATION_H_
#define _DISTANCEESTIMATION_H_

#include "SubstitutionModel.h"
#include "AbstractTreeLikelihood.h"
#include "DRHomogeneousTreeLikelihood.h"
#include "PseudoNewtonOptimizer.h"

// From Utils:
#include <Utils/Clonable.h>

// From NumCalc:
#include <NumCalc/ParameterList.h>
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/Optimizer.h>
#include <NumCalc/SimpleMultiDimensions.h>
#include <NumCalc/MetaOptimizer.h>

// From SeqLib:
#include <Seq/SiteContainer.h>

namespace bpp
{

/**
 * @brief This class is a simplified version of DRHomogeneousTreeLikelihood for 2-Trees.
 */
class TwoTreeLikelihood:
  public AbstractDiscreteRatesAcrossSitesTreeLikelihood  
{
	protected:
		SiteContainer * _shrunkData;
		vector<string> _seqnames;
		SubstitutionModel * _model;
		ParameterList _brLenParameters;
		
		mutable VVVdouble _pxy;

		mutable VVVdouble _dpxy;

		mutable VVVdouble _d2pxy;

		/**
		 * @brief As previous, but for the global container.
		 *
		 * The size of this vector is equal to the number of sites in the container,
		 * each element corresponds to a site in the container and points to the
		 * corresponding column in the likelihood array of the root node.
		 * If the container contains no repeated site, there will be a strict
		 * equivalence between each site and the likelihood array of the root node.
		 * However, if this is not the case, some pointers may point toward the same
		 * element in the likelihood array.
		 */
		vector<unsigned int> _rootPatternLinks;

		/**
		 * @brief The frequency of each site.
		 */
		vector<unsigned int> _rootWeights;

		//some values we'll need:
		unsigned int _nbSites,         //the number of sites in the container
		             _nbClasses,       //the number of rate classes
		             _nbStates,        //the number of states in the alphabet
		             _nbDistinctSites; //the number of distinct sites in the container

		double _brLen;
		mutable VVVdouble _rootLikelihoods;
		mutable VVdouble _rootLikelihoodsS;
		mutable Vdouble _rootLikelihoodsSR;
		mutable Vdouble _dLikelihoods;
		mutable Vdouble _d2Likelihoods;
		mutable VVdouble _leafLikelihoods1, _leafLikelihoods2;
	
    double _minimumBrLen;
    Constraint * _brLenConstraint;

	public:
		TwoTreeLikelihood(
			const string & seq1, const string & seq2,	
			const SiteContainer & data,
			SubstitutionModel * model,
			DiscreteDistribution * rDist,
			bool verbose)	throw (Exception);

    TwoTreeLikelihood(const TwoTreeLikelihood & lik);
    
    TwoTreeLikelihood & operator=(const TwoTreeLikelihood & lik);

#ifndef NO_VIRTUAL_COV
    TwoTreeLikelihood*
#else
    Clonable*
#endif
    clone() const { return new TwoTreeLikelihood(*this); } 

		virtual ~TwoTreeLikelihood();

	public:

		/**
		 * @name The TreeLikelihood interface.
		 *
		 * Other methods are implemented in the AbstractTreeLikelihood class.
		 *
		 * @{
		 */
		double getLikelihood () const;
		double getLogLikelihood() const;
		double getLikelihoodForASite (unsigned int site) const;
		double getLogLikelihoodForASite(unsigned int site) const;
		ParameterList getBranchLengthsParameters() const;
		ParameterList getSubstitutionModelParameters() const;
    SubstitutionModel * getSubstitutionModelForNode(int nodeId) throw (NodeNotFoundException) { return _model; }
    const SubstitutionModel * getSubstitutionModelForNode(int nodeId) const throw (NodeNotFoundException) { return _model; }
    vector<double> getRootFrequencies() const { return _model->getFrequencies(); }
    /**
     * @brief This method is not applicable for this object.
     */
    const VVVdouble & getTransitionProbabilitiesForNode(const Node* node) const { return _pxy; }
    void setData(const SiteContainer & sites) throw (Exception) {}
    void initialize() throw(Exception);
		/** @} */

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
		 * @brief Get the substitution model used for the computation.
		 *
		 * @return A const pointer toward the substitution model of this instance.
		 */
		virtual const SubstitutionModel * getSubstitutionModel() const { return _model; }
		
		/**
		 * @brief Get the substitution model used for the computation.
		 *
		 * @return A pointer toward the substitution model of this instance.
		 */
		virtual SubstitutionModel * getSubstitutionModel() { return _model; }

		/**
		 * @brief Implements the Function interface.
		 *
		 * Update the parameter list and call the applyParameters() method.
		 * Then compute the likelihoods at each node (computeLikelihood() method)
		 * and call the getLogLikelihood() method.
		 *
		 * If a subset of the whole parameter list is passed to the function,
		 * only these parameters are updated and the other remain constant (i.e.
		 * equal to their last value).
		 *
		 * @param parameters The parameter list to pass to the function.
		 */
		void setParameters(const ParameterList & parameters) throw (ParameterNotFoundException, ConstraintException);
		double getValue() const throw(Exception);
		
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
		 * @brief This method initializes the leaves according to a sequence container.
		 *
		 * Here the container _shrunkData is used.
		 * Likelihood is set to 1 for the state corresponding to the sequence site,
		 * otherwise it is set to 0.
		 *
		 * The two likelihood arrays are initialized according to alphabet
		 * size and sequences length, and filled with 1.
		 *
		 * NB: This method is recursive.
		 *
		 * @param sequences The sequence container to use.
		 */
		virtual void initTreeLikelihoods(const SequenceContainer & sequences) throw (Exception);

		void fireParameterChanged(const ParameterList & params);
		virtual void computeTreeLikelihood();
		virtual void computeTreeDLikelihood();
		virtual void computeTreeD2Likelihood();
		/**
		 * @brief This builds the <i>parameters</i> list from all parametrizable objects,
		 * <i>i.e.</i> substitution model, rate distribution and tree.
		 */
		virtual void initParameters();

		/**
		 * @brief All parameters are stores in a parameter list.
		 *
		 * This function apply these parameters to the substitution model,
		 * to the rate distribution and to the branch lengths.
		 */
		virtual void applyParameters() throw (Exception);	

};

/**
 * @brief Estimate a distance matrix from sequence data, according to a given model.
 *
 * By default, the parameters of the model are fixed to there given values.
 * It is possible to estimate one or several parameters by setting them with the
 * setAdditionalParameters() method.
 * Parameters will be estimated separately for each pair of sequence.
 *
 * For now it is not possible to retrieve estimated values.
 * You'll have to specify a 'profiler' to the optimizer and then look at the file
 * if you want to do so.
 */
class DistanceEstimation:
  public virtual Clonable
{
	protected:
		SubstitutionModel * _model;
		DiscreteDistribution * _rateDist;
		const SiteContainer * _sites;
		DistanceMatrix * _dist;
		Optimizer * _optimizer;
		MetaOptimizer * _defaultOptimizer;
		unsigned int _verbose;
		ParameterList parameters_;

	public:
		
		/**
		 * @brief Create a new DistanceEstimation object and compute distances
		 * according to a given substitution model and a rate distribution.
		 *
		 * @param model    The substitution model to use.
		 * @param rateDist The discrete rate distribution to use.
		 * @param sites    The sequence data.
		 * @param verbose  The verbose level:
		 *  - 0=Off,
		 *  - 1=one * by row computation
		 *  - 2=one * by row computation and one . by column computation
		 *  - 3=2 + optimization verbose enabled
		 *  - 4=3 + likelihood object verbose enabled
		 *  @param computeMat if true the computeMatrix() method is called.
		 */
		DistanceEstimation(SubstitutionModel * model, DiscreteDistribution * rateDist, const SiteContainer * sites, unsigned int verbose = 1, bool computeMat = true):
      _model(model), _rateDist(rateDist), _sites(sites), _dist(NULL), _verbose(verbose)
    {
	    _init();
      if(computeMat) computeMatrix();
    }
		
    /**
     * @brief Copy constructor.
     *
     * Only the distance matrix is hard-copied, if there is one.
     *
     * @param distanceEstimation The object to copy.
     */
    DistanceEstimation(const DistanceEstimation & distanceEstimation):
      _model(distanceEstimation._model),
      _rateDist(distanceEstimation._rateDist),
      _sites(distanceEstimation._sites),
      _dist(NULL),
      _optimizer(dynamic_cast<Optimizer *>(distanceEstimation._optimizer->clone())),
      _defaultOptimizer(dynamic_cast<MetaOptimizer *>(_defaultOptimizer->clone())),
      _verbose(distanceEstimation._verbose),
      parameters_(distanceEstimation.parameters_)
    {
      if(distanceEstimation._dist != NULL)
        _dist = new DistanceMatrix(*distanceEstimation._dist);
      else
        _dist = NULL;
    }

    /**
     * @brief Assigment operator.
     *
     * Only the distance matrix is hard-copied, if there is one.
     * 
     * @param distanceEstimation The object to copy.
     * @return A reference toward this object.
     */
    DistanceEstimation & operator=(const DistanceEstimation & distanceEstimation)
    {
      _model = distanceEstimation._model;
      _rateDist = distanceEstimation._rateDist;
      _sites = distanceEstimation._sites;
      if(distanceEstimation._dist != NULL)
        _dist = new DistanceMatrix(*distanceEstimation._dist);
      else
        _dist = NULL;
      _optimizer = dynamic_cast<Optimizer *>(distanceEstimation._optimizer->clone());
      // _defaultOptimizer has already been initialized since the default constructor has been called.
      _verbose = distanceEstimation._verbose;
      parameters_ = distanceEstimation.parameters_;
      return *this;
    }

		virtual ~DistanceEstimation()
		{
			if(_dist != NULL) delete _dist;
			delete _defaultOptimizer;
      delete _optimizer;
		}

#ifndef NO_VIRTUAL_COV
    DistanceEstimation *
#else
    Clonable *
#endif
    clone() const { return new DistanceEstimation(*this); }
		
  private:
    void _init()
    {
      MetaOptimizerInfos* desc = new MetaOptimizerInfos();
      vector<string> name;
      name.push_back("BrLen");
      desc->addOptimizer("Branch length", new PseudoNewtonOptimizer(NULL), name, 2, MetaOptimizerInfos::IT_TYPE_FULL);
      ParameterList tmp = _model->getParameters();
      tmp.addParameters(_rateDist->getParameters());
      desc->addOptimizer("substitution model and rate distribution", new SimpleMultiDimensions(NULL), tmp.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);
    	_defaultOptimizer = new MetaOptimizer(NULL, desc);
      _defaultOptimizer->setMessageHandler(NULL);
	    _defaultOptimizer->setProfiler(NULL);
      _defaultOptimizer->getStopCondition()->setTolerance(0.0001);
	    _optimizer = dynamic_cast<Optimizer *>(_defaultOptimizer->clone());
    }

	public:

		/**
		 * @brief Perform the distance computation.
		 *
		 * Result can be called by the getMatrix() method.
		 *
		 * @throw NullPointerException if at least one of the model,
		 * rate distribution or data are not initialized.
		 */
		void computeMatrix() throw (NullPointerException);
		
		/**
		 * @brief Get the distance matrix.
		 *
		 * @return A pointer toward the computed distance matrix.
		 */
		DistanceMatrix * getMatrix() const { return _dist == NULL ? NULL : new DistanceMatrix(* _dist); }

		SubstitutionModel * getModel() const { return _model; }
		void resetModel() { _model = NULL; }

		DiscreteDistribution * getRateDistribution() const { return _rateDist; }
		void resetRateDistribution() { _rateDist = NULL; }

		void setData(const SiteContainer * sites) { _sites = sites; }
		const SiteContainer * getData() const { return _sites; }
		void resetData() { _sites = NULL; }
		
		void setOptimizer(const Optimizer * optimizer)
    { 
      if(_optimizer) delete _optimizer;
      _optimizer = dynamic_cast<Optimizer *>(optimizer->clone());
    }
		const Optimizer * getOptimizer() const { return _optimizer; }
		Optimizer * getOptimizer() { return _optimizer; }
		void resetOptimizer() { _optimizer = dynamic_cast<Optimizer *>(_defaultOptimizer->clone()); }

		/**
		 * @brief Specify a list of parameters to be estimated.
		 *
		 * Parameters will be estimated separately for each distance.
		 *
		 * @param parameters A list of parameters to estimate.
		 */
		void setAdditionalParameters(const ParameterList & parameters)
		{
      parameters_ = parameters;
		}

    /**
     * @brief Reset all additional parameters.
     */
    void resetAdditionalParameters()
    {
      parameters_.reset();
    }

    /**
     * @param verbose Verbose level.
     */
    void setVerbose(unsigned int verbose) { _verbose = verbose; }
    /**
     * @return Verbose level.
     */
    unsigned int getVerbose() const { return _verbose; }
};

} //end of namespace bpp.

#endif //_DISTANCEESTIMATION_H_

