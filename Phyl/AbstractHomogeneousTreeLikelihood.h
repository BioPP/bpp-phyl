//
// File: AbstractHomogeneousTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Thr Dec 23 12:03 2004
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

#ifndef _ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H_
#define _ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H_

#include "AbstractDiscreteRatesAcrossSitesTreeLikelihood.h"
#include "SubstitutionModel.h"

/**
 * @brief Partial implementation for homogeneous model of the TreeLikelihood interface.
 *
 * This class provides a pointer toward a single substitution model + several utilitary variables.
 */
class AbstractHomogeneousTreeLikelihood:
  public AbstractDiscreteRatesAcrossSitesTreeLikelihood
{
	protected:
		SubstitutionModel * _model;
		ParameterList _brLenParameters;
		
		mutable map<int, VVVdouble> _pxy;

		mutable map<int, VVVdouble> _dpxy;

		mutable map<int, VVVdouble> _d2pxy;
				
		/**
		 * @brief Pointer toward all nodes in the tree.
		 *
		 * The order of the nodes in the vector if the order of the named branches.
		 */
		vector<Node *> _nodes;

		//some values we'll need:
		unsigned int _nbSites,         //the number of sites in the container
                 _nbDistinctSites, //the number of distinct sites
		             _nbClasses,       //the number of rate classes
		             _nbStates,        //the number of states in the alphabet
		             _nbNodes;         //the number of nodes in the tree

    bool _verbose;

    double _minimumBrLen;
    Constraint * _brLenConstraint;

	public:
		AbstractHomogeneousTreeLikelihood(
			const Tree & tree,
			SubstitutionModel * model,
			DiscreteDistribution * rDist,
      bool checkRooted = true,
			bool verbose = true)
      throw (Exception);

    /**
     * @brief Copy constructor
     *
     * This constructor is to be called by the derived class copy constructor.
     */
    AbstractHomogeneousTreeLikelihood(const AbstractHomogeneousTreeLikelihood & lik);
    
    /**
     * @brief Assignation operator
     *
     * This operator is to be called by the derived class operator.
     */
    AbstractHomogeneousTreeLikelihood & operator=(const AbstractHomogeneousTreeLikelihood & lik);
 
		virtual ~AbstractHomogeneousTreeLikelihood();
		
	public:

    /**
     * @brief Method called by constructor.
     */
    void init(const Tree & tree,
			SubstitutionModel * model,
			DiscreteDistribution * rDist,
      bool checkRooted,
			bool verbose) throw (Exception);

		
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
		
	public: //Specific methods:

		/**
		 * @brief This builds the <i>parameters</i> list from all parametrizable objects,
		 * <i>i.e.</i> substitution model, rate distribution and tree.
		 */
		virtual void initParameters();

		/**
		 * @brief This removes a particular parameter from the list.
		 *
		 * This method may be used to not estimate a parameter after having
		 * fixed its value. The previous method reset all calls of thos one.
		 *
		 * @param name The name of the parameter to ignore.
		 */
		virtual void ignoreParameter(const string & name) throw (ParameterNotFoundException);

		/**
		 * @brief All parameters are stored in a parameter list.
		 * This function apply these parameters to the substitution model,
		 * to the rate distribution and to the branch lengths.
		 */
		virtual void applyParameters() throw (Exception);	

		virtual void initBranchLengthsParameters();

    void setMinimumBranchLength(double minimum)
    {
      _minimumBrLen = minimum;
      if(_brLenConstraint != NULL) delete _brLenConstraint;
      _brLenConstraint = new IncludingPositiveReal(_minimumBrLen);
      initBranchLengthsParameters();
    }

    double getMinimumBranchLength() const { return _minimumBrLen; }

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

#endif //_ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H_

