//
// File: AbstractHomogeneousTreeLikelihood.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Thr Dec 23 12:03 2004
//

#ifndef _ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H_
#define _ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H_

#include "AbstractTreeLikelihood.h"
#include "DiscreteRatesAcrossSites.h"
#include "SubstitutionModel.h"

class AbstractHomogeneousTreeLikelihood: public AbstractTreeLikelihood, public DiscreteRatesAcrossSites
{
	protected:
		SubstitutionModel * _model;
		DiscreteDistribution * _rateDistribution;
		ParameterList _brLenParameters;
		
		mutable map<const Node *, VVVdouble> _pxy;

		mutable map<const Node *, VVVdouble> _dpxy;

		mutable map<const Node *, VVVdouble> _d2pxy;
		
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
		 * @brief Pointer toward all nodes in the tree.
		 *
		 * This is used for parameters estimation only, not for likelihood computation.
		 * The order of the nodes in the vector if the order of the named branches and
		 * is initalized once for all in the constructor. It then depends of the tree
		 * topology. This may lead to some problems when we'll act on tree topology...
		 */
		 vector<Node *> _nodes;

		//some values we'll need:
		unsigned int _nbSites,         //the number of sites in the container
		             _nbClasses,       //the number of rate classes
		             _nbStates,        //the number of states in the alphabet
		             _nbNodes;         //the number of nodes in the tree


	public:
		AbstractHomogeneousTreeLikelihood(
			Tree & tree,
			const SiteContainer & data,
			SubstitutionModel * model,
			DiscreteDistribution * rDist,
			bool verbose = true)
			throw (Exception);

		virtual ~AbstractHomogeneousTreeLikelihood();
		
	public:
		/**
		 * @name The TreeLikelihood interface.
		 *
		 * Other methods are implemented in the AbstractTreeLikelihood class.
		 *
		 * @{
		 */
		ParameterList getBranchLengthsParameters() const;
		ParameterList getSubstitutionModelParameters() const;
		/** @} */

		/**
		 * @name The DiscreteRatesAcrossSites interface implementation:
		 *
		 * @{
		 */
		const DiscreteDistribution * getRateDistribution() const;
		      DiscreteDistribution * getRateDistribution();
		ParameterList getRateDistributionParameters() const;
		/** @} */

		/**
		 * @brief Get the substitution model used for the computation.
		 *
		 * @return A const pointer toward the substitution model of this instance.
		 */
		virtual const SubstitutionModel * getSubstitutionModel() const;
		
		/**
		 * @brief Get the substitution model used for the computation.
		 *
		 * @return A pointer toward the substitution model of this instance.
		 */
		virtual SubstitutionModel * getSubstitutionModel();
		
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
		 * @brief All parameters are stores in a parameter list.
		 *
		 * This function apply these parameters to the substitution model,
		 * to the rate distribution and to the branch lengths.
		 */
		virtual void applyParameters() throw (Exception);	

		virtual void initBranchLengthsParameters();

		void resetLikelihoodArray(VVVdouble & likelihoodArray);

		static void displayLikelihoodArray(const VVVdouble & likelihoodArray);

};

#endif //_ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H_

