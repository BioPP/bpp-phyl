//
// File: AbstractTreeLikelihood.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Fri Oct 17 17:57:21 2003
//

#ifndef _ABSTRACTTREELIKELIHOOD_H_
#define _ABSTRACTTREELIKELIHOOD_H_

#include "TreeLikelihood.h"

//From SeqLib:
#include <Seq/SiteContainer.h>

/**
 * @brief <p>This class implements a few methods useful for most of likelihood
 * computation methods.</p>
 * <ul><li>The Parametrizable interface;</li>
 *     <li>The getTree() method;</li></ul>
 * </p>It also adds an abstract method for recursive computations.</p>
 */
class AbstractTreeLikelihood : public TreeLikelihood
{
	protected:
		const SiteContainer * _data;
		mutable ParameterList _parameters;
		mutable Tree *        _tree;

	public:
		virtual ~AbstractTreeLikelihood();
	
	public:
		void computeTreeLikelihood();
		const SiteContainer * getData() const;
	
		/**
		 * @name The Parametrizable interface.
		 *
		 * @{
		 */
		ParameterList getParameters() const throw (Exception);
	
		double getParameter(const string & name) const
			throw (ParameterNotFoundException);

		void setAllParametersValues(const ParameterList & params)
			throw (ParameterNotFoundException, ConstraintException);
	
		void setParameterValue(const string & name, double value)
			throw (ParameterNotFoundException, ConstraintException);
	
		void setParametersValues(const ParameterList & params)
			throw (ParameterNotFoundException, ConstraintException);
	
		void matchParametersValues(const ParameterList & params)
			throw (ConstraintException);
			
		/** @} */
	
		Vdouble getLikelihoodForEachSite()    const;
		Vdouble getLogLikelihoodForEachSite() const;

		Tree * getTree() const;
	
	protected:
		
		/**
		 * @brief Compute the likelihood for a subtree defined by the Tree::Node <i>node</i>.
		 *
		 * @param node The root of the subtree.
		 */
		virtual void computeSubtreeLikelihood(Node * node) = 0; //Recursive method.			

		/**
		 * @brief Recompute _pxy, _dpxy and _d2pxy arrays.
		 *
		 * This method is called when some parameter has changed.
		 *
		 * @param params The parameters that changed.
		 */
		virtual void fireParameterChanged(const ParameterList & params) = 0;

};


#endif	//_ABSTRACTTREELIKELIHOOD_H_
