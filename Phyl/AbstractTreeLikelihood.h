//
// File: AbstractTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Oct 17 17:57:21 2003
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

#ifndef _ABSTRACTTREELIKELIHOOD_H_
#define _ABSTRACTTREELIKELIHOOD_H_

#include "TreeLikelihood.h"
#include "Tree.h"
#include "TreeTemplate.h"

//From SeqLib:
#include <Seq/SiteContainer.h>

/**
 * @brief Low-level implementatoin of the TreeLikelihood interface. 
 *
 * This class implements a few methods useful for most of likelihood
 * computation methods.
 * 
 * - The Parametrizable interface;
 * - The getTree() method;
 * 
 * It also adds an abstract method for recursive computations.
 */
class AbstractTreeLikelihood :
	public virtual TreeLikelihood,
	public virtual AbstractParametrizable
{
	protected:
		const Alphabet * _alphabet;
		const SiteContainer * _data;
		mutable TreeTemplate<Node> *  _tree;
		bool _computeDerivatives;

	public:
		AbstractTreeLikelihood() {}
		virtual ~AbstractTreeLikelihood() {}
	
	public:
		/**
		 * @name The TreeLikelihood interface.
		 *
		 * @{
		 */
		const SiteContainer * getData() const { return _data; }
		const Alphabet * getAlphabet() const { return _data -> getAlphabet(); }	
		Vdouble getLikelihoodForEachSite()                 const;
		Vdouble getLogLikelihoodForEachSite()              const;
		VVdouble getLikelihoodForEachSiteForEachState()    const;
		VVdouble getLogLikelihoodForEachSiteForEachState() const;
		unsigned int getNumberOfSites() const { return _data -> getNumberOfSites(); }
		unsigned int getNumberOfStates() const { return _data -> getAlphabet() -> getSize(); }
		Tree * getTree() const { return _tree; }
		void setComputeDerivatives(bool yn) { _computeDerivatives = yn; }
		bool computeDerivatives() const { return _computeDerivatives; }
		/** @} */

	protected:
		
		/**
		 * @brief Recompute _pxy, _dpxy and _d2pxy arrays, and derivatives if needed.
		 *
		 * This method is called when some parameter has changed.
		 *
		 * @param params The parameters that changed.
		 */
		virtual void fireParameterChanged(const ParameterList & params) = 0;
		
};


#endif	//_ABSTRACTTREELIKELIHOOD_H_

