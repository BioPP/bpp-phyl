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

//From NumCalc:
#include <NumCalc/AbstractParametrizable.h>

//From SeqLib:
#include <Seq/SiteContainer.h>

namespace bpp
{

/**
 * @brief Partial implementation of the TreeLikelihood interface. 
 *
 * This class implements a few methods useful for most of likelihood
 * computation methods.
 *
 * It includes a tree_ and a data_ pointers.
 * This objects are owned by the class, and hence hard copied when cloning, and destroyed by the destructor.
 * 
 * - The Parametrizable interface;
 * - The getTree() method;
 * 
 * It also adds an abstract method for recursive computations.
 */
class AbstractTreeLikelihood :
	public virtual TreeLikelihood,
	public AbstractParametrizable
{
	protected:
		const SiteContainer* data_;
		mutable TreeTemplate<Node>* tree_;
		bool computeFirstOrderDerivatives_;
		bool computeSecondOrderDerivatives_;
    bool initialized_;

	public:
		AbstractTreeLikelihood():
      AbstractParametrizable(""),
      data_(0),
      tree_(0),
      computeFirstOrderDerivatives_(true),
      computeSecondOrderDerivatives_(true),
      initialized_(false) {}

    AbstractTreeLikelihood(const AbstractTreeLikelihood & lik):
      AbstractParametrizable(lik),
      data_(0),
      tree_(0),
      computeFirstOrderDerivatives_(lik.computeFirstOrderDerivatives_),
      computeSecondOrderDerivatives_(lik.computeSecondOrderDerivatives_),
      initialized_(lik.initialized_) 
    {
      if(lik.data_) data_ = dynamic_cast<SiteContainer *>(lik.data_->clone());
      if(lik.tree_) tree_ = lik.tree_->clone();
    }

    AbstractTreeLikelihood & operator=(const AbstractTreeLikelihood & lik)
    {
      AbstractParametrizable::operator=(lik);
      if(data_) delete data_;
      if(lik.data_) data_ = dynamic_cast<SiteContainer *>(lik.data_->clone());
      else          data_ = 0;
      if(tree_) delete tree_;
      if(lik.tree_) tree_ = lik.tree_->clone();
      else          tree_ = 0;
      computeFirstOrderDerivatives_ = lik.computeFirstOrderDerivatives_;
      computeSecondOrderDerivatives_ = lik.computeSecondOrderDerivatives_;
      initialized_        = lik.initialized_;
      return *this;
    }

    /**
     * @brief Abstract class destructor
     *
     * This destructor is empty.
     */
		virtual ~AbstractTreeLikelihood()
    {
      if(data_) delete data_;
      if(tree_) delete tree_;
    }
	
	public:
		/**
		 * @name The TreeLikelihood interface.
		 *
		 * @{
		 */
		const SiteContainer * getData() const { return data_; }
		const Alphabet * getAlphabet() const { return data_->getAlphabet(); }	
		Vdouble getLikelihoodForEachSite()                 const;
		Vdouble getLogLikelihoodForEachSite()              const;
		VVdouble getLikelihoodForEachSiteForEachState()    const;
		VVdouble getLogLikelihoodForEachSiteForEachState() const;
		unsigned int getNumberOfSites() const { return data_->getNumberOfSites(); }
		unsigned int getNumberOfStates() const { return data_->getAlphabet()->getSize(); }
		const Tree * getTree() const { return tree_; }
		void enableDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
		void enableFirstOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = yn; }
		void enableSecondOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
		bool enableFirstOrderDerivatives() const { return computeFirstOrderDerivatives_; }
		bool enableSecondOrderDerivatives() const { return computeSecondOrderDerivatives_; }
    bool isInitialized() const { return initialized_; }
    void initialize() throw (Exception) { initialized_ = true; }
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

} //end of namespace bpp.

#endif	//_ABSTRACTTREELIKELIHOOD_H_

