//
// File: AbstractTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Oct 17 17:57:21 2003
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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
#include "SubstitutionProcess.h"
#include "../Tree.h"
#include "../TreeTemplate.h"

#include <Bpp/Numeric/AbstractParametrizable.h>

//From SeqLib:
#include <Bpp/Seq/Container/SiteContainer.h>

namespace bpp
{
namespace newlik //To avoid name conflick with old likelihood
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
    std::auto_ptr<const SiteContainer> data_;
    std::auto_ptr<SubstitutionProcess> process_;
		bool computeFirstOrderDerivatives_;
		bool computeSecondOrderDerivatives_;
    bool initialized_;
    bool verbose_;
    size_t nbSites_;
    size_t nbDistinctSites_; //Cross-check with SitePartition ! TODO
    size_t nbStates_;
    size_t nbClasses_;


	public:
		AbstractTreeLikelihood():
      AbstractParametrizable(""),
      data_(0),
      process_(0),
      computeFirstOrderDerivatives_(true),
      computeSecondOrderDerivatives_(true),
      initialized_(false),
      verbose_(true),
      nbSites_(0),
      nbDistinctSites_(0),
      nbStates_(0),
      nbClasses_(0)
    {}

    AbstractTreeLikelihood(
        const SiteContainer* data,
        SubstitutionProcess* process,
        bool verbose = true):
      AbstractParametrizable(""),
      data_(data),
      process_(process),
      computeFirstOrderDerivatives_(true),
      computeSecondOrderDerivatives_(true),
      initialized_(false),
      verbose_(verbose),
      nbSites_(data->getNumberOfSites()),
      nbDistinctSites_(0),
      nbStates_(process->getNumberOfStates()),
      nbClasses_(process->getNumberOfClasses())
    {}

    AbstractTreeLikelihood(const AbstractTreeLikelihood& lik):
      AbstractParametrizable(lik),
      data_(0),
      process_(0),
      computeFirstOrderDerivatives_(lik.computeFirstOrderDerivatives_),
      computeSecondOrderDerivatives_(lik.computeSecondOrderDerivatives_),
      initialized_(lik.initialized_), 
      verbose_(lik.verbose_),
      nbSites_(lik.nbSites_),
      nbDistinctSites_(lik.nbDistinctSites_),
      nbStates_(lik.nbStates_),
      nbClasses_(lik.nbClasses_)
    {
      if (lik.data_.get()) data_.reset(lik.data_->clone());
      if (lik.process_.get()) process_.reset(lik.process_->clone());
    }

    AbstractTreeLikelihood& operator=(const AbstractTreeLikelihood& lik)
    {
      AbstractParametrizable::operator=(lik);
      if (lik.data_.get())    data_.reset(lik.data_->clone());
      else                    data_.reset();
      if (lik.process_.get()) process_.reset(lik.process_->clone());
      else                    process_.reset();
      computeFirstOrderDerivatives_  = lik.computeFirstOrderDerivatives_;
      computeSecondOrderDerivatives_ = lik.computeSecondOrderDerivatives_;
      initialized_                   = lik.initialized_;
      verbose_                       = lik.verbose_;
      nbSites_                       = lik.nbSites_;
      nbDistinctSites_               = lik.nbDistinctSites_;
      nbStates_                      = lik.nbStates_;
      nbClasses_                     = lik.nbClasses_;
      return *this;
    }

    /**
     * @brief Abstract class destructor
     *
     * This destructor is empty.
     */
		virtual ~AbstractTreeLikelihood()
    {
      //Auto pointers take care of everything!
    }
	
	public:
		/**
		 * @name The TreeLikelihood interface.
		 *
		 * @{
		 */
		const SiteContainer* getData() const { return data_.get(); }
		const Alphabet* getAlphabet() const { return data_->getAlphabet(); }	
		
    size_t getNumberOfSites() const { return data_->getNumberOfSites(); }
		size_t getNumberOfStates() const { return data_->getAlphabet()->getSize(); }
		size_t getNumberOfClasses() const { return process_->getNumberOfClasses(); }

    Vdouble getLikelihoodForEachSite() const;
		VVdouble getLikelihoodForEachSiteForEachState() const;
		VVdouble getLikelihoodForEachSiteForEachClass() const;
		VVVdouble getLikelihoodForEachSiteForEachClassForEachState() const;
		
    VVdouble getPosteriorProbabilitiesOfEachClass() const;
    std::vector<size_t> getClassWithMaxPostProbOfEachSite() const;
    
		const Tree& getTree() const { return process_->getTree(); }
		void enableDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
		void enableFirstOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = yn; }
		void enableSecondOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
		bool enableFirstOrderDerivatives() const { return computeFirstOrderDerivatives_; }
		bool enableSecondOrderDerivatives() const { return computeSecondOrderDerivatives_; }
    bool isInitialized() const { return initialized_; }
    void initialize() throw (Exception) { initialized_ = true; }
		
    ParameterList getSubstitutionProcessParameters() const { return process_->getParameters(); }
    ParameterList getBranchLengthsParameters() const { return process_->getParametrizableTree().getParameters(); }
    ParameterList getTransitionProbabilitiesParameters() const { return process_->getTransitionProbabilitiesParameters(); }
    //TODO: this has to be modified to deal with special cases...
    ParameterList getDerivableParameters() const { return getBranchLengthsParameters(); }
    ParameterList getNonDerivableParameters() const { return getSubstitutionProcessParameters(); }    
		/** @} */

    //TODO jdutheil on 08.04.13 we drop model iterators for now
    //ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const {
    //  return new ConstNoPartitionBranchModelIterator(model_.get(), sitePartition_->getNumberOfPatternsForPartition(0));
    //}

    //ConstSiteModelIterator* getNewSiteModelIterator(size_t siteIndex) const {
    //  return new ConstHomogeneousSiteModelIterator(*pTree_, model_.get());
    //}

  protected:
		/**
     * @name Generic tools to deal with likelihood arrays
     *
     * @{
     */
    
    /**
     * @brief Set all conditional likelihoods to 1.
     *
     * @param likelihoodArray the likelihood array.
     */
    static void resetLikelihoodArray(VVVdouble& likelihoodArray);

    /**
     * @brief Print the likelihood array to terminal (debugging tool).
     * 
     * @param likelihoodArray the likelihood array.
     */
		static void displayLikelihoodArray(const VVVdouble& likelihoodArray);

};

} //end of namespace newlik.
} //end of namespace bpp.

#endif	//_ABSTRACTTREELIKELIHOOD_H_

