//
// File: AbstractTreeLikelihoodCalculation.h
// Created by: Julien Dutheil
// Created on: Tue July 23 10:50 2013
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

#ifndef _ABSTRACTTREELIKELIHOODCALCULATION_H_
#define _ABSTRACTTREELIKELIHOODCALCULATION_H_

#include "TreeLikelihoodCalculation.h"
#include "SubstitutionProcess.h"

namespace bpp
{
namespace newlik
{

/**
 * @brief Partial implementation of the TreeLikelihoodCalculation interface.
 */
class AbstractTreeLikelihoodCalculation:
  public virtual TreeLikelihoodCalculation
{

protected:
  const SubstitutionProcess* process_;
  std::auto_ptr<const SiteContainer> data_;
  size_t nbSites_;
  size_t nbDistinctSites_;
  size_t nbStates_;
  size_t nbClasses_;
  bool initialized_;
  bool verbose_;

public:
  AbstractTreeLikelihoodCalculation(SubstitutionProcess* process, bool verbose = true):
    process_(process),
    data_(0),
    nbSites_(0),
    nbDistinctSites_(0),
    nbStates_(process->getNumberOfStates()),
    nbClasses_(process->getNumberOfClasses()),
    initialized_(false),
    verbose_(verbose)
  {}
  
  AbstractTreeLikelihoodCalculation(const AbstractTreeLikelihoodCalculation& tlc):
    process_(tlc.process_),
    data_(0),
    nbSites_(tlc.nbSites_),
    nbDistinctSites_(tlc.nbDistinctSites_),
    nbStates_(tlc.nbStates_),
    nbClasses_(tlc.nbClasses_),
    initialized_(tlc.initialized_),
    verbose_(tlc.verbose_)
  {
    if (tlc.data_.get()) data_.reset(tlc.data_->clone());
  }
  
  AbstractTreeLikelihoodCalculation& operator=(const AbstractTreeLikelihoodCalculation& tlc)
  {
    process_ = tlc.process_;
    if (tlc.data_.get()) data_.reset(tlc.data_->clone());
    else data_.reset();
    nbSites_                       = tlc.nbSites_;
    nbDistinctSites_               = tlc.nbDistinctSites_;
    nbStates_                      = tlc.nbStates_;
    nbClasses_                     = tlc.nbClasses_;
    initialized_                   = tlc.initialized_;
    verbose_                       = tlc.verbose_;
    return *this;
  }

  virtual ~AbstractTreeLikelihoodCalculation() {}

public:

  bool isInitialized() const { return initialized_; }
  
  const SubstitutionProcess* getSubstitutionProcess() const { return process_;}

  double getLogLikelihood() const;

  double getDLogLikelihood() const;
  
  double getD2LogLikelihood() const;
  

  /**
   * @brief Print the likelihood array to terminal (debugging tool).
   *
   * @param likelihoodArray the likelihood array.
   */

  static void displayLikelihoodArray(const VVVdouble& likelihoodArray);

};

} // end of namespace newlik.
} // end of namespace bpp.

#endif  // _ABSTRACTTREELIKELIHOODCALCULATION_H_

