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

#ifndef _TREELIKELIHOODCALCULATION_H_
#define _TREELIKELIHOODCALCULATION_H_

#include "TreeLikelihoodData.h"
#include "SubstitutionProcess.h"

#include <cstddef>

#include <Bpp/Clonable.h>

namespace bpp
{
namespace newlik
{

/**
 * @brief Exception thrown in case a TreeLikelihoodCalculation object was not properly initialized.
 */
class TreeLikelihoodCalculationNotInitializedException:
  public virtual Exception
{
public:
  TreeLikelihoodCalculationNotInitializedException(const std::string& msg):
    Exception("TreeLikelihoodCalculation not initialized. " + msg) {}

};

/**
 * @brief The TreeLikelihoodCalculation interface.
 */
class TreeLikelihoodCalculation:
  public virtual Clonable
{

public:
  virtual ~TreeLikelihoodCalculation() {}

  virtual TreeLikelihoodCalculation* clone() const = 0;

public:

  /**
   * @return The alphabet of the data set for which this object is initialized.
   * @throw TreeLikelihoodCalculationNotInitializedException If this instance was not initialized.
   */
  virtual const Alphabet* getAlphabet() const throw (TreeLikelihoodCalculationNotInitializedException) = 0;

  /**
   * @return The process used for calculation.
   *
   */

  virtual const SubstitutionProcess* getSubstitutionProcess() const = 0;

  /**
   * @return The size of the data set for which this object is initialized.
   * @throw TreeLikelihoodCalculationNotInitializedException If this instance was not initialized.
   */
  virtual size_t getNumberOfSites() const throw (TreeLikelihoodCalculationNotInitializedException) = 0;

  /**
   * @brief Get the pattern index for a given site position in the original data.
   *
   * @return The pattern index in the associated data set, given the position in the original data set used to initialize this instance.
   * @throw TreeLikelihoodCalculationNotInitializedException If this instance was not initialized.
   * @throw IndexOutOfBoundsException If the input position is invalid.
   */
  virtual size_t getSiteIndex(size_t site) const throw (TreeLikelihoodCalculationNotInitializedException, IndexOutOfBoundsException) = 0;

  /**
   * @brief Tell is this instance is properly insitialized, that is, if the setData method has been called once.
   *
   * @return True if a data set is associated to this instance.
   */
  virtual bool isInitialized() const  = 0;

  /**
   * @brief Initialize the object according to a data set.
   *
   * @param sites A sequence alignment to initialize the object with.
   */
  virtual void setData(const SiteContainer& sites) throw (Exception) = 0;
  
  /**
   * @return The data set used to initialize this object.
   * @throw TreeLikelihoodCalculationNotInitializedException In this instance was not initialized.
   */
  virtual const SiteContainer* getData() const throw (TreeLikelihoodCalculationNotInitializedException) = 0;

  virtual TreeLikelihoodData* getLikelihoodData() = 0;
  
  virtual const TreeLikelihoodData* getLikelihoodData() const = 0;
  
  /**
   * @brief Get the log-likelihood for the data set.
   *
   * @return The log-likelihood for the data set.
   */
  virtual double getLogLikelihood() const = 0;
  
  /**
   * @brief Get the likelihood for a site.
   *
   * @param site The site index to analyse.
   * @return The likelihood for site <i>site</i>.
   */
  virtual double getLikelihoodForASite(size_t site) const = 0;

  /**
   * @brief Get the likelihood for a site and for a state.
   *
   * @param site The site index to analyse.
   * @param state The state to consider.
   * @return The likelihood for site <i>site</i> and state <i>state</i>.
   */
  virtual double getLikelihoodForASiteForAState(size_t site, int state) const = 0;

  /**
   * @brief Get the logarithm of the likelihood for a site knowing its model class.
   *
   * @param site       The site index.
   * @param classIndex The model class index.
   * @return The likelihood for the specified site and model class.
   */
  virtual double getLikelihoodForASiteForAClass(size_t site, size_t classIndex) const = 0;

  /**
   * @brief Get the likelihood for a site knowing its model class and its ancestral state.
   *
   * @param site       The site index.
   * @param classIndex The model class index.
   * @param state      The ancestral state.
   * @return The likelihood for the specified site and model class and ancestral state..
   */
  virtual double getLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state) const = 0;

  /**
   * @brief Get the derivate of log-likelihood for the data set.
   *
   * Derivation is performed according to a given branch length,
   * which was specified by the last call to computeDLikelihood.
   *
   * @return The derivate of log-likelihood for the data set.
   */
  virtual double getDLogLikelihood() const = 0;
  
  /**
   * @brief Get the second order derivate of log-likelihood for the data set.
   *
   * Derivation is performed according to a given branch length,
   * which was specified by the last call to computeD2Likelihood.
   *
   * @return The second order derivate of log-likelihood for the data set.
   */
  virtual double getD2LogLikelihood() const = 0;

  /**
   * @brief Get the derivative of the likelihood for a site.
   *
   * @param site The site index to analyse.
   * @return The derivative of likelihood for site <i>site</i>.
   */
  virtual double getDLikelihoodForASite(size_t site) const = 0;
  
  /**
   * @brief Get the second-order derivative of the likelihood for a site.
   *
   * @param site The site index to analyse.
   * @return The second-order derivative of likelihood for site <i>site</i>.
   */
  virtual double getD2LikelihoodForASite(size_t site) const = 0;  

  /**
   * @brief Initiate a likelihood computation.
   */
  virtual void computeTreeLikelihood() = 0;
  
  /**
   * @brief Initiate a derivative likelihood computation.
   *
   * @param variable The name of a parameter according to which
   * derivatives can been computed. If variable is not valid, the
   * derivatives set to 0.
   *
   */
  virtual void computeTreeDLikelihood(const std::string& variable) = 0;

  /**
   * @brief Initiate a second-order derivative likelihood computation.
   *
   * @param variable The name of a parameter according to which
   * second-order derivatives can been computed. If variable is not
   * valid, the derivatives set to 0.
   */
  
  virtual void computeTreeD2Likelihood(const std::string& variable) = 0;

};

} // end of namespace newlik.
} // end of namespace bpp.

#endif //_TREELIKELIHOODCALCULATION_H_

