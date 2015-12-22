//
// File: LikelihoodTreeCalculation.h
// Created by: Laurent Guéguen
// Created on: mardi 23 juin 2015, à 14h 01
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

#ifndef _LIKELIHOOD_TREE_CALCULATION_H_
#define _LIKELIHOOD_TREE_CALCULATION_H_

#include "LikelihoodTree.h"
#include "SubstitutionProcess.h"

#include <cstddef>

#include <Bpp/Clonable.h>

namespace bpp
{
/**
 * @brief Exception thrown in case a LikelihoodTreeCalculation object was not properly initialized.
 */
class LikelihoodTreeCalculationNotInitializedException:
  public virtual Exception
{
public:
  LikelihoodTreeCalculationNotInitializedException(const std::string& msg):
    Exception("LikelihoodTreeCalculation not initialized. " + msg) {}

};

/**
 * @brief The LikelihoodTreeCalculation interface.
 */
class LikelihoodTreeCalculation:
  public virtual Clonable
{

public:
  virtual ~LikelihoodTreeCalculation() {}

  virtual LikelihoodTreeCalculation* clone() const = 0;

public:

  /**
   * @return The alphabet of the data set for which this object is initialized.
   * @throw LikelihoodTreeCalculationNotInitializedException If this instance was not initialized.
   */
  virtual const Alphabet* getAlphabet() const throw (LikelihoodTreeCalculationNotInitializedException) = 0;

  /**
   * @return The process used for calculation.
   *
   */

  virtual const SubstitutionProcess* getSubstitutionProcess() const = 0;

  /**
   * @return The size of the data set for which this object is initialized.
   * @throw LikelihoodTreeCalculationNotInitializedException If this instance was not initialized.
   */
  
  virtual size_t getNumberOfSites() const = 0;

  /**
   * @brief Get the pattern index for a given site position in the original data.
   *
   * @return The pattern index in the associated data set, given the position in the original data set used to initialize this instance.
   * @throw LikelihoodTreeCalculationNotInitializedException If this instance was not initialized.
   * @throw IndexOutOfBoundsException If the input position is invalid.
   */
  virtual size_t getSiteIndex(size_t site) const throw (LikelihoodTreeCalculationNotInitializedException, IndexOutOfBoundsException) = 0;

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
  
  virtual void setData(const SiteContainer& sites) = 0;
  
  /**
   * @return The data set used to initialize this object.
   * @throw LikelihoodTreeCalculationNotInitializedException In this instance was not initialized.
   */
  virtual const SiteContainer* getData() const = 0;

  virtual LikelihoodTree& getLikelihoodData() = 0;
  
  virtual const LikelihoodTree& getLikelihoodData() const = 0;

  /**
   * @return if log-likelihood is used at root in a given class
   *
   */
  
  virtual bool usesLogAtRoot(size_t nClass) const = 0;

  /**
   * @brief sets using log in all likelihood arrays.
   *
   */
  
  virtual void setAllUseLog(bool useLog) = 0;

  /**
   * @brief Get the log-likelihood for the data set.
   *
   * @return The log-likelihood for the data set.
   */

  virtual double getLogLikelihood() = 0;
  
  /**
   * @brief Get the likelihood for a site.
   *
   * @param site The site index to analyse.
   * @return The likelihood for site <i>site</i>.
   */

  virtual double getLikelihoodForASite(size_t site) = 0;
      
  virtual double getLogLikelihoodForASite(size_t site) = 0;
  
  double getDLikelihoodForASite(size_t site);

  double getD2LikelihoodForASite(size_t site);


  /**
   * @brief Get the likelihood for a site and for a state.
   *
   * @param site The site index to analyse.
   * @param state The state to consider.
   * @return The likelihood for site <i>site</i> and state <i>state</i>.
   */

  virtual double getLikelihoodForASiteForAState(size_t site, int state) = 0;

  virtual double getLogLikelihoodForASiteForAState(size_t site, int state) = 0;

  /**
   * @brief Get the logarithm of the likelihood for a site knowing its model class.
   *
   * @param site       The site index.
   * @param classIndex The model class index.
   * @return The likelihood for the specified site and model class.
   */
  virtual double getLikelihoodForASiteForAClass(size_t site, size_t classIndex) = 0;

  virtual double getLogLikelihoodForASiteForAClass(size_t site, size_t classIndex) = 0;

  /**
   * @brief Get the likelihood for a site knowing its model class and its ancestral state.
   *
   * @param site       The site index.
   * @param classIndex The model class index.
   * @param state      The ancestral state.
   * @return The likelihood for the specified site and model class and ancestral state..
   */

  virtual double getLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state) = 0;

  virtual double getLogLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state) = 0;
  
  /**
   * @brief Get the derivate of log-likelihood for the data set.
   *
   * Derivation is performed according to a given branch length,
   * which was specified by the last call to computeTreeDLogLikelihood.
   *
   * @return The derivate of log-likelihood for the data set.
   */
  
  virtual double getDLogLikelihood() = 0;
  
  /**
   * @brief Get the second order derivate of log-likelihood for the data set.
   *
   * Derivation is performed according to a given branch length,
   * which was specified by the last call to computeTreeD2LogLikelihood.
   *
   * @return The second order derivate of log-likelihood for the data set.
   */
  virtual double getD2LogLikelihood() = 0;

  /**
   * @brief Get the derivative of the loglikelihood for a site.
   *
   * @param site The site index to analyse.
   * @return The derivative of likelihood for site <i>site</i>.
   */

  virtual double getDLogLikelihoodForASite(size_t site) = 0;
  
  /**
   * @brief Get the second-order derivative of the loglikelihood for a site.
   *
   * @param site The site index to analyse.
   * @return The second-order derivative of likelihood for site <i>site</i>.
   */

  virtual double getD2LogLikelihoodForASite(size_t site) = 0;  

  /**
   * @brief update the likelihood dependencies (but does not compute).
   *  Is necessary befire computation.
   *
   */
  
  virtual void updateLikelihood() = 0;
  
  /**
   * @brief Perform a likelihood computation.
   */
  
  virtual void computeTreeLikelihood() = 0;

  
  virtual void computeLikelihoodsAtNode(int nodeId) = 0;

public:
  /**
   * @brief Initiate a derivative log-likelihood computation.
   *
   * @param variable The name of a parameter according to which
   * derivatives can been computed. If variable is not valid, the
   * derivatives set to 0.
   *
   */
  
  virtual void computeTreeDLogLikelihood(const std::string& variable) = 0;

  /**
   * @brief Initiate a second-order derivative log-likelihood
   * computation.
   *
   * @param variable The name of a parameter according to which
   * second-order derivatives can been computed. If variable is not
   * valid, the derivatives set to 0.
   */
  
  virtual void computeTreeD2LogLikelihood(const std::string& variable) = 0;

};

} // end of namespace bpp.

#endif //_LIKELIHOOD_TREE_CALCULATION_H_


