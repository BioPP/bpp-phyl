//
// File: AbstractLikelihoodTreeCalculation.h
// Created by: Julien Dutheil, Laurent Guéguen
// Created on: mardi 23 juin 2015, à 14h 04
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

#ifndef _ABSTRACT_LIKELIHOOD_TREE_CALCULATION_H_
#define _ABSTRACT_LIKELIHOOD_TREE_CALCULATION_H_

#include "LikelihoodTreeCalculation.h"
#include "SubstitutionProcess.h"
#include "AbstractLikelihoodTree.h"

namespace bpp
{
/**
 * @brief Partial implementation of the LikelihoodTreeCalculation
 * interface.
 *
 */
  class AbstractLikelihoodTreeCalculation:
    public virtual LikelihoodTreeCalculation
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

    // if this is up to date as this knows it, but beware the process could
    // change without any sign!!
    
    bool up2date_;
    
    // booleans to say if the Dlikelihoods are null
  
    bool nullDLikelihood_;
    bool nullD2Likelihood_;
  
  private:

    /*
     * @brief for computation purpose
     *
     */
    mutable std::vector<double> vSites_;
      
  public:
    AbstractLikelihoodTreeCalculation(const SubstitutionProcess* process, bool verbose = true):
      process_(process),
      data_(0),
      nbSites_(0),
      nbDistinctSites_(0),
      nbStates_(process->getNumberOfStates()),
      nbClasses_(process->getNumberOfClasses()),
      initialized_(false),
      verbose_(verbose),
      up2date_(false),
      nullDLikelihood_(true),
      nullD2Likelihood_(true),
      vSites_()
    {
    }
  
    AbstractLikelihoodTreeCalculation(const AbstractLikelihoodTreeCalculation& tlc):
    process_(tlc.process_),
    data_(0),
    nbSites_(tlc.nbSites_),
    nbDistinctSites_(tlc.nbDistinctSites_),
    nbStates_(tlc.nbStates_),
    nbClasses_(tlc.nbClasses_),
    initialized_(tlc.initialized_),
    verbose_(tlc.verbose_),
    up2date_(tlc.up2date_),
    nullDLikelihood_(tlc.nullDLikelihood_),
    nullD2Likelihood_(tlc.nullD2Likelihood_),
    vSites_(tlc.vSites_)
    {
      if (tlc.data_.get()) data_.reset(tlc.data_->clone());
    }
  
    AbstractLikelihoodTreeCalculation& operator=(const AbstractLikelihoodTreeCalculation& tlc)
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
      up2date_                       = tlc.up2date_;
      nullDLikelihood_               = tlc.nullDLikelihood_;
      nullD2Likelihood_              = tlc.nullD2Likelihood_;
      vSites_                        = tlc.vSites_;
        
      return *this;
    }

    virtual ~AbstractLikelihoodTreeCalculation() {}

  public:

    bool isInitialized() const { return initialized_; }

    const Alphabet* getAlphabet() const throw (LikelihoodTreeCalculationNotInitializedException)
    {
      if (!initialized_)
        throw new LikelihoodTreeCalculationNotInitializedException("DoubleRecursiveLikelihoodTreeCalculation::getAlphabet().");
      return data_->getAlphabet();
    }

    size_t getSiteIndex(size_t site) const throw (LikelihoodTreeCalculationNotInitializedException, IndexOutOfBoundsException) {
      if (!initialized_)
        throw new LikelihoodTreeCalculationNotInitializedException("SingleRecursiveLikelihoodTreeCalculation::getSiteIndex().");
      return getLikelihoodData().getRootArrayPosition(site);
    }

    void setData(const SiteContainer& sites);

    const SiteContainer* getData() const
    {
      if (!initialized_)
        throw new LikelihoodTreeCalculationNotInitializedException("SingleRecursiveLikelihoodTreeCalculation::getData().");
      return data_.get();
    }
  
    const SubstitutionProcess* getSubstitutionProcess() const { return process_;}


    size_t getNumberOfDistinctSites() const {
      return nbDistinctSites_;
    }

    size_t getNumberOfSites() const {
      return nbSites_;
    }
      
    size_t getNumberOfStates() const {
      return nbStates_;
    }

    size_t getNumberOfClasses() const {
      return getLikelihoodData().getNumberOfClasses();
    }

    int getRootId() const
    {
      return process_->getTree().getRootNode()->getId();
    }
      
    /**
     * @return if the log-likelihood is stored at root in class classIndex.
     *
     */
       
    bool usesLogAtRoot(size_t classIndex) const
    {
      return getLikelihoodData().usesLogAtRoot(classIndex);
    }
      
    /**
     * @brief sets using log in all likelihood arrays.
     *
     */
    
    void setAllUseLog(bool useLog)
    {
      getLikelihoodData().setAllUseLog(useLog);
    }

    /*
     * @brief Retrieve the likelihood data.
     *
     */
       
      
    AbstractLikelihoodTree& getLikelihoodData() = 0;

    const AbstractLikelihoodTree& getLikelihoodData() const = 0;

      
    /*
     * @brief get DX(log-)likelihoods
     *
     * !!!! These methods do not check that computations are up to
     * date.
     *
     *
     */

    double getLogLikelihood();

    double getDLogLikelihood();

    double getD2LogLikelihood();
      
    /*
     * @brief get DX(log-)likelihoods at sites.
     *
     * !!!! These methods do not check that computations are up to
     * date.
     *
     *
     */
      
    double getLikelihoodForASite(size_t site);

    double getLikelihoodForASiteForAState(size_t site, int state);

    double getLikelihoodForASiteForAClass(size_t site, size_t classIndex);

    double getLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state);

      

    double getLogLikelihoodForASite(size_t site);

    double getLogLikelihoodForASiteForAState(size_t site, int state);

    double getLogLikelihoodForASiteForAClass(size_t site, size_t classIndex);

    double getLogLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state);


    double getDLikelihoodForASite(size_t site);

    double getD2LikelihoodForASite(size_t site);

      
    double getDLogLikelihoodForASite(size_t site)
    {
      // d(f(g(x)))/dx = dg(x)/dx . df(g(x))/dg :
      return getDLikelihoodForASite(site) / getLikelihoodForASite(site);
    }

    double getD2LogLikelihoodForASite(size_t site)
    {
      return getD2LikelihoodForASite(site) / getLikelihoodForASite(site)
        - pow( getDLikelihoodForASite(site) / getLikelihoodForASite(site), 2);
    }

    /**
     * @brief Print the likelihood array to terminal (debugging tool).
     *
     * @param likelihoodArray the likelihood array.
     */

    static void displayLikelihoodArray(const VVVdouble& likelihoodArray);

    /**
     * @brief compute ancestral frequencies
     *
     */

  public:
      
    /**
     * @brief Compute the expected ancestral frequencies of all
     * states at all (inner) nodes according to a Markov process
     * defined by a given substitution model.
     *
     * The computation is averaged over all sites. If the likelihood
     * object has no site partition, then the method will return the
     * same result as all single site numbers.
     *
     * @param frequencies [out] A map where to store the results, as
     *        a vector of double (the size of which being equal to
     *        the number of states in the model), and with nodes id
     *        as keys.

     * @param alsoForLeaves [opt] Tell if frequencies should also be
     *        estimated for terminal nodes.
     */

    void getAncestralFrequencies(
      std::map<int, std::vector<double> >& frequencies,
      bool alsoForLeaves = false);

    /**
     * @brief Compute the expected ancestral frequencies of all
     *        states at all (inner) nodes according to a given
     *        Markov process.
     *
     * The computation is performed for a given site. If the
     *       likelihood object has no site partition, then the
     *       method will return the same result for all positions.
     *
     * @param site the given site.
     *
     * @param frequencies [out] A map where to store the results, as
     *       a vector of double (the size of which being equal to
     *       the number of states in the model), and with nodes id
     *       as keys.
     *
     * @param alsoForLeaves [opt] Tell if frequencies should also be
     *       estimated for terminal nodes.
     *
     */

    void getAncestralFrequencies(
      size_t site,
      std::map<int, std::vector<double> >& frequencies,
      bool alsoForLeaves = false);

  private:
    /**
     * @brief Recursive method, for internal use only.
     *
     * @see getAncestralFrequencies()
     */

    void getAncestralFrequencies_(
      size_t siteIndex,
      size_t classIndex,
      int parentId,
      const std::vector<double>& ancestralFrequencies,
      std::map<int, std::vector<double> >& frequencies,
      bool alsoForLeaves);
 
  };

} // end of namespace bpp.

#endif  // _ABSTRACT_LIKELIHOOD_TREE_CALCULATION_H_

