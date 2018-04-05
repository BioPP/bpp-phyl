//
// File: OneProcessSequencePhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: mardi 28 avril 2015, à 12h 17
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

#ifndef _ONE_PROCESS_SEQUENCE_PHYLOLIKELIHOOD_H_
#define _ONE_PROCESS_SEQUENCE_PHYLOLIKELIHOOD_H_

#include "../../Tree/PhyloTree.h"
#include "../SitePartition.h"

#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/VectorTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignedValuesContainer.h>

#include "SequencePhyloLikelihood.h"
#include "../SubstitutionProcess.h"
#include "../OneProcessSequenceEvolution.h"
#include "../RecursiveLikelihoodTreeCalculation.h"

namespace bpp
{
/**
 * @brief The OneProcessSequencePhyloLikelihood class: phylogenetic
 * likelihood computation with a single process.
 *
 * This class implements likelihood calculation with a single
 * process/tree. It uses a unique LikelihoodTreeCalculation instance,
 * and implements the Function interface, dealing with parameters from
 * the associated SubstitutionProcess.
 */

  class OneProcessSequencePhyloLikelihood :
    public AbstractSequencePhyloLikelihood
  {
  private:
    /**
     * @brief to avoid the dynamic casts
     *
     */

    OneProcessSequenceEvolution& mSeqEvol_;

  protected:
    mutable std::unique_ptr<LikelihoodTreeCalculation> tlComp_;

  public:
    OneProcessSequencePhyloLikelihood(
      OneProcessSequenceEvolution& evol,
      size_t nSeqEvol = 0,
      bool verbose = true,
      bool patterns = true);

    OneProcessSequencePhyloLikelihood(
      const AlignedValuesContainer& data,
      OneProcessSequenceEvolution& evol,
      size_t nSeqEvol = 0,
      size_t nData = 0,
      bool verbose = true,
      bool patterns = true);

    OneProcessSequencePhyloLikelihood(const OneProcessSequencePhyloLikelihood& lik) :
      AbstractPhyloLikelihood(lik),
      AbstractAlignedPhyloLikelihood(lik),
      AbstractSequencePhyloLikelihood(lik),
      mSeqEvol_(lik.mSeqEvol_),
      tlComp_()
    {
      if (lik.tlComp_.get()) tlComp_.reset(lik.tlComp_->clone());
    }

    OneProcessSequencePhyloLikelihood& operator=(const OneProcessSequencePhyloLikelihood& lik)
    {
      AbstractSequencePhyloLikelihood::operator=(lik);
      mSeqEvol_ = lik.mSeqEvol_;

      if (lik.tlComp_.get())
        tlComp_.reset(lik.tlComp_->clone());
      else
        tlComp_.reset();

      return *this;
    }

    virtual ~OneProcessSequencePhyloLikelihood() {}

    OneProcessSequencePhyloLikelihood* clone() const { return new OneProcessSequencePhyloLikelihood(*this); }

  public:
    /**
     * @name Handling of data
     *
     * @{
     */
    void setData(const AlignedValuesContainer& sites, size_t nData = 0) 
    {
      AbstractSequencePhyloLikelihood::setData(sites, nData);
      tlComp_->setData(sites);
    }

    /**
     * @brief return a pointer to the compressed data.
     *
     */
    const AlignedValuesContainer* getData() const
    {
      return tlComp_->getData();
    }

    const Alphabet* getAlphabet() const
    {
      return tlComp_->getAlphabet();
    }

    /** @} */

    /**
     * @name Handling of substitution process
     *
     * @{
     */

    /**
     * @brief Get the substition process
     */

    const SubstitutionProcess& getSubstitutionProcess() const { return mSeqEvol_.getSubstitutionProcess(); }

    /**
     * @brief Get the number of model classes.
     *
     */

    size_t getNumberOfClasses() const { return mSeqEvol_.getSubstitutionProcess().getNumberOfClasses(); }

    /**
     * @brief Get the tree (topology and branch lengths).
     *
     * @return The tree of this OneProcessSequencePhyloLikelihood object.
     */

    const ParametrizablePhyloTree& getTree() const { return mSeqEvol_.getSubstitutionProcess().getParametrizablePhyloTree(); }


    /** @} */

  protected:
    void computeDLogLikelihood_(const std::string& variable) const;

    void computeD2LogLikelihood_(const std::string& variable) const;

  public:
    /**
     * @brief set it arrays should be computed in log.
     *
     */
    void setUseLog(bool useLog)
    {
      tlComp_->setAllUseLog(useLog);
    }

    /**
     * @return The underlying likelihood computation structure.
     */
    LikelihoodTreeCalculation* getLikelihoodCalculation() { return tlComp_.get(); }

    /**
     * @return The underlying likelihood data structure.
     */
    LikelihoodTree& getLikelihoodData() { return tlComp_->getLikelihoodData(); }

    /**
     * @return The underlying likelihood data structure.
     */
    const LikelihoodTree& getLikelihoodData() const { return tlComp_->getLikelihoodData(); }

  public:
    void updateLikelihood() const
    {
      if (computeLikelihoods_)
        tlComp_->updateLikelihood();
    }

    void computeLikelihood() const
    {
      if (computeLikelihoods_)
      {
        tlComp_->computeTreeLikelihood();
        computeLikelihoods_ = false;
      }
    }

    double getLogLikelihood() const
    {
      updateLikelihood();
      computeLikelihood();

      return tlComp_->getLogLikelihood();
    }

    double getDLogLikelihood(const std::string& variable) const
    {
      // check it is a "BrLen" variable

      if (!hasParameter(variable) || (variable.compare(0,5,"BrLen")!=0))
        return 0;

      return tlComp_->getDLogLikelihood();
    }

    double getD2LogLikelihood(const std::string& variable) const
    {
      // check it is a "BrLen" variable

      if (!hasParameter(variable) || (variable.compare(0,5,"BrLen")!=0))
        return 0;

      return tlComp_->getD2LogLikelihood();
    }

    double getLikelihoodForASite(size_t site) const
    {
      updateLikelihood();
      computeLikelihood();

      return tlComp_->getLikelihoodForASite(site);
    }

    double getLogLikelihoodForASite(size_t site) const
    {
      updateLikelihood();
      computeLikelihood();

      return tlComp_->getLogLikelihoodForASite(site);
    }

    double getDLogLikelihoodForASite(const std::string& variable, size_t site) const
    {
      // check it is a "BrLen" variable

      if (!hasParameter(variable) || (variable.compare(0,5,"BrLen")!=0))
        return 0;

      return tlComp_->getDLogLikelihoodForASite(site);
    }

    double getD2LogLikelihoodForASite(const std::string& variable, size_t site) const
    {
      // check it is a "BrLen" variable

      if (!hasParameter(variable) || (variable.compare(0,5,"BrLen")!=0))
        return 0;

      return tlComp_->getD2LogLikelihoodForASite(site);
    }

    /**
     * @brief Get the likelihood for each site and for each state.
     *
     * @return A 2d vector with all likelihoods for each site and for each state.
     */
    VVdouble getLikelihoodPerSitePerState() const;

    /**
     * @brief Get the likelihood for each site and each model class.
     *
     * @return A two-dimension vector with all likelihoods:
     * <code>V[i][j] =</code> likelihood of site i and model class j.
     */
    VVdouble getLikelihoodPerSitePerClass() const;

    /**
     * @brief Get the likelihood for each site and each model class and each state.
     *
     * @return A three-dimension vector with all likelihoods:
     * <code>V[i][j][k} =</code> likelihood of site i and model class j and state k.
     */

    VVVdouble getLikelihoodPerSitePerClassPerState() const;

    /** @} */

    /**
     * @brief Get the index (used for inner computations) of a given site (original alignment column).
     *
     * @param site An alignment position.
     * @return The site index corresponding to the given input alignment position.
     */
    // size_t getSiteIndex(size_t site) const throw (IndexOutOfBoundsException)
    // {
    //   return tlComp_->getSiteIndex(site);
    // }

    /**
     * Utilities
     *
     */

    VVdouble getPosteriorProbabilitiesPerClass() const;

    /**
     * @brief Get the posterior model class (the one with maximum posterior
     * probability) for each site.
     *
     * @return A vector with all model classes indexes.
     */

    std::vector<size_t> getClassWithMaxPostProbPerSite() const;

    Vdouble getPosteriorRatePerSite() const;

    /* @} */

  };
} // end of namespace bpp.

#endif  // _ONE_PROCESS_SEQUENCE_PHYLOLIKELIHOOD_H_
