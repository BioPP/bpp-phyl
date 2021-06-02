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

#include "../DataFlow/LikelihoodCalculationSingleProcess.h"

#include "SequencePhyloLikelihood.h"
#include "../SubstitutionProcess.h"
#include "../OneProcessSequenceEvolution.h"

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

    /**
     * @brief For Dataflow computing
     *
     */
      
    mutable std::unordered_map<std::string, ValueRef<Eigen::RowVectorXd>> firstOrderDerivativeVectors_;

    mutable std::unordered_map<std::pair<std::string, std::string>, ValueRef<Eigen::RowVectorXd>,
                               StringPairHash>
    secondOrderDerivativeVectors_;

  protected:
    mutable std::shared_ptr<LikelihoodCalculationSingleProcess> likCal_;

  public:
    OneProcessSequencePhyloLikelihood(
      Context& context,
      OneProcessSequenceEvolution& evol,
      size_t nSeqEvol = 0,
      bool verbose = true,
      bool patterns = true);

    OneProcessSequencePhyloLikelihood(
      Context& context,
      const AlignedValuesContainer& data,
      OneProcessSequenceEvolution& evol,
      size_t nSeqEvol = 0,
      size_t nData = 0,
      bool verbose = true,
      bool patterns = true);

    OneProcessSequencePhyloLikelihood(
      const AlignedValuesContainer& data,
      OneProcessSequenceEvolution& evol,
      CollectionNodes& collNodes,
      size_t nSeqEvol = 0,
      size_t nData = 0,
      bool verbose = true,
      bool patterns = true);

    OneProcessSequencePhyloLikelihood(const OneProcessSequencePhyloLikelihood& lik) :
      AbstractPhyloLikelihood(lik),
      AbstractAlignedPhyloLikelihood(lik),
      AbstractSequencePhyloLikelihood(lik),
      mSeqEvol_(lik.mSeqEvol_),
      likCal_(lik.likCal_)
    {
    }

    virtual ~OneProcessSequencePhyloLikelihood() {}

    OneProcessSequencePhyloLikelihood* clone() const override { return new OneProcessSequencePhyloLikelihood(*this); }

  public:
    /**
     * @name Handling of data
     *
     * @{
     */

    void setData(const AlignedValuesContainer& sites, size_t nData = 0) override
    {
      AbstractSequencePhyloLikelihood::setData(sites, nData);
      getLikelihoodCalculationSingleProcess()->setData(sites);
    }

    bool isInitialized() const override {
      return getLikelihoodCalculationSingleProcess()->getData();
    };

    /**
     * @brief return a pointer to the compressed data.
     *
     */
    const AlignedValuesContainer* getData() const override
    {
      return getLikelihoodCalculationSingleProcess()->getData();
    }

    const Alphabet* getAlphabet() const override
    {
      return getLikelihoodCalculationSingleProcess()->getStateMap().getAlphabet();
    }

    /** @} */

    /**
     * @name Handling of substitution process
     *
     * @{
     */

    /*
     * @brief Get the ParametrizablePhyloTree.
     *
     * Warning: the branch lengths may not be up to date with those of
     * the LikelihoodCalculationSingleProcess.
     *
     */

    const SubstitutionProcess& getSubstitutionProcess() const { return mSeqEvol_.getSubstitutionProcess(); }

    /**
     * @brief Get the number of model classes.
     *
     */

    size_t getNumberOfClasses() const { return mSeqEvol_.getSubstitutionProcess().getNumberOfClasses(); }

    /*
     * @brief Return the ref to the SubstitutionProcess used to build
     * the phylolikelihood.
     *
     * Warning; the process parameter values may not be up to date
     * with some of the LikelihoodCalculationSingleProcess
     *
     */

    const ParametrizablePhyloTree& getTree() const { return mSeqEvol_.getSubstitutionProcess().getParametrizablePhyloTree(); }


    /** @} */

  public:
    /**
     * @return The underlying likelihood computation structure.
     */
    
    std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation() const override
    {
      return likCal_;
    } 
    
    std::shared_ptr<AlignedLikelihoodCalculation> getAlignedLikelihoodCalculation() const override
    {
      return likCal_;
    } 

    std::shared_ptr<LikelihoodCalculationSingleProcess> getLikelihoodCalculationSingleProcess() const
    {
      return likCal_;
    }

    // Get nodes of derivatives directly
      
    ValueRef<Eigen::RowVectorXd> getFirstOrderDerivativeVector (const std::string & variable) const  {
      return firstOrderDerivativeVector(variable);
    }

    ValueRef<Eigen::RowVectorXd> firstOrderDerivativeVector (const std::string & variable) const {
      const auto it = firstOrderDerivativeVectors_.find (variable);
      if (it != firstOrderDerivativeVectors_.end ()) {
        return it->second;
      } else {
        auto vector = getLikelihoodCalculationSingleProcess()->getSiteLikelihoods(true)->deriveAsValue (context_, accessVariableNode (variable));
        firstOrderDerivativeVectors_.emplace (variable, vector);
        return vector;
      }
    }
      
    ValueRef<Eigen::RowVectorXd> getSecondOrderDerivativeVector (const std::string & variable) const {
      return getSecondOrderDerivativeVector (variable, variable);
    }

    ValueRef<Eigen::RowVectorXd> getSecondOrderDerivativeVector (const std::string & variable1,
                                                                 const std::string & variable2) const {
      return secondOrderDerivativeVector (variable1, variable2);
    }

    ValueRef<Eigen::RowVectorXd> secondOrderDerivativeVector (const std::string & variable1,
                                                              const std::string & variable2) const {
      const auto key = std::make_pair (variable1, variable2);
      const auto it = secondOrderDerivativeVectors_.find (key);
      if (it != secondOrderDerivativeVectors_.end ()) {
        return it->second;
      } else {
        // Reuse firstOrderDerivative() to generate the first derivative with caching
        auto vector =
          firstOrderDerivativeVector (variable1)->deriveAsValue (context_, accessVariableNode (variable2));
        secondOrderDerivativeVectors_.emplace (key, vector);
        return vector;
      }
    }

  public:
    /**
     * Utilities
     *
     */

    /*
     *@brief return the posterior probabilities of rate classes on each site.
     *
     *@return 2D-vector sites x classes
     */

    VVdouble getPosteriorProbabilitiesPerSitePerClass() const;

    Vdouble getPosteriorProbabilitiesForSitePerClass(size_t pos) const;

    /*
     *@brief return the likelihood of rate classes on each site.
     *
     *@return 2D-vector sites x classes
     */
    
    VVdouble getLikelihoodPerSitePerClass() const;
    
    std::vector<size_t> getClassWithMaxPostProbPerSite() const;

    Vdouble getPosteriorRatePerSite() const;

    Vdouble getPosteriorStateFrequencies(uint nodeId);

    /* @} */

  };
} // end of namespace bpp.

#endif  // _ONE_PROCESS_SEQUENCE_PHYLOLIKELIHOOD_H_
