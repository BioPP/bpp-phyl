// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_ONEPROCESSSEQUENCEPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_ONEPROCESSSEQUENCEPHYLOLIKELIHOOD_H

#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/VectorTools.h>

#include "../../Tree/PhyloTree.h"
#include "../SitePartition.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignmentData.h>

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
  public AbstractParametrizableSequencePhyloLikelihood
{
private:
  /**
   * @brief to avoid the dynamic casts
   */
  std::shared_ptr<OneProcessSequenceEvolution> mSeqEvol_;

  /**
   * @brief For Dataflow computing
   */
  mutable std::unordered_map<std::string, ValueRef<RowLik>> firstOrderDerivativeVectors_;

  mutable std::unordered_map<std::pair<std::string, std::string>, ValueRef<RowLik>,
      StringPairHash>
  secondOrderDerivativeVectors_;

protected:
  mutable std::shared_ptr<LikelihoodCalculationSingleProcess> likCal_;

public:
  OneProcessSequencePhyloLikelihood(
      Context& context,
      std::shared_ptr<OneProcessSequenceEvolution> evol,
      size_t nSeqEvol = 0);

  OneProcessSequencePhyloLikelihood(
      Context& context,
      std::shared_ptr<const AlignmentDataInterface> data,
      std::shared_ptr<OneProcessSequenceEvolution> evol,
      size_t nSeqEvol = 0,
      size_t nData = 0);

  OneProcessSequencePhyloLikelihood(
      std::shared_ptr<const AlignmentDataInterface> data,
      std::shared_ptr<OneProcessSequenceEvolution> evol,
      std::shared_ptr<CollectionNodes> collNodes,
      size_t nSeqEvol = 0,
      size_t nData = 0);

protected:
  OneProcessSequencePhyloLikelihood(const OneProcessSequencePhyloLikelihood& lik) = default;

  OneProcessSequencePhyloLikelihood& operator=(const OneProcessSequencePhyloLikelihood& lik) = default;

  OneProcessSequencePhyloLikelihood* clone() const override
  {
    return new OneProcessSequencePhyloLikelihood(*this);
  }

public:
  virtual ~OneProcessSequencePhyloLikelihood() {}

public:
  /**
   * @name Handling of data
   *
   * @{
   */
  void setData(std::shared_ptr<const AlignmentDataInterface> sites, size_t nData = 0) override
  {
    AbstractSequencePhyloLikelihood::setData(sites, nData);
    likelihoodCalculationSingleProcess().setData(sites);
  }

  bool isInitialized() const override
  {
    return likelihoodCalculationSingleProcess().isInitialized();
  }

  /**
   * @brief return a pointer to the compressed data.
   */
  std::shared_ptr<const AlignmentDataInterface> getData() const override
  {
    return getLikelihoodCalculationSingleProcess()->getData();
  }

  std::shared_ptr<const Alphabet> getAlphabet() const override
  {
    return getLikelihoodCalculationSingleProcess()->stateMap().getAlphabet();
  }

  /** @} */

  /**
   * @name Handling of substitution process
   *
   * @{
   */

  /**
   * @brief Get the ParametrizablePhyloTree.
   *
   * Warning: the branch lengths may not be up to date with those of
   * the LikelihoodCalculationSingleProcess.
   */
  const SubstitutionProcessInterface& substitutionProcess() const { return mSeqEvol_->substitutionProcess(); }

  std::shared_ptr<const SubstitutionProcessInterface> getSubstitutionProcess() const { return mSeqEvol_->getSubstitutionProcess(); }

  /**
   * @brief Get the number of model classes.
   */
  size_t getNumberOfClasses() const { return mSeqEvol_->substitutionProcess().getNumberOfClasses(); }

  /**
   * @brief Return the ref to the SubstitutionProcess used to build
   * the phylolikelihood.
   *
   * Warning; the process parameter values may not be up to date
   * with some of the LikelihoodCalculationSingleProcess
   */
  std::shared_ptr<const ParametrizablePhyloTree> tree() const
  {
    return mSeqEvol_->substitutionProcess().getParametrizablePhyloTree();
  }

  /** @} */

public:
  /**
   * @return The underlying likelihood computation structure.
   */
  LikelihoodCalculation& likelihoodCalculation() const override
  {
    return *likCal_;
  }

  std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation() const override
  {
    return likCal_;
  }

  AlignedLikelihoodCalculation& alignedLikelihoodCalculation() const override
  {
    return *likCal_;
  }

  std::shared_ptr<AlignedLikelihoodCalculation> getAlignedLikelihoodCalculation() const override
  {
    return likCal_;
  }

  LikelihoodCalculationSingleProcess& likelihoodCalculationSingleProcess() const
  {
    return *likCal_;
  }

  std::shared_ptr<LikelihoodCalculationSingleProcess> getLikelihoodCalculationSingleProcess() const
  {
    return likCal_;
  }

  // Get nodes of derivatives directly

  ValueRef<RowLik> getFirstOrderDerivativeVector (const std::string& variable) const
  {
    return firstOrderDerivativeVector(variable);
  }

  ValueRef<RowLik> firstOrderDerivativeVector (const std::string& variable) const
  {
    const auto it = firstOrderDerivativeVectors_.find (variable);
    if (it != firstOrderDerivativeVectors_.end ())
    {
      return it->second;
    }
    else
    {
      auto vector = getLikelihoodCalculationSingleProcess()->getSiteLikelihoods(true)->deriveAsValue (context_, accessVariableNode (variable));
      firstOrderDerivativeVectors_.emplace (variable, vector);
      return vector;
    }
  }

  ValueRef<RowLik> getSecondOrderDerivativeVector (const std::string& variable) const
  {
    return getSecondOrderDerivativeVector (variable, variable);
  }

  ValueRef<RowLik> getSecondOrderDerivativeVector (const std::string& variable1,
      const std::string& variable2) const
  {
    return secondOrderDerivativeVector (variable1, variable2);
  }

  ValueRef<RowLik> secondOrderDerivativeVector (const std::string& variable1,
      const std::string& variable2) const
  {
    const auto key = std::make_pair (variable1, variable2);
    const auto it = secondOrderDerivativeVectors_.find (key);
    if (it != secondOrderDerivativeVectors_.end ())
    {
      return it->second;
    }
    else
    {
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
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_ONEPROCESSSEQUENCEPHYLOLIKELIHOOD_H
