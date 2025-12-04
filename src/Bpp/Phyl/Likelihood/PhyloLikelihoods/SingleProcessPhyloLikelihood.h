// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SINGLEPROCESSPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SINGLEPROCESSPHYLOLIKELIHOOD_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>
#include <unordered_map>

#include "../DataFlow/DataFlowNumeric.h"
#include "../DataFlow/LikelihoodCalculationSingleProcess.h"
#include "../DataFlow/Parameter.h"
#include "SingleDataPhyloLikelihood.h"

/* This file contains wrappers.
 * They are used to bridge the gap between bpp::dataflow stuff and the rest of bpp.
 *
 */

namespace bpp
{
/**
 * @brief Wraps a dataflow graph as a function: resultNode = f(variableNodes).
 */
class SingleProcessPhyloLikelihood :
  public AbstractSingleDataPhyloLikelihood,
  public AbstractParametrizable
{
protected:
  // Store nodes
  mutable std::shared_ptr<LikelihoodCalculationSingleProcess> likCal_;

  /**
   * @brief the Substitution Process number
   */
  size_t nProc_;

  /**
   * @brief For Dataflow computing
   */
  mutable std::unordered_map<std::string, ValueRef<RowLik>> firstOrderDerivativeVectors_;

  mutable std::unordered_map<std::pair<std::string, std::string>, ValueRef<RowLik>, StringPairHash> secondOrderDerivativeVectors_;

public:
  SingleProcessPhyloLikelihood (Context& context,
      std::shared_ptr<LikelihoodCalculationSingleProcess> likCal,
      const ParameterList& variableNodes,
      size_t nProc = 0, size_t nData = 0) :
    AbstractPhyloLikelihood(context),
    AbstractAlignedPhyloLikelihood(context, likCal->getNumberOfSites()),
    AbstractSingleDataPhyloLikelihood(context, likCal->getNumberOfSites(), likCal->stateMap().getNumberOfModelStates(), nData),
    AbstractParametrizable(""),
    likCal_(likCal),
    nProc_(nProc),
    firstOrderDerivativeVectors_(),
    secondOrderDerivativeVectors_()
  {
    shareParameters_(variableNodes);
  }

  /**
   * @brief: the parameters the independent parameters of the LikelihoodCalculation
   */
  SingleProcessPhyloLikelihood (Context& context,
      std::shared_ptr<LikelihoodCalculationSingleProcess> likCal,
      size_t nProc = 0, size_t nData = 0) :
    AbstractPhyloLikelihood(context),
    AbstractAlignedPhyloLikelihood(context, likCal->getNumberOfSites()),
    AbstractSingleDataPhyloLikelihood(context, likCal->getNumberOfSites(), likCal->stateMap().getNumberOfModelStates(), nData),
    AbstractParametrizable(""),
    likCal_(likCal),
    nProc_(nProc),
    firstOrderDerivativeVectors_(),
    secondOrderDerivativeVectors_()
  {
#ifdef DEBUG
    std::cerr << "SingleProcessPhyloLikelihood(context, LikelihoodCalculationSingleProcess)" << std::endl;
#endif
    shareParameters_(likCal_->getIndependentParameters());
#ifdef DEBUG
    std::cerr << "singleprocessphylolikelihood(context, likelihoodcalculationsingleprocess)" << std::endl;
#endif
  }

  // Legacy boilerplate
  SingleProcessPhyloLikelihood* clone () const override
  {
    throw Exception("SingleProcessPhyloLikelihood::clone should not be called.");
    return new SingleProcessPhyloLikelihood (*this);
  }

  void setData(std::shared_ptr<const AlignmentDataInterface> sites, size_t nData = 0) override
  {
    AbstractSingleDataPhyloLikelihood::setData(sites, nData);
    likelihoodCalculationSingleProcess().setData(sites);
  }

  /**
   * @brief return a pointer to the compressed data.
   *
   */
  std::shared_ptr<const AlignmentDataInterface> getShrunkData() const
  {
    return likelihoodCalculationSingleProcess().getShrunkData();
  }

  /**
   * @brief return a pointer to the original  data.
   *
   */
  std::shared_ptr<const AlignmentDataInterface> getData() const override
  {
    return likelihoodCalculationSingleProcess().getData();
  }

  size_t getNumberOfSites() const override
  {
    return likelihoodCalculationSingleProcess().getNumberOfSites();
  }

  size_t getNumberOfDistinctSites() const
  {
    return likelihoodCalculationSingleProcess().getNumberOfDistinctSites();
  }

  /**
   * @brief Get the number of model classes.
   *
   */
  size_t getNumberOfClasses() const { return substitutionProcess().getNumberOfClasses(); }


  std::shared_ptr<const Alphabet> getAlphabet() const override
  {
    return likelihoodCalculationSingleProcess().stateMap().getAlphabet();
  }

  /*
   * @brief Get the ParametrizablePhyloTree.
   *
   * Warning: the branch lengths may not be up to date with those of
   * the LikelihoodCalculationSingleProcess.
   *
   */
  std::shared_ptr<const ParametrizablePhyloTree> tree() const
  {
    return likelihoodCalculationSingleProcess().substitutionProcess().getParametrizablePhyloTree();
  }

  /**
   * @brief Return the ref to the SubstitutionProcess used to build
   * the phylolikelihood.
   *
   * Warning; the process parameter values may not be up to date
   * with some of the LikelihoodCalculationSingleProcess
   */
  const SubstitutionProcessInterface& substitutionProcess() const
  {
    return likelihoodCalculationSingleProcess().substitutionProcess();
  }

  /**
   * @brief Return a smarter pointer to the SubstitutionProcess used to build
   * the phylolikelihood.
   *
   * Warning; the process parameter values may not be up to date
   * with some of the LikelihoodCalculationSingleProcess
   */
  std::shared_ptr<const SubstitutionProcessInterface> getSubstitutionProcess() const
  {
    return likelihoodCalculationSingleProcess().getSubstitutionProcess();
  }

  size_t getSubstitutionProcessNumber() const { return nProc_; }

  /**
   * @return 'true' is the likelihood function has been initialized.
   */
  bool isInitialized() const override
  {
    return likelihoodCalculationSingleProcess().isInitialized();
  }

  /**
   * @name Retrieve some particular independent parameters subsets.
   *
   * @{
   */
  ParameterList getNonDerivableParameters() const override
  {
    return likelihoodCalculationSingleProcess().substitutionProcess().getNonDerivableParameters();
  }

  ParameterList getDerivableParameters() const override
  {
    return likelihoodCalculationSingleProcess().substitutionProcess().getBranchLengthParameters(true);
  }

  /**
   * @brief Get the independent branch lengths parameters.
   *
   * @return A ParameterList with all branch lengths.
   */
  ParameterList getBranchLengthParameters() const override
  {
    return likelihoodCalculationSingleProcess().substitutionProcess().getBranchLengthParameters(true);
  }

  /**
   * @brief Get the independent parameters associated to substitution model(s).
   *
   * @return A ParameterList.
   */
  ParameterList getSubstitutionModelParameters() const override
  {
    return likelihoodCalculationSingleProcess().substitutionProcess().getSubstitutionModelParameters(true);
  }

  /**
   * @brief Get the independent parameters associated to the rate distribution(s).
   *
   * @return A ParameterList.
   */
  ParameterList getRateDistributionParameters() const override
  {
    return likelihoodCalculationSingleProcess().substitutionProcess().getRateDistributionParameters(true);
  }

  /**
   * @brief Get the independent parameters associated to the root
   * frequencies(s).
   *
   * @return A ParameterList.
   */
  ParameterList getRootFrequenciesParameters() const override
  {
    return likelihoodCalculationSingleProcess().substitutionProcess().getRootFrequenciesParameters(true);
  }

  /**
   *
   * @}
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

  /**
   * @brief Get the posterior probabilities of each class, for
   * each site.
   *
   * @return A 2D-vector (Sites X Class) of all posterior
   * probabilities.
   */

  VVdouble getPosteriorProbabilitiesPerSitePerClass() const;

  Vdouble getPosteriorProbabilitiesForSitePerClass(size_t pos) const;

  /*
   *@brief return the likelihood of rate classes on each site.
   *
   *@return 2D-vector sites x classes
   */

  VVDataLik getLikelihoodPerSitePerClass() const;

  std::vector<size_t> getClassWithMaxPostProbPerSite() const;

  Vdouble getPosteriorRatePerSite() const;

  Vdouble getPosteriorStateFrequencies(unsigned int nodeId);
};
} // namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_SINGLEPROCESSPHYLOLIKELIHOOD_H
