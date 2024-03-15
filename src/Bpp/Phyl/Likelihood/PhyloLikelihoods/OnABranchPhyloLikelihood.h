// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_ONABRANCHPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_ONABRANCHPHYLOLIKELIHOOD_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>
#include <unordered_map>

#include "../DataFlow/DataFlowNumeric.h"
#include "../DataFlow/LikelihoodCalculationOnABranch.h"
#include "../DataFlow/Parameter.h"
#include "AlignedPhyloLikelihood.h"

/* This file contains wrappers.
 * They are used to bridge the gap between bpp::dataflow stuff and the rest of bpp.
 *
 */

namespace bpp
{
/**
 * @brief Wraps a dataflow graph as a function: resultNode = f(variableNodes).
 */
  class OnABranchPhyloLikelihood :
    public AbstractAlignedPhyloLikelihood,
    public AbstractParametrizable
  {
  protected:
    // Store nodes
    mutable std::shared_ptr<LikelihoodCalculationOnABranch> likCal_;

    /**
     * @brief For Dataflow computing
     */
    mutable std::unordered_map<std::string, ValueRef<RowLik> > firstOrderDerivativeVectors_;

    mutable std::unordered_map<std::pair<std::string, std::string>, ValueRef<RowLik>,
                               StringPairHash>
    secondOrderDerivativeVectors_;

  public:
    OnABranchPhyloLikelihood (Context& context,
                              std::shared_ptr<LikelihoodCalculationSingleProcess> likCal,
                              uint edgeId,
                              const ParameterList& variableNodes) :
      AbstractPhyloLikelihood(context),
      AbstractAlignedPhyloLikelihood(context, likCal->getNumberOfSites()),
      AbstractParametrizable(""),
      likCal_(std::make_shared<LikelihoodCalculationOnABranch>(context, *likCal, edgeId))
    {
      shareParameters_(variableNodes);
    }
              
    OnABranchPhyloLikelihood (Context& context,
                              std::shared_ptr<LikelihoodCalculationSingleProcess> likCal,
                              uint edgeId) :
      AbstractPhyloLikelihood(context),
      AbstractAlignedPhyloLikelihood(context, likCal->getNumberOfSites()),
      AbstractParametrizable(""),
      likCal_(std::make_shared<LikelihoodCalculationOnABranch>(context, *likCal, edgeId))
    {
    }


    OnABranchPhyloLikelihood (Context& context,
                              std::shared_ptr<LikelihoodCalculationOnABranch> likCal,
                              const ParameterList& variableNodes) :
      AbstractPhyloLikelihood(context),
      AbstractAlignedPhyloLikelihood(context, likCal->getNumberOfSites()),
      AbstractParametrizable(""),
      likCal_(likCal)
    {
      shareParameters_(variableNodes);
    }

    /**
     * @brief: the parameters the independent parameters of the LikelihoodCalculation
     */
    OnABranchPhyloLikelihood (Context& context,
                              std::shared_ptr<LikelihoodCalculationOnABranch> likCal) :
      AbstractPhyloLikelihood(context),
      AbstractAlignedPhyloLikelihood(context, likCal->getNumberOfSites()),
      AbstractParametrizable(""),
      likCal_(likCal)
    {
#ifdef DEBUG
      std::cerr << "OnABranchPhyloLikelihood(context, OnABranchPhyloLikelihood)" << std::endl;
#endif
      shareParameters_(likCal_->getIndependentParameters());
#ifdef DEBUG
      std::cerr << "OnABranchPhyloLikelihood(context, OnABranchPhyloLikelihood)" << std::endl;
#endif
    }

    // Legacy boilerplate
    OnABranchPhyloLikelihood* clone () const override
    {
      throw Exception("OnABranchPhyloLikelihood::clone should not be called.");
      return new OnABranchPhyloLikelihood (*this);
    }


    size_t getNumberOfSites() const override
    {
      return likelihoodCalculationOnABranch().getNumberOfSites();
    }

    size_t getNumberOfDistinctSites() const
    {
      return likelihoodCalculationOnABranch().getNumberOfDistinctSites();
    }

    /**
     * @brief Get the number of model classes.
     *
     */
    // size_t getNumberOfClasses() const { return substitutionProcess().getNumberOfClasses(); }


    // std::shared_ptr<const Alphabet> getAlphabet() const override
    // {
    //   return likelihoodCalculationOnABranch.getAlphabet();
    // }

    /**
     * @return 'true' is the likelihood function has been initialized.
     */
  
    bool isInitialized() const override
    {
      return likelihoodCalculationOnABranch().isInitialized();
    }

    void setModel(std::shared_ptr<ConfiguredModel> model)
    {
      likCal_->setModel(model);
    }
    
    /**
     * @name Retrieve some particular independent parameters subsets.
     *
     * @{
     */
    ParameterList getNonDerivableParameters() const override
    {
      return ParameterList();
    }

    ParameterList getDerivableParameters() const override
    {
      return ParameterList();
    }

    /**
     * @brief Get the independent branch lengths parameters.
     *
     * @return A ParameterList with all branch lengths.
     */
    ParameterList getBranchLengthParameters() const override
    {
      return ParameterList();
    }

    /**
     * @brief Get the independent parameters associated to substitution model(s).
     *
     * @return A ParameterList.
     */
    ParameterList getSubstitutionModelParameters() const override
    {
      return ParameterList();
    }

    /**
     * @brief Get the independent parameters associated to the rate distribution(s).
     *
     * @return A ParameterList.
     */
    ParameterList getRateDistributionParameters() const override
    {
      return ParameterList();
    }

    /**
     * @brief Get the independent parameters associated to the root
     * frequencies(s).
     *
     * @return A ParameterList.
     */
    ParameterList getRootFrequenciesParameters() const override
    {
      return ParameterList();
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

    LikelihoodCalculationOnABranch& likelihoodCalculationOnABranch() const
    {
      return *likCal_;
    }
  
    std::shared_ptr<LikelihoodCalculationOnABranch> getLikelihoodCalculationOnABranch() const
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
        auto vector = getLikelihoodCalculationOnABranch()->getSiteLikelihoods(true)->deriveAsValue (context_, accessVariableNode (variable));
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

    // VVdouble getPosteriorProbabilitiesPerSitePerClass() const;

    // Vdouble getPosteriorProbabilitiesForSitePerClass(size_t pos) const;

    /*
     *@brief return the likelihood of rate classes on each site.
     *
     *@return 2D-vector sites x classes
     */

    VVDataLik getLikelihoodPerSitePerClass() const;

    std::vector<size_t> getClassWithMaxPostProbPerSite() const;

    // Vdouble getPosteriorRatePerSite() const;

    // Vdouble getPosteriorStateFrequencies(uint nodeId);
  };
} // namespace bpp

#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_ONABRANCHPHYLOLIKELIHOOD_H
