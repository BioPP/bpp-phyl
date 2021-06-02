//
// File: SingleProcessPhyloLikelihood.h
// Authors: François Gindraud, Laurent Guéguen (2017)
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef SINGLE_PROCESS_PHYLOLIKELIHOOD_H
#define SINGLE_PROCESS_PHYLOLIKELIHOOD_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>

#include "../DataFlow/DataFlowNumeric.h"
#include "../DataFlow/Parameter.h"
#include "../DataFlow/LikelihoodCalculationSingleProcess.h"

#include "SingleDataPhyloLikelihood.h"

#include <unordered_map>

/* This file contains wrappers.
 * They are used to bridge the gap between bpp::dataflow stuff and the rest of bpp.
 *
 */

namespace bpp {

  /* Wraps a dataflow graph as a function: resultNode = f(variableNodes).
   *
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
     *
     **/

    size_t nProc_;
      
    /**
     * @brief For Dataflow computing
     *
     */
      
    mutable std::unordered_map<std::string, ValueRef<Eigen::RowVectorXd>> firstOrderDerivativeVectors_;

    mutable std::unordered_map<std::pair<std::string, std::string>, ValueRef<Eigen::RowVectorXd>,
                               StringPairHash>
    secondOrderDerivativeVectors_;

  public:
    SingleProcessPhyloLikelihood (Context & context,
                                  std::shared_ptr<LikelihoodCalculationSingleProcess> likCal,
                                  const ParameterList & variableNodes,
                                  size_t nProc = 0, size_t nData=0) :
      AbstractPhyloLikelihood(context),
      AbstractAlignedPhyloLikelihood(context, likCal->getNumberOfSites()),
      AbstractSingleDataPhyloLikelihood(context, likCal->getNumberOfSites(), likCal->getStateMap().getNumberOfModelStates(), nData),
      AbstractParametrizable(""),
      likCal_(likCal), nProc_(nProc)
    {
      shareParameters_(variableNodes);
    }

    /*
     * @brief: the parameters are those of the LikelihoodCalculation
     */
      
    SingleProcessPhyloLikelihood (Context & context,
                                  std::shared_ptr<LikelihoodCalculationSingleProcess> likCal,
                                  size_t nProc = 0, size_t nData=0)
      :
      AbstractPhyloLikelihood(context),
      AbstractAlignedPhyloLikelihood(context, likCal->getNumberOfSites()),
      AbstractSingleDataPhyloLikelihood(context, likCal->getNumberOfSites(), likCal->getStateMap().getNumberOfModelStates(), nData),
      AbstractParametrizable(""),
      likCal_(likCal), nProc_(nProc)
    {
      shareParameters_(likCal_->getParameters());
    }

    // Legacy boilerplate
    SingleProcessPhyloLikelihood * clone () const override {
      throw Exception("SingleProcessPhyloLikelihood::clone should not be called.");
      return new SingleProcessPhyloLikelihood (*this); }

    void setData(const AlignedValuesContainer& sites, size_t nData = 0) override
    {
      AbstractSingleDataPhyloLikelihood::setData(sites, nData);  
      getLikelihoodCalculationSingleProcess()->setData(sites);
    }

    /**
     * @brief return a pointer to the compressed data. 
     *
     */
      
    const AlignedValuesContainer* getShrunkData() const {
      return getLikelihoodCalculationSingleProcess()->getShrunkData();
    }

    /**
     * @brief return a pointer to the original  data. 
     *
     */
      
    const AlignedValuesContainer* getData() const override {
      return getLikelihoodCalculationSingleProcess()->getData();
    }

    size_t getNumberOfSites() const override {
      return getLikelihoodCalculationSingleProcess()->getNumberOfSites();
    }

    size_t getNumberOfDistinctSites() const {
      return getLikelihoodCalculationSingleProcess()->getNumberOfDistinctSites();
    }

    /**
     * @brief Get the number of model classes.
     *
     */

    size_t getNumberOfClasses() const { return getSubstitutionProcess().getNumberOfClasses(); }

    
    const Alphabet* getAlphabet() const override {
      return getLikelihoodCalculationSingleProcess()->getStateMap().getAlphabet();
    }

    /*
     * @brief Get the ParametrizablePhyloTree.
     *
     * Warning: the branch lengths may not be up to date with those of
     * the LikelihoodCalculationSingleProcess.
     *
     */
    
    const ParametrizablePhyloTree& getTree() const {
      return getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getParametrizablePhyloTree(); }

    /*
     * @brief Return the ref to the SubstitutionProcess used to build
     * the phylolikelihood.
     *
     * Warning; the process parameter values may not be up to date
     * with some of the LikelihoodCalculationSingleProcess
     *
     */
    
    const SubstitutionProcess& getSubstitutionProcess() const {
      return getLikelihoodCalculationSingleProcess()->getSubstitutionProcess();
    }

    size_t getSubstitutionProcessNumber() const { return nProc_; }

    /**
     * @return 'true' is the likelihood function has been initialized.
     */

    bool isInitialized() const override {
      return getLikelihoodCalculationSingleProcess()->getData();
    };

    /**
     * @}
     */

    /**
     * @name Retrieve some particular independent parameters subsets.
     *
     * @{
     */
    

    ParameterList getNonDerivableParameters() const override
    {
      return getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getNonDerivableParameters();
    }

    ParameterList getDerivableParameters() const override
    {
      return getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getBranchLengthParameters(true);
    }

    /**
     * @brief Get the independent branch lengths parameters.
     *
     * @return A ParameterList with all branch lengths.
     */

    ParameterList getBranchLengthParameters() const override
    {
      return getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getBranchLengthParameters(true);
    }
      
    /**
     * @brief Get the independent parameters associated to substitution model(s).
     *
     * @return A ParameterList.
     */

    ParameterList getSubstitutionModelParameters() const override
    {
      return getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getSubstitutionModelParameters(true);
    }

    /**
     * @brief Get the independent parameters associated to the rate distribution(s).
     *
     * @return A ParameterList.
     */

    ParameterList getRateDistributionParameters() const override
    {
      return getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getRateDistributionParameters(true);
    }
      
    /**
     * @brief Get the independent parameters associated to the root
     * frequencies(s).
     *
     * @return A ParameterList.
     */
      
    ParameterList getRootFrequenciesParameters() const override
    {
      return getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getRootFrequenciesParameters(true);
    }
      
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
    
    VVdouble getLikelihoodPerSitePerClass() const;
    
    std::vector<size_t> getClassWithMaxPostProbPerSite() const;

    Vdouble getPosteriorRatePerSite() const;

    Vdouble getPosteriorStateFrequencies(uint nodeId);
  };
  
} // namespace bpp

#endif // SINGLE_PROCESS_PHYLOLIKELIHOOD_H
