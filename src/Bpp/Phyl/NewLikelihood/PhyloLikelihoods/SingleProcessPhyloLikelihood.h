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
  private:
    // Cache generated nodes representing derivatives, to avoid recreating them every time.
    // Using the mutable keyword because the table must be changed even in const methods.
    struct StringPairHash {
      std::size_t operator() (const std::pair<std::string, std::string> & p) const {
        std::hash<std::string> strHash{};
        return strHash (p.first) ^ (strHash (p.second) << 1);
      }
    };
      
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
                                  size_t nProc = 0, size_t nData=0)
      :
      AbstractPhyloLikelihood(context),
      AbstractAlignedPhyloLikelihood(context, likCal->getNumberOfSites()),
      AbstractSingleDataPhyloLikelihood(context, likCal->getNumberOfSites(), likCal->getSubstitutionProcess().getNumberOfStates(), nData),
      AbstractParametrizable(""),
      likCal_(likCal), nProc_(nProc)
    {
      shareParameters(variableNodes);
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
      AbstractSingleDataPhyloLikelihood(context, likCal->getNumberOfSites(), likCal->getSubstitutionProcess().getNumberOfStates(), nData),
      AbstractParametrizable(""),
      likCal_(likCal), nProc_(nProc)
    {
      shareParameters(likCal_->getParameters());
    }

    // Legacy boilerplate
    SingleProcessPhyloLikelihood * clone () const override {
      throw Exception("SingleProcessPhyloLikelihood::clone should not be called.");
      return new SingleProcessPhyloLikelihood (*this); }

    void setData(const AlignedValuesContainer& sites, size_t nData = 0)
    {
      AbstractSingleDataPhyloLikelihood::setData(sites, nData);  
      getLikelihoodCalculation()->setData(sites);
    }

    /**
     * @brief return a pointer to the compressed data. 
     *
     */
      
    const AlignedValuesContainer* getShrunkData() const {
      return getLikelihoodCalculation()->getShrunkData();
    }

    /**
     * @brief return a pointer to the original  data. 
     *
     */
      
    const AlignedValuesContainer* getData() const {
      return getLikelihoodCalculation()->getData();
    }

    size_t getNumberOfSites() const {
      return getLikelihoodCalculation()->getNumberOfSites();
    }

    size_t getNumberOfDistinctSites() const {
      return getLikelihoodCalculation()->getNumberOfDistinctSites();
    }
    
    const Alphabet* getAlphabet() const {
      return getLikelihoodCalculation()->getStateMap().getAlphabet();
    }

    const ParametrizablePhyloTree& getTree() const {
      return getLikelihoodCalculation()->getSubstitutionProcess().getParametrizablePhyloTree(); }

    const SubstitutionProcess& getSubstitutionProcess() const {
      return getLikelihoodCalculation()->getSubstitutionProcess();
    }

    size_t getSubstitutionProcessNumber() const { return nProc_; }

    /**
     * @return initialize the likelihood function.
     */
      
    void initialize() {};
      
    /**
     * @return 'true' is the likelihood function has been initialized.
     */

    bool isInitialized() const {
      return getLikelihoodCalculation()->getData();
    };

    /**
     * @}
     */

    /**
     * @name The likelihood functions.
     *
     * @{
     */
      
    /**
     * @brief Get the logarithm of the likelihood for the whole dataset.
     *
     * @return The logarithm of the likelihood of the dataset.
     */
    
    double getLogLikelihood() const {
      return -getLikelihoodCalculation()->getLikelihood()->getTargetValue ();
    }
    
    /**
     * @brief Compute the derivates of the LogLikelihood.
     *
     */

    void computeDLogLikelihood_(const std::string& variable) const {};

    
    void computeD2LogLikelihood_(const std::string& variable) const {};
    
    /**
     * @brief Get the derivates of the LogLikelihood.
     *
     */

    double getDLogLikelihood(const std::string& variable) const
    {
      return getFirstOrderDerivative(variable);
    }

    double getD2LogLikelihood(const std::string& variable) const
    {
      return getSecondOrderDerivative(variable);
    }

    /** @} */

    /**
     * @name Retrieve some particular independent parameters subsets.
     *
     * @{
     */
    
    /**
     * @brief Get the independent branch lengths parameters.
     *
     * @return A ParameterList with all branch lengths.
     */

    ParameterList getBranchLengthParameters() const
    {
      return getLikelihoodCalculation()->getSubstitutionProcess().getBranchLengthParameters(true);
    }
      
    /**
     * @brief Get the independent parameters associated to substitution model(s).
     *
     * @return A ParameterList.
     */

    ParameterList getSubstitutionModelParameters() const
    {
      return getLikelihoodCalculation()->getSubstitutionProcess().getSubstitutionModelParameters(true);
    }

    /**
     * @brief Get the independent parameters associated to the rate distribution(s).
     *
     * @return A ParameterList.
     */

    ParameterList getRateDistributionParameters() const
    {
      return getLikelihoodCalculation()->getSubstitutionProcess().getRateDistributionParameters(true);
    }
      
    /**
     * @brief Get the independent parameters associated to the root
     * frequencies(s).
     *
     * @return A ParameterList.
     */
      
    ParameterList getRootFrequenciesParameters() const
    {
      return getLikelihoodCalculation()->getSubstitutionProcess().getRootFrequenciesParameters(true);
    }
      
    /**
     * @brief All independent non derivable parameters.
     *
     * Usually, this contains all substitution model parameters and rate distribution.
     *
     * @return A ParameterList.
     */

    ParameterList getNonDerivableParameters() const
    {
      return getLikelihoodCalculation()->getSubstitutionProcess().getNonDerivableParameters();
    }
     

    /** @} */


    std::shared_ptr<LikelihoodCalculationSingleProcess> getLikelihoodCalculation() const
    {
      return likCal_;
    }

    ValueRef<double> getLikelihoodNode() const
    {
      return getLikelihoodCalculation()->getLikelihood();
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
        auto vector = getLikelihoodCalculation()->getSiteLikelihoods(true)->deriveAsValue (context_, accessVariableNode (variable));
        firstOrderDerivativeVectors_.emplace (variable, vector);
        return vector;
      }
    }
      
    ValueRef<Eigen::RowVectorXd> getSecondOrderDerivativeVector (const std::string & variable) const {
      return getSecondOrderDerivativeVector (variable, variable);
    }

    ValueRef<Eigen::RowVectorXd>  getSecondOrderDerivativeVector (const std::string & variable1,
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
     * @brief Get the likelihood for a site.
     *
     * @param site The site index to analyse.
     * @return The likelihood for site <i>site</i>.
     */

    double getLikelihoodForASite(size_t site) const
    {
      return getLikelihoodCalculation()->getLikelihoodForASite(site);
    }
    
    /**
     * @brief Get the log likelihood for a site, and its derivatives.
     *
     * @param site The site index to analyse.
     * @return The (D)log likelihood for site <i>site</i>.
     */

    double getLogLikelihoodForASite(size_t site) const
    {
      return std::log(getLikelihoodForASite(site));
    }
  
    double getDLogLikelihoodForASite(const std::string& variable, size_t site) const
    {
      throw Exception("SingleProcessPhyloLikelihood::getDLogLikelihoodForASite not finished : ask developpers.");
      return getFirstOrderDerivativeVector(variable)->getTargetValue()[site];
    }
    
    double getD2LogLikelihoodForASite(const std::string& variable, size_t site) const
    {
      throw Exception("SingleProcessPhyloLikelihood::getD2LogLikelihoodForASite not finished : ask developpers.");
      return getSecondOrderDerivativeVector(variable)->getTargetValue()[site];
    }
  
    /**
     * @brief Get the likelihood for each site.
     *
     *@return A vector with all likelihoods for each site.
     *
     */

    Vdouble getLikelihoodPerSite() const;

    /**
     * @brief Get the posterior probabilities of each class, for
     * each site.
     *
     * @return A 2D-vector (Sites X Class) of all posterior
     * probabilities.
     */

    VVdouble getPosteriorProbabilitiesPerClass() const;
      
    Vdouble getPosteriorProbabilitiesForSitePerClass(size_t pos) const;

  };
  
} // namespace bpp

#endif // SINGLE_PROCESS_PHYLOLIKELIHOOD_H
