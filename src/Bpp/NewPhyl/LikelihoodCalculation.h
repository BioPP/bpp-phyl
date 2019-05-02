//
// File: LikelihoodCalculation.h
// Authors: François Gindraud, Laurent Guéguen (2018)
// Created: jeudi 28 février 2019, à 07h 22
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

#ifndef LIKELIHOOD_CALCULATION_H
#define LIKELIHOOD_CALCULATION_H

#include "Bpp/NewPhyl/Model.h"
#include "Bpp/NewPhyl/DiscreteDistribution.h"
#include "Bpp/NewPhyl/FrequenciesSet.h"
#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowCWise.h>

#include <Bpp/Seq/Container/AlignedValuesContainer.h>
#include "Bpp/Phyl/NewLikelihood/SubstitutionProcess.h"

/* This file contains temporary helpers and wrappers.
 * They are used to bridge the gap between bpp::dataflow stuff and the rest of bpp.
 * They have only been used (and thus tested) for a single likelihood example.
 * They do not deal with all of bpp features, which is why they are only temporary.
 *
 * Ultimately, stuff in this file should be changed to a new system to describe phylogenic computations, which
 * would generate dataflow graphs to do the actual computations.
 */

namespace bpp {

  /** Conditional likelihoods are stored in a matrix of sizes (nbState, nbSite).
   * Rows represents states (nucleotides, proteins or codon).
   * Columns represents sites (one site for each column).
   * Conditional likelihood is thus accessed by m(state,site) for an eigen matrix.
   * Eigen defaults to ColMajor matrices, meaning data is stored column by column.
   * So with the current layout, values for a site are grouped together for locality.
   *
   * A Transition matrix is a (nbState,nbState) matrix.
   * tm(fromState, toState) = probability of going to toState from fromState.
   * This matches the convention from TransitionModel::getPij_t().
   *
   * Equilibrium frequencies are stored as a RowVector(nbState) : matrix with 1 row and n columns.
   * This choice allows to reuse the MatrixProduct numeric node directly.
   *
   * Initial conditional likelihood for leaves (sequences on the tree) should be computed outside of the
   * dataflow graph, and provided as NumericConstant<MatrixXd>.
   */

  inline MatrixDimension conditionalLikelihoodDimension (std::size_t nbState, std::size_t nbSite) {
    return {Eigen::Index (nbState), Eigen::Index (nbSite)};
  }
  inline MatrixDimension transitionMatrixDimension (std::size_t nbState) {
    return {Eigen::Index (nbState), Eigen::Index (nbState)};
  }
  inline MatrixDimension equilibriumFrequenciesDimension (std::size_t nbState) {
    return rowVectorDimension (Eigen::Index (nbState));
  }

  namespace dataflow {

    class PhyloTree_BrRef;
    class ForwardLikelihoodTree;
    class BackwardLikelihoodTree;
    
    /** @brief likelihood = f(equilibriumFrequencies, rootConditionalLikelihood).
     * - likelihood: RowVector(site).
     * - equilibriumFrequencies: RowVector(state).
     * - rootConditionalLikelihood: Matrix(state, site).
     *
     * likelihood(site) = sum_state equFreqs(state) * rootConditionalLikelihood(state, site).
     * Using matrix multiply: likelihood = equilibriumFrequencies * rootConditionalLikelihood.
     */

    using LikelihoodFromRootConditional =
      MatrixProduct<Eigen::RowVectorXd, Eigen::RowVectorXd, Eigen::MatrixXd>;
    
    /** @brief totalLikelihood = product_site likelihood(site).
     * - likelihood: RowVector (site).
     * - totalLikelihood: Extended float.
     */
    
    using TotalLogLikelihood = SumOfLogarithms<Eigen::RowVectorXd>;

    /** @brief Conditionallikelihood = AboveConditionalLikelihood * BelowConditionalLikelihood
     *
     * lik(state, site) = abova(state, site) * below(state,site)
     * Using member wise multiply
     */

    using BuildConditionalLikelihood =
      CWiseMul<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

    using ConditionalLikelihood = Value<Eigen::MatrixXd>;
    
    using ConditionalLikelihoodTree = AssociationTreeGlobalGraphObserver<ConditionalLikelihood,uint>;
    
    class LikelihoodCalculation :
      public AbstractParametrizable
    {
    private:
      
      class RateCategoryTrees {
      public:
        std::shared_ptr<PhyloTree_BrRef> phyloTree;
        std::shared_ptr<ForwardLikelihoodTree> flt;
        std::shared_ptr<BackwardLikelihoodTree> blt;
        std::shared_ptr<ConditionalLikelihoodTree> lt;
      };
      
      Context& context_;

      /************************************/
      /* Dependencies */
      
      const SubstitutionProcess& process_;
      const AlignedValuesContainer* psites_;

      /************************************/
      /* DataFlow objects */
      
         
      std::shared_ptr<PhyloTree_BrRef> treeNode_;

      /* link towards a model (any) to get stateMap & nb of states */
      
      std::shared_ptr<ConfiguredModel> modelNode_;
      
      std::shared_ptr<ConfiguredFrequenciesSet> rootFreqsNode_;
      
      std::shared_ptr<ConfiguredDistribution> ratesNode_;

      ValueRef<Eigen::RowVectorXd> rFreqs_;

      /******************************************/
      /** Likelihoods  **/
      
          
      ValueRef<double> likelihood_;

      ValueRef<Eigen::RowVectorXd> siteLikelihoods_;

      /* Likelihood Trees with for all rate categories */
      std::vector<RateCategoryTrees> vRateCatTrees_;

    public:
      LikelihoodCalculation(Context & context,
                            const AlignedValuesContainer & sites,
                            const SubstitutionProcess& process);

      LikelihoodCalculation(Context & context,
                            const SubstitutionProcess& process);

      LikelihoodCalculation(const LikelihoodCalculation& lik);

      LikelihoodCalculation* clone() const
      {
        return new LikelihoodCalculation(*this);
      }
    
      ValueRef<double> getLikelihood() 
      {
        if (psites_ && !likelihood_)
          makeLikelihoodAtRoot_();
        
        return likelihood_;
      }

      ValueRef<Eigen::RowVectorXd> getSiteLikelihoods()
      {
        if (psites_ && !siteLikelihoods_)
          makeLikelihoodAtRoot_();
        
        return siteLikelihoods_;
      }

      double getLikelihoodForASite(size_t pos)
      {
        if (psites_ && !siteLikelihoods_)
            makeLikelihoodAtRoot_();

        return siteLikelihoods_->getTargetValue()[pos];
      }
      
      // ValueRef<double> getLikelihoodAtNode(uint nodeId) 
      // {
      //   makeLikelihoodsAtNode_(nodeId);
      // }
      
      void setData(const AlignedValuesContainer& sites)
      {
        psites_=&sites;
      }

      void setNumericalDerivateConfiguration(double delta, const NumericalDerivativeType& config);

      /**
       * Set Tree ClockLike :
       *  - add a RateNode parameter for multiplying all branch lengths
       *  - remove all branch lengths parameters from the parameters
       *
       */
      
      void setClockLike(double rate=1);

      const SubstitutionProcess& getSubstitutionProcess() const
      {
        return process_;
      }
        
      const AlignedValuesContainer* getData() const
      {
        return psites_;
      }

      std::size_t getNumberOfSites() const 
      {
        if (psites_!=0)
          return psites_->getNumberOfSites();
        else
          return 0;
      }
      
      const StateMap& getStateMap() const
      {
        return modelNode_->getTargetValue()->getStateMap();
      }
      
      std::shared_ptr<PhyloTree_BrRef> getTreeNode()
      {
        return treeNode_;
      }

      ValueRef<Eigen::RowVectorXd> getRootFreqs()
      {
        return rFreqs_;
      }

      std::shared_ptr<ForwardLikelihoodTree> getForwardTree(size_t nCat)
      {
        if (nCat>=vRateCatTrees_.size())
          throw Exception("LikelihoodCalculation::getForwardTree : Bad ForwardTree number " + TextTools::toString(nCat));
        return vRateCatTrees_[nCat].flt;
      }

      
    private:
      void makeForwardLikelihoodTree_();

      void makeProcessNodes_();

      void makeRootFreqs_();

      void makeLikelihoodAtRoot_();

      void makeLikelihoodsAtNode_(uint nodeId);

    };

  }
  
} // namespace bpp

#endif // LIKELIHOOD_CALCULATION

