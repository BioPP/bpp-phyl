//
// File: LikelihoodCalculationSingleProcess.h
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

#ifndef LIKELIHOOD_CALCULATION_SINGLE_PROCESS_H
#define LIKELIHOOD_CALCULATION_SINGLE_PROCESS_H

#include "Bpp/Phyl/NewLikelihood/DataFlow/Model.h"
#include "Bpp/Phyl/NewLikelihood/DataFlow/DiscreteDistribution.h"
#include "Bpp/Phyl/NewLikelihood/DataFlow/FrequenciesSet.h"
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlow.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowCWise.h>

#include <Bpp/Seq/Container/AlignedValuesContainer.h>
#include <Bpp/Phyl/SitePatterns.h>
#include <Bpp/Graph/AssociationDAGraphImplObserver.h>
#include "Bpp/Phyl/NewLikelihood/SubstitutionProcess.h"


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
   * Initial conditional likelihood for leaves (sequences on the tree)
   * should be computed outside of the dataflow graph, and provided as
   * NumericConstant<MatrixXd>.
   *
   * Default Derivate Method is set to
   * NumericalDerivativeType::ThreePoints with delta=0.001;
   *  
   */


    class ProcessTree;
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
     * lik(state, site) = above(state, site) * below(state,site)
     * Using member wise multiply
     */

    using BuildConditionalLikelihood =
      CWiseMul<Eigen::MatrixXd, std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>>;


    
    using ConditionalLikelihood = Value<Eigen::MatrixXd>;
    
    using SiteLikelihoods = Value<Eigen::RowVectorXd>;
    using SiteLikelihoodsRef = ValueRef<Eigen::RowVectorXd>;
    
    using AllRatesSiteLikelihoods = Eigen::MatrixXd;
    
    using SiteWeights = NumericConstant<Eigen::RowVectorXi>;

    /*
     * @brief DAG of the conditional likelihoods (product of above and
     * below likelihoods), with same topology as forward & backward
     * likelihood DAGs.
     *
     */
     
    using ConditionalLikelihoodTree = AssociationDAGlobalGraphObserver<ConditionalLikelihood,uint>;

    using SiteLikelihoodsTree = AssociationTreeGlobalGraphObserver<SiteLikelihoods, uint>;

    class LikelihoodCalculationSingleProcess :
      public AbstractParametrizable
    {
    private:
      
      class RateCategoryTrees {
      public:
        std::shared_ptr<ProcessTree> phyloTree;
        std::shared_ptr<ForwardLikelihoodTree> flt;

        /*
         * @brief backward likelihood tree (only computed when needed)
         *
         */
        
        std::shared_ptr<BackwardLikelihoodTree> blt;

        /*
         * @brief for each node n:  clt[n] = flt[n] * dlt[n]
         *
         */
        std::shared_ptr<ConditionalLikelihoodTree> clt;
        
        /*
         * Site Likelihoods on the tree, with shrunked positions.
         * 
         * @brief for each node n:  lt[n] = sum_{state s} flt[n][s]
         *
         * This is only computed when node specific likelihood is
         * computed.
         */

        std::shared_ptr<SiteLikelihoodsTree> lt;
      };

      class ProcessNodes {
      public:
        std::shared_ptr<ProcessTree> treeNode_;
        std::shared_ptr<ConfiguredModel> modelNode_; // Used for
                                                     // StateMap and
                                                     // root
                                                     // frequencies if
                                                     // stationarity
        std::shared_ptr<ConfiguredFrequenciesSet> rootFreqsNode_;
        std::shared_ptr<ConfiguredDistribution> ratesNode_;
      };
        
      Context& context_;

      /************************************/
      /* Dependencies */
      
      const SubstitutionProcess& process_;
      const AlignedValuesContainer* psites_;
      
      /*****************************
       ****** Patterns
       *
       * @brief Links between sites and patterns.
       * 
       * The size of this vector is equal to the number of sites in the container,
       * each element corresponds to a site in the container and points to the
       * corresponding column in the likelihood array of the root node.
       * If the container contains no repeated site, there will be a strict
       * equivalence between each site and the likelihood array of the root node.
       * However, if this is not the case, some pointers may point toward the same
       * element in the likelihood array.
       */
      
      ValueRef<PatternType> rootPatternLinks_;

      /**
       * @brief The frequency of each site.
       */

      std::shared_ptr<SiteWeights> rootWeights_;
      std::shared_ptr<AlignedValuesContainer> shrunkData_;

      /************************************/
      /* DataFlow objects */

      ProcessNodes processNodes_;

      ValueRef<Eigen::RowVectorXd> rFreqs_;

      /******************************************/
      /** Likelihoods  **/
          
      ValueRef<double> likelihood_;

      ValueRef<Eigen::RowVectorXd> siteLikelihoods_;

      ValueRef<Eigen::RowVectorXd> patternedSiteLikelihoods_;
      
      /* Likelihood Trees with for all rate categories */
      std::vector<RateCategoryTrees> vRateCatTrees_;

    public:
      LikelihoodCalculationSingleProcess(Context & context,
                                         const AlignedValuesContainer & sites,
                                         const SubstitutionProcess& process);

      LikelihoodCalculationSingleProcess(Context & context,
                                         const SubstitutionProcess& process);

      LikelihoodCalculationSingleProcess(Context & context,
                                         const AlignedValuesContainer & sites,
                                         const SubstitutionProcess& process,
                                         ParameterList& paramList);

      LikelihoodCalculationSingleProcess(const LikelihoodCalculationSingleProcess& lik);

      LikelihoodCalculationSingleProcess* clone() const
      {
        return new LikelihoodCalculationSingleProcess(*this);
      }
    
      ValueRef<double> getLikelihood() 
      {
        if (!shrunkData_)
          throw Exception("LikelihoodCalculationSingleProcess::getLikelihood : data not set.");
        
        if (!likelihood_)
          makeLikelihoodAtRoot_();
        
        return likelihood_;
      }

      double getLogLikelihood() 
      {
        return getLikelihood()->getTargetValue();
      }

      size_t getNumberOfSites() const
      {
        return (psites_?psites_->getNumberOfSites():0);
      }

      size_t getNumberOfDistinctSites() const
      {
        return shrunkData_?shrunkData_->getNumberOfSites():0;
      }

      
      /*
       * @brief Return the Sitelikelihoods_ vector on shrunked data
       *
       * @brief shrunk: bool true if vector on shrunked data
       *
       */
      
      ValueRef<Eigen::RowVectorXd> getSiteLikelihoods(bool shrunk)
      {
        if (!shrunkData_)
          throw Exception("LikelihoodCalculationSingleProcess::getSiteLikelihoods : data not set.");
        
        if (!siteLikelihoods_)
          makeLikelihoodAtRoot_();

        if (shrunk)
          return siteLikelihoods_;
        else
          return patternedSiteLikelihoods_;
      }

      /************************************************
       *** Patterns
       ****************************/
      
      /*
       * @brief the relations between real position and shrunked data
       * positions.
       *
       */
     
      size_t getRootArrayPosition(size_t currentPosition) const
      {
        return rootPatternLinks_->getTargetValue()(currentPosition);
      }
    
      const PatternType& getRootArrayPositions() const { return rootPatternLinks_->getTargetValue(); }

      const AlignedValuesContainer* getShrunkData() const {
        return shrunkData_.get();
      }

      /*
       * @brief Expands (ie reverse of shrunkage) a vector computed of
       * shrunked data (ie from one value per distinct site to one
       * value per site).
       *
       */
      
      ValueRef<Eigen::RowVectorXd> expandVector(ValueRef<Eigen::RowVectorXd> vector)
      {
        return CWisePattern<Eigen::RowVectorXd>::create(context_,{vector,rootPatternLinks_}, rowVectorDimension (Eigen::Index (getData()->getNumberOfSites())));
      }
      
      /*
       * @brief: Get the weight of a position in the shrunked data (ie
       * the number of sites corresponding to this site)
       *
       */
      
      unsigned int getWeight(size_t pos) const
      {
        return rootWeights_->getTargetValue()(pos);
      }

      ValueRef<double> getLikelihoodAtNode(uint nodeId) 
      {
        const auto pt=process_.getParametrizablePhyloTree();
        if (!pt.hasNode(nodeId))
          throw Exception("LikelihoodCalculationSingleProcess::getLikelihoodAtNode : node " + TextTools::toString(nodeId) + " does not exist in Phylo Tree.");

        // No need to recompute with backward if at root
        if (pt.getRootIndex()==nodeId)
          return likelihood_;

        return makeLikelihoodsAtNode_(nodeId);
      }
      
      void setData(const AlignedValuesContainer& sites)
      {
        if (psites_)
          throw Exception("LikelihoodCalculationSingleProcess::setData : data already assigned.");
        psites_=&sites;
        setPatterns_();
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

      bool isInitialized() const {
        return getData();
      };


      const StateMap& getStateMap() const
      {
        return processNodes_.modelNode_->getTargetValue()->getStateMap();
      }
      
      std::shared_ptr<ProcessTree> getTreeNode()
      {
        return processNodes_.treeNode_;
      }

      ValueRef<Eigen::RowVectorXd> getRootFreqs()
      {
        return rFreqs_;
      }

      //std::shared_ptr<ConditionalLikelihoodTree> getConditionalLikelihoodTree(size_t nCat);

      /*
       *@brief Get tree of site likelihoods, on shrunk data.
       *
       */
      
      std::shared_ptr<SiteLikelihoodsTree> getSiteLikelihoodsTree(size_t nCat);
      
      /*
       *@ brief make DF nodes of the process, using
       *ConfiguredParameters defined in a ParameterList.
       */

      void makeProcessNodes(ParameterList& pl);

            
      /*********************************/
      /*@brief Methods for external usage (after lik computation) */

      /*
       *@brief Get site likelihoods for a rate category
       *
       *@param nCat : index of the rate category
       *@param shrunk : if returns on shrunked data (default: false)
       */
      
      
            
      const Eigen::RowVectorXd& getSiteLikelihoodsForAClass(size_t nCat, bool shrunk = false);

      /*
       *@brief Output array (Classes X Sites) of likelihoods for all
       *sites & classes.
       *
       *@param shrunk : if returns on shrunked data (default: false)
       */

      const AllRatesSiteLikelihoods& getSiteLikelihoodsForAllClasses(bool shrunk = false);
      
      double getLikelihoodForASite(size_t pos)
      {
        if (!shrunkData_)
          throw Exception("LikelihoodCalculationSingleProcess::getLikelihoodForASite : data not set.");
        
        if (!siteLikelihoods_)
          makeLikelihoodAtRoot_();

        return patternedSiteLikelihoods_->getTargetValue()(pos);
      }

    private:
      void setPatterns_();
      
      void makeForwardLikelihoodTree_();

      void makeProcessNodes_();
      
      void makeRootFreqs_();

      void makeLikelihoodAtRoot_();

      /*
       * @brief Compute the likelihood at a given node in the tree,
       * which number may not be the same number number in the DAG.
       * Several node in the DAG may be related to this tree node, in
       * which case a sum is computer.
       *
       * @param nodeId : index of the node in the phylo Tree
       */
      
      ValueRef<double> makeLikelihoodsAtNode_(uint nodeId);

    };

} // namespace bpp

#endif // LIKELIHOOD_CALCULATION_SINGLE_PROCESS_H

