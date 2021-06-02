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

#include "Bpp/Phyl/NewLikelihood/DataFlow/LikelihoodCalculation.h"

#include "Bpp/Phyl/NewLikelihood/DataFlow/CollectionNodes.h"
#include "Bpp/Phyl/NewLikelihood/DataFlow/Model.h"
#include "Bpp/Phyl/NewLikelihood/DataFlow/DiscreteDistribution.h"
#include "Bpp/Phyl/NewLikelihood/DataFlow/FrequencySet.h"
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlow.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowCWiseComputing.h>

#include "Bpp/Phyl/NewLikelihood/SubstitutionProcess.h"

#include <Bpp/Seq/Container/AlignedValuesContainer.h>
#include <Bpp/Phyl/SitePatterns.h>
#include <Bpp/Graph/AssociationDAGraphImplObserver.h>


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
   * NumericalDerivativeType::ThreePoints with delta=0.0001;
   *  
   */


  class ProcessTree;
  class ForwardLikelihoodTree;
  class BackwardLikelihoodTree;

  using RowLik = Eigen::Matrix<double, 1, Eigen::Dynamic>;
  using MatrixLik = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  
  /** @brief likelihood = f(equilibriumFrequencies, rootConditionalLikelihood).
   * - likelihood: RowVector(site).
   * - equilibriumFrequencies: RowVector(state).
   * - rootConditionalLikelihood: Matrix(state, site).
   *
   * likelihood(site) = sum_state equFreqs(state) * rootConditionalLikelihood(state, site).
   * Using matrix multiply: likelihood = equilibriumFrequencies * rootConditionalLikelihood.
   */

  using LikelihoodFromRootConditional =
    MatrixProduct<RowLik, RowLik, MatrixLik>;
    
  /** @brief totalLikelihood = product_site likelihood(site).
   * - likelihood: RowVector (site).
   * - totalLikelihood: Extended float.
   */
    
  using TotalLogLikelihood = SumOfLogarithms<RowLik>;

  /** @brief Conditionallikelihood = AboveConditionalLikelihood * BelowConditionalLikelihood
   *
   * lik(state, site) = above(state, site) * below(state,site)
   * Using member wise multiply
   */

  using BuildConditionalLikelihood =
    CWiseMul<MatrixLik, std::tuple<MatrixLik, MatrixLik>>;
    
  using ConditionalLikelihood = Value<MatrixLik>;
  using ConditionalLikelihoodRef = ValueRef<MatrixLik>;
    
  using SiteLikelihoods = Value<RowLik>;
  using SiteLikelihoodsRef = ValueRef<RowLik>;
    
  using AllRatesSiteLikelihoods = MatrixLik;
    
  using SiteWeights = NumericConstant<Eigen::RowVectorXi>;

  /*
   * @brief DAG of the conditional likelihoods (product of above and
   * below likelihoods), with same topology as forward & backward
   * likelihood DAGs.
   *
   */
     
  using ConditionalLikelihoodDAG = AssociationDAGlobalGraphObserver<ConditionalLikelihood,uint>;

  using ConditionalLikelihoodTree = AssociationTreeGlobalGraphObserver<ConditionalLikelihood,uint>;

  using SiteLikelihoodsDAG = AssociationDAGlobalGraphObserver<SiteLikelihoods, uint>;

  using SiteLikelihoodsTree = AssociationTreeGlobalGraphObserver<SiteLikelihoods, uint>;

  using DAGindexes = std::vector<uint>;

  class LikelihoodCalculationSingleProcess :
    public AlignedLikelihoodCalculation
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
      std::shared_ptr<ConditionalLikelihoodDAG> clt;

      /*
       * @brief for each node n:  clt[n] = sum_states(clt[n][s])
       *
       */
      
      std::shared_ptr<SiteLikelihoodsDAG> lt;

      /*
       * Site Likelihoods on the tree, with shrunked positions, summed
       * on all paths.
       * 
       * @brief for each node n:  lt[n] = sum_{state s} flt[n][s]
       *
       * This is only computed when node specific likelihood is
       * computed.
       */

      std::shared_ptr<SiteLikelihoodsTree> speciesLt;

    };

    /*
     * @brief DF Nodes used in the process. ProcessTree is used
     * without any rate multiplier.
     *
     */
    
    class ProcessNodes {
    public:
      std::shared_ptr<ProcessTree> treeNode_;
      std::shared_ptr<ConfiguredModel> modelNode_; // Used for
      // StateMap and root frequencies if stationarity
      
      std::shared_ptr<ConfiguredFrequencySet> rootFreqsNode_;
      std::shared_ptr<ConfiguredDistribution> ratesNode_;
    };
        
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

    /* Likelihood Trees with for all rate categories */
    std::vector<RateCategoryTrees> vRateCatTrees_;

    /* Likelihood tree on mean likelihoods on rate categories */
    std::shared_ptr<ConditionalLikelihoodTree> condLikelihoodTree_;
      
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

    LikelihoodCalculationSingleProcess(Context & context,
                                       const SubstitutionProcess& process,
                                       ParameterList& paramList);

    /*
     * @brief Build using Nodes of CollectionNodes.
     *
     * @param collection The CollectionNodes
     * @param nProcess the process Number in the collection
     * @param nData the data Number in the collection
     */

    LikelihoodCalculationSingleProcess(CollectionNodes& collection,
                                       const AlignedValuesContainer & sites,
                                       size_t nProcess);

    
    LikelihoodCalculationSingleProcess(CollectionNodes& collection,
                                       size_t nProcess);

    
    LikelihoodCalculationSingleProcess(const LikelihoodCalculationSingleProcess& lik);

    LikelihoodCalculationSingleProcess* clone() const
    {
      return new LikelihoodCalculationSingleProcess(*this);
    }

    void setData(const AlignedValuesContainer& sites)
    {
      if (psites_)
        throw Exception("LikelihoodCalculationSingleProcess::setData : data already assigned.");
      psites_=&sites;
      setPatterns_();
    }

    /**
     * @brief Set derivation procedure (see DataFlowNumeric.h)
     *
     */
      
    void setNumericalDerivateConfiguration(double delta, const NumericalDerivativeType& config);

    /**
     * Set Tree ClockLike :
     *  - add a RateNode parameter for multiplying all branch lengths
     *  - remove all branch lengths parameters from the parameters
     *
     */
      
    void setClockLike(double rate=1);

    /**************************************************/

    /*
     * @brief Return the loglikehood (see as a function, ie
     * computation node).
     *
     */
      
    void makeLikelihoods()
    {
      if (!psites_)
        throw Exception("LikelihoodCalculationSingleProcess::makeLikelihoods : data not set.");
        
      makeLikelihoodsAtRoot_();
    }


    /*
     * @brief Get indexes of the nodes in the Likelihood DAG that have
     * a given species index.
     *
     * @param speciesId  Looked species Index
     *
     */
    
    const DAGindexes& getNodesIds(uint speciesId) const;

    /*
     * @brief Get indexes of the non-empty edges in the Likelihood DAG
     * that have a given species index for a given rate class index.
     *
     * @param speciesId  Looked species Index
     * @param nCat  Rate class category
     *
     */
    
    const DAGindexes& getEdgesIds(uint speciesId, size_t nCat) const;

    /*
     * @brief Return the loglikehood computation at a given node.
     *
     */

    // ValueRef<double> getLogLikelihoodAtNode(uint nodeId) 
    // {
    //   const auto pt=process_.getParametrizablePhyloTree();
    //   if (!pt.hasNode(nodeId))
    //     throw Exception("LikelihoodCalculationSingleProcess::getLogLikelihoodAtNode : node " + TextTools::toString(nodeId) + " does not exist in Phylo Tree.");

    //   // No need to recompute with backward if at root
    //   if (pt.getRootIndex()==nodeId)
    //     return likelihood_;

    //   makeLikelihoodsAtNode_(nodeId);

      
    // }

      
    size_t getNumberOfSites() const
    {
      return (psites_?psites_->getNumberOfSites():0);
    }

    size_t getNumberOfDistinctSites() const
    {
      return shrunkData_?shrunkData_->getNumberOfSites():getNumberOfSites();
    }
 
    /*
     * @brief Return the ref to the SubstitutionProcess
     *
     * Warning; the process parameter values may not be up to date
     * with some of the LikelihoodCalculationSingleProcess
     *
     */
    
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
      
    /************************************************
     *** Patterns
     ****************************/
      
    /*
     * @brief the relations between real position and shrunked data
     * positions.
     *
     * @param currentPosition : position in real data
     *
     * @return matching position in shrunked data
     *
     */
     
    size_t getRootArrayPosition(size_t currentPosition) const
    {
      return rootPatternLinks_?rootPatternLinks_->getTargetValue()(Eigen::Index(currentPosition)):currentPosition;
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
      
    ValueRef<RowLik> expandVector(ValueRef<RowLik> vector)
    {
      if (!rootPatternLinks_)
        return vector;
      else
        return CWisePattern<RowLik>::create(getContext_(),{vector,rootPatternLinks_}, RowVectorDimension (Eigen::Index (getData()->getNumberOfSites())));
    }

    /*
     * @brief Expands (ie reverse of shrunkage) a matrix computed of
     * shrunked data (ie from one value per distinct site to one
     * value per site). Columns are sites.
     *
     */
      
    ValueRef<MatrixLik> expandMatrix(ValueRef<MatrixLik> matrix)
    {
      if (!rootPatternLinks_)
        return matrix;
      else
        return CWisePattern<MatrixLik>::create(getContext_(),{matrix,rootPatternLinks_}, MatrixDimension (matrix->getTargetValue().rows(), Eigen::Index (getData()->getNumberOfSites())));
    }

    /*
     * @brief: Get the weight of a position in the shrunked data (ie
     * the number of sites corresponding to this site)
     *
     */
      
    unsigned int getWeight(size_t pos) const
    {
      return (uint)(rootWeights_->getTargetValue()(Eigen::Index(pos)));
    }

    std::shared_ptr<SiteWeights> getRootWeights()
    {
      return rootWeights_;
    }

    ValueRef<Eigen::RowVectorXd> getRootFreqs()
    {
      return rFreqs_;
    }

    /********************************************************
     * @Likelihoods
     *
     *****************************************************/

    /*
     * @brief Get Matrix of Conditional Likelihoods at Node *
     *
     * @param nodeId  Id of the node in PhyloTree, ie species id
     * @param shrunk if matrix is on shrunked data (default: false)
     *
     */
    
    ConditionalLikelihoodRef getLikelihoodsAtNode(uint nodeId, bool shrunk = false)
    {
      if (!(condLikelihoodTree_ && condLikelihoodTree_->hasNode(nodeId)))
        makeLikelihoodsAtNode_(nodeId);

      auto vv = condLikelihoodTree_->getNode(nodeId);

      return shrunk?vv:expandMatrix(vv);
    }

    /*
     * @brief Get forward shrunked likelihood matrix at Node (ie just
     * above the node), for a given rate class.
     *
     * @param nodeId Node Index in the forward tree (! ie in the
     * computation tree, not the species tree).
     *
     * @param nCat  Rate class category
     *
     */
    
    ConditionalLikelihoodRef getForwardLikelihoodsAtNodeForClass(uint nodeId, size_t nCat);

    /*
     * @brief Get backward shrunked likelihood matrix at Edge (ie at
     * the top of the edge), for a given rate class.
     *
     * @param edgeId Edge Index in the backward tree (! ie in the
     * computation tree, not the species tree).
     *
     * @param nCat  Rate class category
     *
     */
    
    ConditionalLikelihoodRef getBackwardLikelihoodsAtNodeForClass(uint nodeId, size_t nCat);

    /*
     * @brief Get backward shrunked likelihood matrix at Node (ie at
     * the top of the edge), for a given rate class.
     *
     * These likelihoods are multiplied by the probability of the node
     *
     * @param edgeId Node Index in the backward tree (! ie in the
     * computation tree, not the species tree).
     *
     * @param nCat  Rate class category
     *
     */
    
    ConditionalLikelihoodRef getBackwardLikelihoodsAtEdgeForClass(uint edgeId, size_t nCat);

    /*
     * @brief Get shrunked conditional likelihood matrix at Node (ie
     * just above the node), for a given rate class.
     *
     * These likelihoods are multiplied by the probability of the node.
     *
     * @param nodeId Node Index in the forward tree (! ie in the
     * computation tree, not the species tree).
     *
     * @param nCat  Rate class category
     *
     */
    
    ConditionalLikelihoodRef getConditionalLikelihoodsAtNodeForClass(uint nodeId, size_t nCat);

    /*
     * @brief Get shrunked conditional likelihood matrix at Node (ie
     * just above the node), for a given rate class.
     *
     * These likelihoods are multiplied by the probability of the node.
     *
     * @param nodeId Node Index in the forward tree (! ie in the
     * computation tree, not the species tree).
     *
     * @param nCat  Rate class category
     *
     */

    SiteLikelihoodsRef getLikelihoodsAtNodeForClass(uint nodeId, size_t nCat);
    
    /*
     * @brief make backward likelihood tree (only computed when needed)
     *
     */
        
    void makeLikelihoodsTree()
    {
      auto allIndex = process_.getParametrizablePhyloTree().getAllNodesIndexes();
      
      for (auto id: allIndex)
         makeLikelihoodsAtNode_(id);
    }

    /*********************************/
    /*@brief Methods for external usage (after lik computation) */

    /*
     *@brief Get site likelihoods for a rate category
     *
     *@param nCat : index of the rate category
     *@param shrunk : if returns on shrunked data (default: false)
     */
      
    RowLik getSiteLikelihoodsForAClass(size_t nCat, bool shrunk = false);

    /*
     *@brief Output array (Classes X Sites) of likelihoods for all
     *sites & classes.
     *
     *@param shrunk : if returns on shrunked data (default: false)
     */

    AllRatesSiteLikelihoods getSiteLikelihoodsForAllClasses(bool shrunk = false);

    
    /*
     *@brief Get process tree for a rate category
     *
     *@param nCat : index of the rate category
     *
     */
    
    std::shared_ptr<ProcessTree> getTreeNode(size_t nCat)
    {
      if (nCat >= vRateCatTrees_.size())
        throw Exception("LikelihoodCalculationSingleProcess::getTreeNode : bad class number " + TextTools::toString(nCat));
      
      return vRateCatTrees_[nCat].phyloTree;
    }

    std::shared_ptr<ForwardLikelihoodTree> getForwardLikelihoodTree(size_t nCat);

  private:
    void setPatterns_();
      
    void makeForwardLikelihoodTree_();

    void makeProcessNodes_();
      
    void makeRootFreqs_();

    void makeLikelihoodsAtRoot_();

    /*
     *@ brief make DF nodes of the process, using
     *ConfiguredParameters defined in a ParameterList.
     */

    void makeProcessNodes_(ParameterList& pl);

    /*
     *@ brief make DF nodes of a process in a collection, using
     *ConfiguredParameters defined in a CollectionNodes.
     */

    void makeProcessNodes_(CollectionNodes& pl, size_t nProc);

    /*
     * @brief Compute the likelihood at a given node in the tree,
     * which number may not be the same number number in the DAG.
     *
     * Several nodes in the DAG may be related to this tree node, in
     * which case a sum is computed.
     *
     * @param nodeId : index of the node in the phylo Tree
     */
      
    void makeLikelihoodsAtNode_(uint nodeId);

    /*
     * @brief Compute the likelihood at a given node in the DAG,
     *
     * This is not enough to compute likelihoods at species nodes, use
     * makeLikelihoodsAtNode_ instead.
     *
     * @param nodeId : index of the node in the DAG
     */
      
    void makeLikelihoodsAtDAGNode_(uint nodeId);

    std::shared_ptr<SiteLikelihoodsTree> getSiteLikelihoodsTree_(size_t nCat);

  };

} // namespace bpp

#endif // LIKELIHOOD_CALCULATION_SINGLE_PROCESS_H

