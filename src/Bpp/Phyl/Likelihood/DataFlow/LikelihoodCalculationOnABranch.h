// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_LIKELIHOODCALCULATION_ON_A_BRANCH_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_LIKELIHOODCALCULATION_ON_A_BRANCH_H

#include <Bpp/Graph/AssociationDAGraphImplObserver.h>
#include <Bpp/Seq/Container/AlignmentData.h>

#include "Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h"
#include "Bpp/Phyl/Likelihood/DataFlow/Model.h"
#include "Bpp/Phyl/Likelihood/SubstitutionProcess.h"

namespace bpp
{
/** Likelihood computed on a branch. This likelihood is dependent on a
 *  LikelihoodCalculationSingleProcess, from which it will get the
 *  conditional likelihoods at both ends of the branch, a Branch Node
 *  (with the parameterized length), as well as a Model DF.
 *
 */


using LikelihoodFromRootConditionalAtRoot =
    MatrixProduct<RowLik, Eigen::RowVectorXd, MatrixLik>;

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

// Lower Conditional Likelihood under nodes
using ForwardLikelihoodBelow = Value<MatrixLik>;
using ForwardLikelihoodBelowRef = ValueRef<MatrixLik>;

// Upper Likelihood at top of edges
using BackwardLikelihoodAbove = Value<MatrixLik>;
using BackwardLikelihoodAboveRef = ValueRef<MatrixLik>;

using SiteLikelihoods = Value<RowLik>;
using SiteLikelihoodsRef = ValueRef<RowLik>;

using AllRatesSiteLikelihoods = MatrixLik;

using SiteWeights = NumericConstant<Eigen::RowVectorXi>;

using ForwardTransition =
    MatrixProduct<MatrixLik, Eigen::MatrixXd, MatrixLik>;

/*
 * @brief DAG of the conditional likelihoods (product of above and
 * below likelihoods), with same topology as forward & backward
 * likelihood DAGs.
 *
 */


class LikelihoodCalculationOnABranch :
  public AlignedLikelihoodCalculation
{
private:
  /* Considered Edge for each rate category */

  class RateCategoryEdge
  {
public:
    /* Vectors are necessary to account for all branches in the DAG matching a tree branch.
     */

    /* vector of Forward Conditional likelihoods at the bottom of the edge */
    std::vector<ForwardLikelihoodBelowRef> vBotLik_;

    /* Backward Conditional likelihood at the bottom of the edge */

    std::vector<BackwardLikelihoodAboveRef> vTopLik_;

    /* Conditional likelihoods of the edge */

    SiteLikelihoodsRef siteLik_;

    /************************************/
    /* Dependencies */

    std::shared_ptr<ConfiguredParameter> brlen_;
  };


  /*
   * @brief Number of sites
   *
   */

  size_t numberOfSites_;

  /*
   * @brief Dimension of the shrunked data
   */

  MatrixDimension likelihoodMatrixDim_;

  /*
   *@brief DF Model used on this branch
   *
   */

  std::shared_ptr<ConfiguredModel> model_;


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

  /*
   * @brief Node for the probabilites of the rate classes
   *
   */

  ValueRef<Eigen::RowVectorXd> catProb_;

  /* Likelihood Edges with for all rate categories */
  std::vector<RateCategoryEdge> vRateCatEdges_;

public:
  LikelihoodCalculationOnABranch(Context& context, LikelihoodCalculationSingleProcess& likcalsp, uint edgeId);


  void setModel(std::shared_ptr<ConfiguredModel> model)
  {
    if (model == model_)
      return;
    model_ = model;
    shareParameters_(model_->getParameters());
    makeLikelihoods();
  }

  LikelihoodCalculationOnABranch(const LikelihoodCalculationOnABranch& lik);

  LikelihoodCalculationOnABranch* clone() const
  {
    throw bpp::Exception("LikelihoodCalculationOnABranch clone should not happen.");
  }

  const StateMapInterface& stateMap() const
  {
    return model_->targetValue()->stateMap();
  }

  std::shared_ptr<const StateMapInterface> getStateMap() const
  {
    return model_->targetValue()->getStateMap();
  }

  std::shared_ptr<const Alphabet> getAlphabet() const
  {
    return stateMap().getAlphabet();
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
    return rootPatternLinks_ ? rootPatternLinks_->targetValue()(Eigen::Index(currentPosition)) : currentPosition;
  }

  const PatternType& getRootArrayPositions() const { return rootPatternLinks_->targetValue(); }

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
      return CWisePattern<RowLik>::create(getContext_(), {vector, rootPatternLinks_}, RowVectorDimension ((int)numberOfSites_));
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
      return CWisePattern<MatrixLik>::create(getContext_(), {matrix, rootPatternLinks_}, MatrixDimension (matrix->targetValue().rows(), Eigen::Index (numberOfSites_)));
  }

  /*
   * @brief: Get the weight of a position in the shrunked data (ie
   * the number of sites corresponding to this site)
   *
   */
  unsigned int getWeight(size_t pos) const
  {
    return (uint)(rootWeights_->targetValue()(Eigen::Index(pos)));
  }

  std::shared_ptr<SiteWeights> getRootWeights()
  {
    return rootWeights_;
  }

  /********************************************************
   * @Likelihoods
   *
   *****************************************************/

  /*********************************/
  /* @brief Methods for external usage (after lik computation) */

  /**
   * @brief Get site likelihoods for a rate category
   *
   * @param nCat : index of the rate category
   * @param shrunk : if returns on shrunked data (default: false)
   */
  RowLik getSiteLikelihoodsForAClass(size_t nCat, bool shrunk = false);

  /**
   * @brief Output array (Classes X Sites) of likelihoods for all
   * sites & classes.
   *
   * @param shrunk : if returns on shrunked data (default: false)
   */
  AllRatesSiteLikelihoods getSiteLikelihoodsForAllClasses(bool shrunk = false);

  size_t getNumberOfDistinctSites() const
  {
    return (size_t)likelihoodMatrixDim_.cols;
  }

  size_t getNumberOfSites() const
  {
    return numberOfSites_;
  }

  std::shared_ptr<ConfiguredModel> getModel()
  {
    return model_;
  }

private:
  void makeLikelihoods();

  std::shared_ptr<SiteLikelihoods> getSiteLikelihoods_(size_t nCat);

  void cleanAllLikelihoods();
};
} // namespace bpp

#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_LIKELIHOODCALCULATION_ON_A_BRANCH_H
