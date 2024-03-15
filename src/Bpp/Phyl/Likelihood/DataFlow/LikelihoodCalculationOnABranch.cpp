// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <list>
#include <numeric>
#include <unordered_map>

#include "Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationOnABranch.h"
#include "Bpp/Phyl/Likelihood/DataFlow/ForwardLikelihoodTree.h"
#include "Bpp/Phyl/Likelihood/DataFlow/BackwardLikelihoodTree.h"
#include "Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h"

using namespace std;
using namespace bpp;
using namespace numeric;
using namespace Eigen;

LikelihoodCalculationOnABranch::LikelihoodCalculationOnABranch(
  Context& context,
  LikelihoodCalculationSingleProcess& likcalsp,
  uint edgeId) :
  AlignedLikelihoodCalculation(context),
  numberOfSites_(likcalsp.getNumberOfSites()),
  likelihoodMatrixDim_({likcalsp.stateMap().getNumberOfModelStates(), likcalsp.getNumberOfDistinctSites()}),
  model_(likcalsp.getNumberOfClasses()?likcalsp.getTreeNode(0)->getEdge(edgeId)->getModel():0),
  rootPatternLinks_(likcalsp.getRootPatternLinks()),
  rootWeights_(likcalsp.getRootWeights()),
  catProb_(likcalsp.catProb_), vRateCatEdges_()
{
  auto nbCl = likcalsp.getNumberOfClasses();
  // ensure that all likelihoods exist
  likcalsp.makeLikelihoodsTree();
  for (size_t ncl = 0; ncl<nbCl;ncl++)
  {
    auto tr = likcalsp.getTreeNode(ncl);

    const auto& dagId = likcalsp.getEdgesIds(edgeId, ncl);
    auto catEdge = make_shared<RateCategoryEdge>();    
    for (const auto& did:dagId)
    {
      auto ft = likcalsp.getForwardLikelihoodTree(ncl);
      catEdge->vBotLik_.push_back(ft->getNode(ft->getSon(did)));
      catEdge->vTopLik_.push_back(likcalsp.getBackwardLikelihoodTree(ncl)->getEdge(did));
    }

    catEdge->brlen_ = likcalsp.getTreeNode(ncl)->getEdge(dagId[0])->getBrLen();
    vRateCatEdges_.push_back(*catEdge);
  }
  shareParameters_(likcalsp.getParameters());
}

LikelihoodCalculationOnABranch::LikelihoodCalculationOnABranch(const LikelihoodCalculationOnABranch& lik) :
  AlignedLikelihoodCalculation(lik),
  numberOfSites_(lik.numberOfSites_),
  likelihoodMatrixDim_(lik.likelihoodMatrixDim_),
  model_(lik.model_),
  rootPatternLinks_(lik.rootPatternLinks_),
  rootWeights_(lik.rootWeights_),
  catProb_(lik.catProb_),
  vRateCatEdges_(lik.vRateCatEdges_)
{
}

RowLik LikelihoodCalculationOnABranch::getSiteLikelihoodsForAClass(size_t nCat, bool shrunk)
{
  if (shrunk)
    return vRateCatEdges_[nCat].siteLik_->targetValue();
  else
    return expandVector(vRateCatEdges_[nCat].siteLik_)->targetValue();
}

AllRatesSiteLikelihoods LikelihoodCalculationOnABranch::getSiteLikelihoodsForAllClasses(bool shrunk)
{
  auto nbCat = vRateCatEdges_.size();
  auto allLk = std::make_shared<AllRatesSiteLikelihoods>(nbCat, shrunk ? getNumberOfDistinctSites() : numberOfSites_);

  for (size_t nCat = 0; nCat < nbCat; nCat++)
  {
    allLk->row(Eigen::Index(nCat)) = getSiteLikelihoodsForAClass(nCat, shrunk);
  }
 
  return *allLk;
}


/****************************************
* Construction methods
****************************************/

void LikelihoodCalculationOnABranch::makeLikelihoods()
{
  if (!model_)
    throw Exception("LikelihoodCalculationOnABranch::makeLikelihoods_: Missing model");

  auto nbDistSite = Eigen::Index(getNumberOfDistinctSites());
  auto nbState = Eigen::Index(stateMap().getNumberOfModelStates());

  std::shared_ptr<ConditionalLikelihood> cond(0);

  SiteLikelihoodsRef distinctSiteLikelihoodsNode;

  auto one = ConstantOne<Eigen::RowVectorXd>::create(getContext_(), RowVectorDimension (nbState));
  
  const auto model = getModel();

  auto zero = getContext_().getZero();

  std::vector<std::shared_ptr<Node_DF> > vCond;

  std::vector<NodeRef> vLikRate;

  for (auto& rateCat: vRateCatEdges_)
  {
    auto transitionMatrix = ConfiguredParametrizable::createMatrix<ConfiguredModel, TransitionMatrixFromModel, Eigen::MatrixXd>(getContext_(), {model, rateCat.brlen_, zero}, transitionMatrixDimension(static_cast<size_t>(nbState)));

    for (size_t i = 0; i < rateCat.vBotLik_.size(); i++)
    {
      auto forwardEdge = ForwardTransition::create (
        getContext_(), {transitionMatrix, rateCat.vBotLik_[i]}, likelihoodMatrixDim_);

      cond = BuildConditionalLikelihood::create (
        getContext_(), {rateCat.vTopLik_[i], forwardEdge}, likelihoodMatrixDim_);

      if (rateCat.vBotLik_.size()>1)
        vCond.push_back(cond);
    }

    /*
     * If several DAG nodes related with this species node, sum the
     * likelihoods of all (already multiplied by their probability).
     *
     */


    if (vCond.size()>1)
      cond = CWiseAdd<MatrixLik, ReductionOf<MatrixLik> >::create(getContext_(), std::move(vCond), likelihoodMatrixDim_);

    rateCat.siteLik_ = LikelihoodFromRootConditionalAtRoot::create (
      getContext_(), {one, cond}, RowVectorDimension (nbDistSite));

    vLikRate.push_back(rateCat.siteLik_);
  }

  if (catProb_)
  {
    auto catProb = Convert<RowLik, Eigen::RowVectorXd>::create(getContext_(), {catProb_}, RowVectorDimension (Eigen::Index (vRateCatEdges_.size())));
    vLikRate.push_back(catProb);
    distinctSiteLikelihoodsNode = CWiseMean<RowLik, ReductionOf<RowLik>, RowLik>::create(getContext_(), std::move(vLikRate), RowVectorDimension (nbDistSite));
  }
  else
    distinctSiteLikelihoodsNode = vRateCatEdges_[0].siteLik_;

  
  // likelihoods per distinct site
  setSiteLikelihoods(distinctSiteLikelihoodsNode, true);
  
  // likelihoods per site
  setSiteLikelihoods(expandVector(patternedSiteLikelihoods_), false);
  
  // global likelihood
  ValueRef<DataLik> val;
  if (rootPatternLinks_)
    val = TotalLogLikelihood::create (getContext_(), {distinctSiteLikelihoodsNode, rootWeights_}, RowVectorDimension (Eigen::Index (nbDistSite)));
  else
    val = TotalLogLikelihood::create (getContext_(), {distinctSiteLikelihoodsNode}, RowVectorDimension (Eigen::Index (nbDistSite)));

  setLikelihoodNode(val);
}


void LikelihoodCalculationOnABranch::cleanAllLikelihoods()
{
  rootPatternLinks_.reset();
  rootWeights_.reset();

  vRateCatEdges_.clear();

  AlignedLikelihoodCalculation::cleanAllLikelihoods();
}
   
