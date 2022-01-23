//
// File: LikelihoodCalculationSingleProcess.cpp
// Authors:
//   FranÃÂ§ois Gindraud, Laurent GuÃÂ©guen (2018)
// Created: jeudi 28 fÃÂ©vrier 2019, ÃÂ  07h 22
//

#include <Bpp/Phyl/Model/MixedTransitionModel.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <list>
#include <numeric>
#include <unordered_map>

#include "Bpp/Phyl/Likelihood/DataFlow/BackwardLikelihoodTree.h"
#include "Bpp/Phyl/Likelihood/DataFlow/ForwardLikelihoodTree.h"
#include "Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h"
#include "Bpp/Phyl/Likelihood/SubstitutionProcessCollectionMember.h"

using namespace std;
using namespace bpp;
using namespace numeric;
using namespace Eigen;

LikelihoodCalculationSingleProcess::LikelihoodCalculationSingleProcess(Context& context,
                                                                       const AlignedValuesContainer& sites,
                                                                       const SubstitutionProcess& process) :
  AlignedLikelihoodCalculation(context), process_(process), psites_(&sites),
  rootPatternLinks_(), rootWeights_(), shrunkData_(),
  processNodes_(), rFreqs_(),
  vRateCatTrees_(), condLikelihoodTree_(0)
{
  setPatterns_();
  makeProcessNodes_();

  // Default Derivate
  setNumericalDerivateConfiguration(0.0001, NumericalDerivativeType::ThreePoints);
}

LikelihoodCalculationSingleProcess::LikelihoodCalculationSingleProcess(Context& context,
                                                                       const SubstitutionProcess& process) :
  AlignedLikelihoodCalculation(context),
  process_(process), psites_(),
  rootPatternLinks_(), rootWeights_(), shrunkData_(),
  processNodes_(), rFreqs_(),
  vRateCatTrees_(), condLikelihoodTree_(0)
{
  makeProcessNodes_();

  // Default Derivate
  setNumericalDerivateConfiguration(0.0001, NumericalDerivativeType::ThreePoints);
}

LikelihoodCalculationSingleProcess::LikelihoodCalculationSingleProcess(CollectionNodes& collection,
                                                                       const AlignedValuesContainer& sites,
                                                                       size_t nProcess) :
  AlignedLikelihoodCalculation(collection.getContext()), process_(collection.getCollection().getSubstitutionProcess(nProcess)), psites_(&sites),
  rootPatternLinks_(), rootWeights_(), shrunkData_(),
  processNodes_(), rFreqs_(),
  vRateCatTrees_(), condLikelihoodTree_(0)
{
  setPatterns_();
  makeProcessNodes_(collection, nProcess);

  // Default Derivate
  setNumericalDerivateConfiguration(0.0001, NumericalDerivativeType::ThreePoints);
}


LikelihoodCalculationSingleProcess::LikelihoodCalculationSingleProcess(CollectionNodes& collection,
                                                                       size_t nProcess) :
  AlignedLikelihoodCalculation(collection.getContext()), process_(collection.getCollection().getSubstitutionProcess(nProcess)), psites_(),
  rootPatternLinks_(), rootWeights_(), shrunkData_(),
  processNodes_(), rFreqs_(),
  vRateCatTrees_(), condLikelihoodTree_(0)
{
  makeProcessNodes_(collection, nProcess);

  // Default Derivate
  setNumericalDerivateConfiguration(0.0001, NumericalDerivativeType::ThreePoints);
}


LikelihoodCalculationSingleProcess::LikelihoodCalculationSingleProcess(const LikelihoodCalculationSingleProcess& lik) :
  AlignedLikelihoodCalculation(lik),
  process_(lik.process_), psites_(lik.psites_),
  rootPatternLinks_(lik.rootPatternLinks_), rootWeights_(), shrunkData_(lik.shrunkData_),
  processNodes_(), rFreqs_(),
  vRateCatTrees_(), condLikelihoodTree_(0)
{
  setPatterns_();
  makeProcessNodes_();

  // Default Derivate
  setNumericalDerivateConfiguration(0.0001, NumericalDerivativeType::ThreePoints);
}

void LikelihoodCalculationSingleProcess::setPatterns_()
{
  SitePatterns patterns(psites_, process_.getParametrizablePhyloTree().getAllLeavesNames());
  shrunkData_.reset(patterns.getSites());
  rootPatternLinks_ = NumericConstant<PatternType>::create(getContext_(), patterns.getIndices());
  size_t nbSites    = shrunkData_->getNumberOfSites();
  Eigen::RowVectorXi weights(nbSites);
  for (std::size_t i = 0; i < nbSites; i++)
  {
    weights(Eigen::Index(i)) = int(patterns.getWeights()[i]);
  }
  rootWeights_ = SiteWeights::create(getContext_(), std::move(weights));
}

void LikelihoodCalculationSingleProcess::makeProcessNodes_()
{
#ifdef DEBUG
  cerr << "LikelihoodCalculationSingleProcess::makeProcessNodes_(){" << endl;
#endif
  // add Independent Parameters
  const auto& paramProc = process_.getIndependentParameters();

  for (size_t i = 0; i < paramProc.size(); i++)
  {
    shareParameter_(ConfiguredParameter::create(getContext_(), paramProc[i]));
  }

  // Share dependencies with aliased parameters
  for (size_t i = 0; i < paramProc.size(); i++)
  {
    auto name = paramProc[i].getName();
    auto vs = process_.getAlias(name);

    auto dep = dynamic_cast<const ConfiguredParameter*>(&getParameters_()[i])->dependency(0);
    for (const auto& s:vs)
    {
      auto newacp = ConfiguredParameter::create(getContext_(), {dep}, process_.getParameter(s));
      shareParameter_(newacp);
      aliasParameters(name, s);
    }
  }

  const auto spcm = dynamic_cast<const SubstitutionProcessCollectionMember*>(&process_);


  auto& pl2(getParameters_());

  // rates node
  std::string suff = spcm ? ("_" + TextTools::toString(spcm->getRateDistributionNumber())) : "";

  auto rates = process_.getRateDistribution();
  if (rates && dynamic_cast<const ConstantRateDistribution*>(rates) == nullptr)
    processNodes_.ratesNode_ = ConfiguredParametrizable::createConfigured<DiscreteDistribution, ConfiguredDistribution>(getContext_(), *rates, pl2, suff);

  ///////
  // tree node
  suff = spcm ? ("_" + TextTools::toString(spcm->getTreeNumber())) : "";
  processNodes_.treeNode_ = ProcessTree::makeProcessTree(getContext_(), process_, pl2, suff);

  ///////////////////////////
  // rootFrequencies node

  auto root = process_.getRootFrequencySet();
  if (root)
  {
    suff = spcm ? ("_" + TextTools::toString(spcm->getRootFrequenciesNumber())) : "";
    processNodes_.rootFreqsNode_ = ConfiguredParametrizable::createConfigured<FrequencySet, ConfiguredFrequencySet>(getContext_(), *root, pl2, suff);
  }

  auto itE = processNodes_.treeNode_->allEdgesIterator();
  // get any modelNode from the map (only for StateMap)
  for (itE->start(); !itE->end(); itE->next())
  {
    if ((*(*itE))->getModel() != 0)
    {
      processNodes_.modelNode_ = (*(*itE))->getModel();
      break;
    }
  }
  if (!processNodes_.modelNode_)
    throw Exception("LikelihoodCalculationSingleProcess::makeProcessNodes_: null modelNode_");

#ifdef DEBUG
  cerr << "likelihoodcalculationsingleprocess::makeprocessnodes_()}" << endl;
#endif
}

void LikelihoodCalculationSingleProcess::makeProcessNodes_(CollectionNodes& collection, size_t nProc)
{
  auto& spcm = collection.getCollection().getSubstitutionProcess(nProc);

  // share process parameters with those of the collection
  const auto& paramProc = spcm.getIndependentParameters();

  for (size_t i = 0; i < paramProc.size(); i++)
  {
    auto name = paramProc[i].getName();
    if (!collection.hasParameter(name))
      throw Exception("LikelihoodCalculationSingleProcess::makeProcessNodes_ : CollectionNodes does not have parameter " + name);

    shareParameter_(collection.getSharedParameter(name));
  }

  // share dependencies

  for (size_t i = 0; i < paramProc.size(); i++)
  {
    auto name = paramProc[i].getName();
    auto vs = spcm.getAlias(name);

    auto dep = dynamic_cast<const ConfiguredParameter*>(&getParameters_()[i])->dependency(0);
    for (const auto& s:vs)
    {
      auto newacp = ConfiguredParameter::create(getContext_(), {dep}, spcm.getParameter(s));
      shareParameter_(newacp);
      aliasParameters(name, s);
    }
  }

  // rates node

  const DiscreteDistribution* rates = spcm.getRateDistribution();

  if (dynamic_cast<const ConstantRateDistribution*>(rates) == nullptr)
    processNodes_.ratesNode_ = collection.getRateDistribution(spcm.getRateDistributionNumber());

  ///////
  // tree node

  processNodes_.treeNode_ = ProcessTree::makeProcessTree(collection, nProc);

  ///////////////////////////
  // rootFrequencies node

  if (!spcm.isStationary())
    processNodes_.rootFreqsNode_ = collection.getFrequencies(spcm.getRootFrequenciesNumber());

  //////////////////////
  // get any modelNode from the map (only for StateMap)

  auto itE = processNodes_.treeNode_->allEdgesIterator();
  for (itE->start(); !itE->end(); itE->next())
  {
    if ((*(*itE))->getModel() != 0)
    {
      processNodes_.modelNode_ = (*(*itE))->getModel();
      break;
    }
  }
  if (!processNodes_.modelNode_)
    throw Exception("LikelihoodCalculationSingleProcess::makeProcessNodes_: null modelNode_");
}


void LikelihoodCalculationSingleProcess::setNumericalDerivateConfiguration(double delta, const NumericalDerivativeType& config)
{
  auto deltaNode = NumericMutable<double>::create(getContext_(), delta);

  if (processNodes_.ratesNode_)
  {
    processNodes_.ratesNode_->config.delta = deltaNode;
    processNodes_.ratesNode_->config.type = config;
  }

  /////////////////
  // model nodes

  vector<shared_ptr<ProcessEdge> > vpn = processNodes_.treeNode_->getAllEdges();

  for (auto& it: vpn)
  {
    auto mN = it->getModel();
    if (mN)
    {
      mN->config.delta = deltaNode;
      mN->config.type = config;
    }

    // auto pN=it->getProba();
    // {
    //   pN->config.delta = deltaNode;
    //   pN->config.type = config;
    // }
  }

  ///////////////////////////
  // rootFrequencies node

  if (processNodes_.rootFreqsNode_)
  {
    processNodes_.rootFreqsNode_->config.delta = deltaNode;
    processNodes_.rootFreqsNode_->config.type = config;
  }
}

void LikelihoodCalculationSingleProcess::setClockLike(double rate)
{
  Parameter pRate("BrLen_rate", rate, Parameter::R_PLUS_STAR);

  auto rateNode = ConfiguredParameter::create(getContext_(), pRate);

  auto rateRef = ValueFromConfiguredParameter::create(getContext_(), {rateNode});

  /////////////////
  // brlen nodes

  vector<shared_ptr<ProcessEdge> > vpn = processNodes_.treeNode_->getAllEdges();

  for (auto& it: vpn)
  {
    auto cp = it->getBrLen();

    if (cp)
    {
      auto mulref = CWiseMul<double, std::tuple<double, double> >::create (getContext_(), {cp->dependency(0), rateRef}, Dimension<double>());

      auto cp2 = ConfiguredParameter::resetDependencies(getContext_(), cp, {mulref});

      it->setBrLen(cp2);
    }
  }

  // Remove all BrLen parameters
  auto parNames = getParameters().getParameterNames();

  for (auto& name:parNames)
  {
    if (name.substr(0, 5) == "BrLen")
      deleteParameter_(name);
  }

  shareParameter_(rateNode);
}

RowLik LikelihoodCalculationSingleProcess::getSiteLikelihoodsForAClass(size_t nCat, bool shrunk)
{
  if (shrunk)
    return getSiteLikelihoodsTree_(nCat)->getRoot()->getTargetValue();
  else
    return expandVector(getSiteLikelihoodsTree_(nCat)->getRoot())->getTargetValue();
}

AllRatesSiteLikelihoods LikelihoodCalculationSingleProcess::getSiteLikelihoodsForAllClasses(bool shrunk)
{
  auto nbCat = vRateCatTrees_.size();
  auto allLk = std::make_shared<AllRatesSiteLikelihoods>(nbCat, shrunk ? getNumberOfDistinctSites() : getNumberOfSites());

  for (size_t nCat = 0; nCat < nbCat; nCat++)
  {
    allLk->row(Eigen::Index(nCat)) = getSiteLikelihoodsForAClass(nCat, shrunk);
  }
  return *allLk;
}

// std::shared_ptr<ConditionalLikelihoodTree> LikelihoodCalculationSingleProcess::getConditionalLikelihoodTree(size_t nCat)
// {
//   if (nCat>=vRateCatTrees_.size())
//     throw Exception("LikelihoodCalculationSingleProcess::getConditionalLikelihoodTree : Bad Class number " + TextTools::toString(nCat));

//   if (shrunkData_ && !likelihood_)
//     makeLikelihoodsAtRoot_();

//   if (vRateCatTrees_[nCat].clt==0)
//     makeLikelihoodsAtNode_(getTreeNode()->getRootIndex());

//   return vRateCatTrees_[nCat].clt;
// }


/****************************************
* Construction methods
****************************************/
void LikelihoodCalculationSingleProcess::makeRootFreqs_()
{
// Set root frequencies

  size_t nbState = getStateMap().getNumberOfModelStates();
  rFreqs_ = processNodes_.rootFreqsNode_ ? ConfiguredParametrizable::createRowVector<ConfiguredFrequencySet, FrequenciesFromFrequencySet, Eigen::RowVectorXd>(
    getContext_(), {processNodes_.rootFreqsNode_}, RowVectorDimension (Eigen::Index (nbState))) :
            ConfiguredParametrizable::createRowVector<ConfiguredModel, EquilibriumFrequenciesFromModel, Eigen::RowVectorXd>(
    getContext_(), {processNodes_.modelNode_}, RowVectorDimension (Eigen::Index (nbState)));
}


void LikelihoodCalculationSingleProcess::makeForwardLikelihoodTree_()
{
  // Build conditional likelihoods up to root recursively.
  if (!processNodes_.treeNode_->isRooted ())
  {
    throw Exception ("LikelihoodCalculationSingleProcess::makeForwardLikelihoodTree_ : PhyloTree must be rooted");
  }

  if (processNodes_.ratesNode_)
  {
    uint nbCat = (uint)processNodes_.ratesNode_->getTargetValue()->getNumberOfCategories();

    vRateCatTrees_.resize(nbCat);

    for (uint nCat = 0; nCat < nbCat; nCat++)
    {
      ValueRef<double> catRef = CategoryFromDiscreteDistribution::create(getContext_(), {processNodes_.ratesNode_}, nCat);

      auto treeCat = std::make_shared<ProcessTree>(*processNodes_.treeNode_, catRef);

      vRateCatTrees_[nCat].phyloTree = treeCat;

      auto flt = std::make_shared<ForwardLikelihoodTree>(getContext_(), treeCat, getStateMap());

      if (getShrunkData())
        flt->initialize(*getShrunkData());
      else
        flt->initialize(*psites_);
      vRateCatTrees_[nCat].flt = flt;
    }
  }
  else
  {
    vRateCatTrees_.resize(1);
    vRateCatTrees_[0].phyloTree = processNodes_.treeNode_;

    auto flt = std::make_shared<ForwardLikelihoodTree >(getContext_(), processNodes_.treeNode_, processNodes_.modelNode_->getTargetValue()->getStateMap());

    if (getShrunkData())
      flt->initialize(*getShrunkData());
    else
      flt->initialize(*psites_);
    vRateCatTrees_[0].flt = flt;
  }
}


void LikelihoodCalculationSingleProcess::makeLikelihoodsAtRoot_()
{
  if (vRateCatTrees_.size() == 0)
    makeForwardLikelihoodTree_();

  size_t nbDistSite = getNumberOfDistinctSites();

  // Set root frequencies
  if (rFreqs_ == 0)
    makeRootFreqs_();

  ValueRef<RowLik> sL;

  if (processNodes_.ratesNode_)
  {
    std::vector<std::shared_ptr<Node_DF> > vLikRoot;

    auto zero = NumericConstant<size_t>::create(getContext_(), size_t(0));

    for (auto& rateCat: vRateCatTrees_)
    {
      vLikRoot.push_back(LikelihoodFromRootConditionalAtRoot::create (
                           getContext_(), {rFreqs_, rateCat.flt->getForwardLikelihoodArrayAtRoot()},
                           RowVectorDimension (Eigen::Index (nbDistSite))));
    }

    auto catProb = ProbabilitiesFromDiscreteDistribution::create(getContext_(), {processNodes_.ratesNode_});

    for (size_t nCat = 0; nCat < vRateCatTrees_.size(); nCat++)
    {
      vLikRoot.push_back(ProbabilityFromDiscreteDistribution::create(getContext_(), {processNodes_.ratesNode_}, (uint)nCat));
    }

    sL = CWiseMean<RowLik, ReductionOf<RowLik>, ReductionOf<double> >::create(getContext_(), std::move(vLikRoot), RowVectorDimension (Eigen::Index(nbDistSite)));
  }
  else
  {
    sL = LikelihoodFromRootConditionalAtRoot::create (
      getContext_(), {rFreqs_, vRateCatTrees_[0].flt->getForwardLikelihoodArrayAtRoot()}, RowVectorDimension (Eigen::Index (nbDistSite)));
  }


  // likelihoods per distinct site
  setSiteLikelihoods(sL, true);

  // likelihoods per site
  setSiteLikelihoods(expandVector(patternedSiteLikelihoods_), false);

  // global likelihood
  ValueRef<DataLik> val;
  if (rootPatternLinks_)
    val = SumOfLogarithms<RowLik>::create (getContext_(), {sL, rootWeights_}, RowVectorDimension (Eigen::Index (nbDistSite)));
  else
    val = SumOfLogarithms<RowLik>::create (getContext_(), {sL}, RowVectorDimension (Eigen::Index (nbDistSite)));

  setLikelihoodNode(val);


#ifdef DEBUG
  using bpp::DotOptions;
  writeGraphToDot(
    "debug_lik.dot", {likelihood_.get()}, DotOptions::DetailedNodeInfo);
#endif
}


void LikelihoodCalculationSingleProcess::makeLikelihoodsAtNode_(uint speciesId)
{
  // Already built
  if (condLikelihoodTree_ && condLikelihoodTree_->hasNode(speciesId))
    return;

  if (vRateCatTrees_.size() == 0)
    makeForwardLikelihoodTree_();

  if (rFreqs_ == 0)
    makeRootFreqs_();

  const auto& stateMap = getStateMap();
  auto nbDistSite = Eigen::Index(getNumberOfDistinctSites());
  auto nbState = Eigen::Index(stateMap.getNumberOfModelStates());
  MatrixDimension likelihoodMatrixDim = conditionalLikelihoodDimension (nbState, nbDistSite);

  const auto& phylotree = process_.getParametrizablePhyloTree();

  ValueRef<RowLik> siteLikelihoodsNode;

  std::shared_ptr<ConditionalLikelihood> cond(0);

  std::vector<NodeRef> vCondRate;

  SiteLikelihoodsRef distinctSiteLikelihoodsNode;
  ConditionalLikelihoodRef conditionalLikelihoodsNode;

  std::vector<std::shared_ptr<Node_DF> > vRoot; // if several rates

  if (!condLikelihoodTree_)
    condLikelihoodTree_ = std::make_shared<ConditionalLikelihoodTree>(phylotree.getGraph());

  auto one = ConstantOne<Eigen::RowVectorXd>::create(getContext_(), RowVectorDimension (nbState));

  for (auto& rateCat: vRateCatTrees_)
  {
    if (!rateCat.blt)
      rateCat.blt = std::make_shared<BackwardLikelihoodTree>(getContext_(), rateCat.flt, rateCat.phyloTree, rFreqs_, stateMap, nbDistSite);

    if (!rateCat.clt)
      rateCat.clt = std::make_shared<ConditionalLikelihoodDAG>(rateCat.flt->getGraph());

    if (!rateCat.lt)
      rateCat.lt = std::make_shared<SiteLikelihoodsDAG>(rateCat.flt->getGraph());


    if (!rateCat.speciesLt)
      rateCat.speciesLt = std::make_shared<SiteLikelihoodsTree>(phylotree.getGraph());

    auto& dagIndexes = rateCat.flt->getDAGNodesIndexes(speciesId);

    std::vector<std::shared_ptr<Node_DF> > vCond;

    for (const auto& index : dagIndexes)
    {
      if (rateCat.clt->hasNode(index))
      {
        cond = rateCat.clt->getNode(index);
        if (dagIndexes.size() > 1) // for sum
          vCond.push_back(cond);
        continue;
      }

      auto condAbove = rateCat.blt->getBackwardLikelihoodArray(index);
      auto condBelow = rateCat.flt->getForwardLikelihoodArray(index);

      cond = BuildConditionalLikelihood::create (
        getContext_(), {condAbove, condBelow}, likelihoodMatrixDim);

      if (dagIndexes.size() > 1) // for sum
        vCond.push_back(cond);

      rateCat.clt->associateNode(cond, rateCat.flt->getNodeGraphid(rateCat.flt->getNode(index)));
      rateCat.clt->setNodeIndex(cond, index);

      // Site Likelihoods on this point
      auto lt = LikelihoodFromRootConditionalAtRoot::create (
        getContext_(), {one, cond}, RowVectorDimension (nbDistSite));

      rateCat.lt->associateNode(lt, rateCat.flt->getNodeGraphid(rateCat.flt->getNode(index)));
      rateCat.lt->setNodeIndex(lt, index);
    }

    /*
     * If several DAG nodes related with this species node, sum the
     * likelihoods of all (already multiplied by their probability).
     *
     */

    if (dagIndexes.size() > 1)
      cond = CWiseAdd<MatrixLik, ReductionOf<MatrixLik> >::create(getContext_(), std::move(vCond), likelihoodMatrixDim);

    // for Lik at Node
    auto siteLikelihoodsCat = LikelihoodFromRootConditionalAtRoot::create (
      getContext_(), {one, cond}, RowVectorDimension (nbDistSite));

    if (!rateCat.speciesLt->hasNode(speciesId))
    {
      rateCat.speciesLt->associateNode(siteLikelihoodsCat, phylotree.getNodeGraphid(phylotree.getNode(speciesId)));
      rateCat.speciesLt->setNodeIndex(siteLikelihoodsCat, speciesId);
    }

    if (!processNodes_.ratesNode_)
    {
      distinctSiteLikelihoodsNode = siteLikelihoodsCat;
      conditionalLikelihoodsNode = cond;
      break;
    }
    else
    {
      // For Conditional at Node
      vCondRate.push_back(cond);
      vRoot.push_back(siteLikelihoodsCat);
    }
  }

  if (processNodes_.ratesNode_)
  {
    auto catProb = ProbabilitiesFromDiscreteDistribution::create(getContext_(), {processNodes_.ratesNode_});
    auto catProbEf = Convert<RowLik, Eigen::RowVectorXd>::create(getContext_(), {catProb}, RowVectorDimension (Eigen::Index (nbDistSite)));
    vRoot.push_back(catProbEf);
    vCondRate.push_back(catProbEf);


    distinctSiteLikelihoodsNode = CWiseMean<RowLik, ReductionOf<RowLik>, RowLik>::create(getContext_(), std::move(vRoot), RowVectorDimension (nbDistSite));

    conditionalLikelihoodsNode = CWiseMean<MatrixLik, ReductionOf<MatrixLik>, RowLik>::create(getContext_(), std::move(vCondRate), MatrixDimension (nbState, nbDistSite));
  }

  condLikelihoodTree_->associateNode(conditionalLikelihoodsNode, phylotree.getNodeGraphid(phylotree.getNode(speciesId)));
  condLikelihoodTree_->setNodeIndex(conditionalLikelihoodsNode, speciesId);

  // We want -log(likelihood)
  // auto totalNegLogLikelihood =
  //   CWiseNegate<double>::create (getContext_(), {totalLogLikelihood}, Dimension<double> ());
  //  return totalNegLogLikelihood;
}

void LikelihoodCalculationSingleProcess::makeLikelihoodsAtDAGNode_(uint nodeId)
{
  if (vRateCatTrees_.size() == 0)
    makeForwardLikelihoodTree_();

  if (rFreqs_ == 0)
    makeRootFreqs_();

  const auto& stateMap = getStateMap();
  auto nbDistSite = Eigen::Index(getNumberOfDistinctSites());
  auto nbState = Eigen::Index(stateMap.getNumberOfModelStates());
  MatrixDimension likelihoodMatrixDim = conditionalLikelihoodDimension (nbState, nbDistSite);

  ValueRef<RowLik> siteLikelihoodsNode;

  auto one = ConstantOne<Eigen::RowVectorXd>::create(getContext_(), RowVectorDimension (Eigen::Index (nbState)));

  for (auto& rateCat: vRateCatTrees_)
  {
    if (!rateCat.clt)
      rateCat.clt = std::make_shared<ConditionalLikelihoodDAG>(rateCat.flt->getGraph());

    if (rateCat.clt->hasNode(nodeId)) // already computed
      continue;

    if (!rateCat.blt)
      rateCat.blt = std::make_shared<BackwardLikelihoodTree>(getContext_(), rateCat.flt, rateCat.phyloTree, rFreqs_, stateMap, nbDistSite);

    if (!rateCat.lt)
      rateCat.lt = std::make_shared<SiteLikelihoodsDAG>(rateCat.flt->getGraph());

    // Conditional Likelihoods on this node

    auto condAbove = rateCat.blt->getBackwardLikelihoodArray(nodeId);
    auto condBelow = rateCat.flt->getForwardLikelihoodArray(nodeId);

    auto cond = BuildConditionalLikelihood::create (
      getContext_(), {condAbove, condBelow}, likelihoodMatrixDim);

    rateCat.clt->associateNode(cond, rateCat.flt->getNodeGraphid(rateCat.flt->getNode(nodeId)));
    rateCat.clt->setNodeIndex(cond, nodeId);

    // Site Likelihoods on this node
    auto lt = LikelihoodFromRootConditionalAtRoot::create (
      getContext_(), {one, cond}, RowVectorDimension (Eigen::Index (nbDistSite)));

    rateCat.lt->associateNode(lt, rateCat.flt->getNodeGraphid(rateCat.flt->getNode(nodeId)));
    rateCat.lt->setNodeIndex(lt, nodeId);
  }
}


std::shared_ptr<SiteLikelihoodsTree> LikelihoodCalculationSingleProcess::getSiteLikelihoodsTree_(size_t nCat)
{
  if (nCat >= vRateCatTrees_.size())
    throw Exception("LikelihoodCalculationSingleProcess::getSiteLikelihoodsTree : Bad Class number " + TextTools::toString(nCat));

  if (!shrunkData_)
    throw Exception("LikelihoodCalculationSingleProcess::getSiteLikelihoodsTree : data not set.");

  if (!getLikelihoodNode_())
    makeLikelihoodsAtRoot_();

  if (vRateCatTrees_[nCat].speciesLt == 0)
    makeLikelihoodsAtNode_(getTreeNode(nCat)->getRoot()->getSpeciesIndex());

  return vRateCatTrees_[nCat].speciesLt;
}


ConditionalLikelihoodRef LikelihoodCalculationSingleProcess::getForwardLikelihoodsAtNodeForClass(uint nodeId, size_t nCat)
{
  // compute forward likelihoods for all nodes (not the quickest, but
  // in pratice they are all needed)

  if (!getLikelihoodNode_())
    makeLikelihoods();

  if (nCat >= vRateCatTrees_.size())
    throw Exception("LikelihoodCalculationSingleProcess::getForwardLikelihoodsAtNodeForClass : bad class number " + TextTools::toString(nCat));

  return vRateCatTrees_[nCat].flt->getNode(nodeId);
}

ConditionalLikelihoodRef LikelihoodCalculationSingleProcess::getConditionalLikelihoodsAtNodeForClass(uint nodeId, size_t nCat)
{
  // compute likelihoods for all edges with similar species index
  // (not the quickest, but in pratice they are all needed)

  makeLikelihoodsAtDAGNode_(nodeId);

  if (nCat >= vRateCatTrees_.size())
    throw Exception("LikelihoodCalculationSingleProcess::getConditionalLikelihoodsAtNodeForClass : bad class number " + TextTools::toString(nCat));

  return vRateCatTrees_[nCat].clt->getNode(nodeId);
}

SiteLikelihoodsRef LikelihoodCalculationSingleProcess::getLikelihoodsAtNodeForClass(uint nodeId, size_t nCat)
{
  // compute likelihoods for all edges with similar species index
  // (not the quickest, but in pratice they are all needed)

  makeLikelihoodsAtDAGNode_(nodeId);

  if (nCat >= vRateCatTrees_.size())
    throw Exception("LikelihoodCalculationSingleProcess::getConditionalLikelihoodsAtNodeForClass : bad class number " + TextTools::toString(nCat));

  return vRateCatTrees_[nCat].lt->getNode(nodeId);
}


ConditionalLikelihoodRef LikelihoodCalculationSingleProcess::getBackwardLikelihoodsAtEdgeForClass(uint edgeId, size_t nCat)
{
  // compute likelihoods for all edges with similar species index
  // (not the quickest, but in pratice they are all needed)

  auto spId = getTreeNode(nCat)->getEdge(edgeId)->getSpeciesIndex();

  if (!(condLikelihoodTree_ && condLikelihoodTree_->hasNode(spId)))
    makeLikelihoodsAtNode_(spId);

  if (nCat >= vRateCatTrees_.size())
    throw Exception("LikelihoodCalculationSingleProcess::getForwardLikelihoodsAtNodeForClass : bad class number " + TextTools::toString(nCat));

  return vRateCatTrees_[nCat].blt->getEdge(edgeId);
}

ConditionalLikelihoodRef LikelihoodCalculationSingleProcess::getBackwardLikelihoodsAtNodeForClass(uint nodeId, size_t nCat)
{
  // compute likelihoods for all nodes with similar species index
  // (not the quickest, but in pratice they are all needed)

  makeLikelihoodsAtDAGNode_(nodeId);

  if (nCat >= vRateCatTrees_.size())
    throw Exception("LikelihoodCalculationSingleProcess::getForwardLikelihoodsAtNodeForClass : bad class number " + TextTools::toString(nCat));

  return vRateCatTrees_[nCat].blt->getNode(nodeId);
}


const DAGindexes& LikelihoodCalculationSingleProcess::getNodesIds(uint speciesId) const
{
  if (vRateCatTrees_.size() == 0)
    throw Exception("LikelihoodCalculationSingleProcess::getNodeIds. ForwardLikelihoodTree not computed.");

  return vRateCatTrees_[0].flt->getDAGNodesIndexes(speciesId);
}

const DAGindexes& LikelihoodCalculationSingleProcess::getEdgesIds(uint speciesId, size_t nCat) const
{
  if (nCat >= vRateCatTrees_.size())
    throw Exception("LikelihoodCalculationSingleProcess::getEdgesIds : bad class number " + TextTools::toString(nCat));

  return vRateCatTrees_[nCat].flt->getDAGEdgesIndexes(speciesId);
}

std::shared_ptr<ForwardLikelihoodTree> LikelihoodCalculationSingleProcess::getForwardLikelihoodTree(size_t nCat)
{
  if (nCat >= vRateCatTrees_.size())
    throw Exception("LikelihoodCalculationSingleProcess::getForwardTree : bad class number " + TextTools::toString(nCat));

  return vRateCatTrees_[nCat].flt;
}
