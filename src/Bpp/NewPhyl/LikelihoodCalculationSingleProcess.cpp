//
// File: LikelihoodCalculationSingleProcess.cpp
// Authors: François Gindraud, Laurent Guéguen (2018)
// Created: jeudi 28 février 2019, à 07h 22
//

#include "Bpp/NewPhyl/LikelihoodCalculationSingleProcess.h"
#include "Bpp/NewPhyl/ForwardLikelihoodTree.h"
#include "Bpp/NewPhyl/BackwardLikelihoodTree.h"

#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/MixedTransitionModel.h>
#include "Bpp/Phyl/NewLikelihood/SubstitutionProcessCollectionMember.h"

#include <unordered_map>
#include <list>
#include <numeric>

using namespace std;
using namespace bpp;
using namespace bpp::dataflow;

LikelihoodCalculationSingleProcess::LikelihoodCalculationSingleProcess(Context& context,
                                                                       const AlignedValuesContainer & sites,
                                                                       const SubstitutionProcess& process):
  AbstractParametrizable(""),
  context_(context), process_(process), psites_(&sites),
  rootPatternLinks_(), rootWeights_(), shrunkData_(),
  processNodes_(), rFreqs_(),
  likelihood_(), siteLikelihoods_(), patternedSiteLikelihoods_(),
  vRateCatTrees_()
{
  setPatterns_();
  makeProcessNodes_();
}

LikelihoodCalculationSingleProcess::LikelihoodCalculationSingleProcess(Context & context,
                                                                       const SubstitutionProcess& process):
  AbstractParametrizable(""),
  context_(context), process_(process), psites_(),
  rootPatternLinks_(), rootWeights_(), shrunkData_(),
  processNodes_(), rFreqs_(),
  likelihood_(), siteLikelihoods_(), patternedSiteLikelihoods_(),
  vRateCatTrees_()
{
  makeProcessNodes_();
}

LikelihoodCalculationSingleProcess::LikelihoodCalculationSingleProcess(Context & context,
                                                                       const AlignedValuesContainer & sites,
                                                                       const SubstitutionProcess& process,
                                                                       ParameterList& paramList):
  AbstractParametrizable(""),
  context_(context), process_(process), psites_(&sites),
  rootPatternLinks_(), rootWeights_(), shrunkData_(),
  processNodes_(), rFreqs_(),
  likelihood_(), siteLikelihoods_(), patternedSiteLikelihoods_(),
  vRateCatTrees_()
{
  setPatterns_();
  makeProcessNodes(paramList);
}


LikelihoodCalculationSingleProcess::LikelihoodCalculationSingleProcess(const LikelihoodCalculationSingleProcess& lik) :
  AbstractParametrizable(lik),
  context_(*std::shared_ptr<Context>().get()), process_(lik.process_), psites_(lik.psites_),
  rootPatternLinks_(lik.rootPatternLinks_), rootWeights_(), shrunkData_(lik.shrunkData_),
  processNodes_(), rFreqs_(),
  likelihood_(), siteLikelihoods_(), patternedSiteLikelihoods_(),
  vRateCatTrees_()
{
  setPatterns_();
  makeProcessNodes_();
}

void LikelihoodCalculationSingleProcess::setPatterns_()
{
  SitePatterns patterns(psites_);
  shrunkData_       = patterns.getSites();
  rootPatternLinks_ = NumericConstant<PatternType>::create(context_, patterns.getIndices());
  size_t nbSites    = shrunkData_->getNumberOfSites();
  Eigen::RowVectorXi weights(nbSites);
  for (std::size_t i=0;i<nbSites;i++)
    weights(Eigen::Index(i))=patterns.getWeights()[i];
  rootWeights_ = SiteWeights::create(context_, std::move(weights));
}

void LikelihoodCalculationSingleProcess::makeProcessNodes_()
{
  ParameterList paramList;
  
  // add Independent Parameters
  const auto& paramProc=process_.getIndependentParameters();
  
  for (size_t i=0;i<paramProc.size();i++)
    paramList.shareParameter(ConfiguredParameter::create(context_, paramProc[i]));
  
  // Share dependencies with aliased parameters

  for (size_t i=0;i<paramProc.size();i++)
  {
    auto vs=process_.getAlias(paramProc[i].getName());
    auto dep=dynamic_cast<const ConfiguredParameter*>(&paramList[i])->dependency(0);
    for (const auto& s:vs)
    {
      auto newacp = ConfiguredParameter::create(context_, {dep}, process_.getParameter(s));
      paramList.shareParameter(newacp);
    }
  }
  
  makeProcessNodes(paramList);
}

void LikelihoodCalculationSingleProcess::makeProcessNodes(ParameterList& paramList)
{
  const auto spcm=dynamic_cast<const SubstitutionProcessCollectionMember*>(&process_);

  // share process_ parameters with those of the paramList
  const auto& paramProc=process_.getParameters();
  
  for (size_t i=0;i<paramProc.size();i++)
  {
    auto name=paramProc[i].getName();
    if (!paramList.hasParameter(name))
      throw Exception("LikelihoodCalculationSingleProcess::makeProcessNodes : paramList does not have parameter " + name);
    auto* confPar=dynamic_cast<ConfiguredParameter*>(&paramList.getParameter(name));
    if (!confPar)
      throw Exception("LikelihoodCalculationSingleProcess::makeProcessNodes : parameter " + name + "is not a ConfiguredParameter.");
      
    shareParameter_(paramList.getSharedParameter(name));
  }
  
  // // Share dependencies with aliased parameters

  // for (size_t i=0;i<paramProc.size();i++)
  // {
  //   auto vs=process_.getAlias(paramProc[i].getName());
  //   auto dep=dynamic_cast<const ConfiguredParameter*>(&paramList.getParameters()[i])->dependency(0);
  //   for (const auto& s:vs)
  //   {
  //     auto newacp = ConfiguredParameter::create(context_, {dep}, process_.getParameter(s));
  //     paramList.shareParameter_(newacp);
  //   }
  // }

  // rates node
  std::string suff=spcm?("_"+TextTools::toString(spcm->getRateDistributionNumber())):"";
  
  auto rates = process_.getRateDistribution(); 
  if (rates && dynamic_cast<const ConstantRateDistribution*>(rates)==nullptr)
    processNodes_.ratesNode_ = ConfiguredParametrizable::createConfigured<DiscreteDistribution, ConfiguredDistribution>(context_, *rates, paramList, suff);


  /////////////////
  // model nodes
  
  auto vnMod=process_.getModelNumbers();
        
  ModelMap modelmap;

  for (auto nMod:vnMod)
  {    
    auto modelNode = ConfiguredParametrizable::createConfigured<TransitionModel, ConfiguredModel>(context_, *process_.getModel(nMod), paramList, "_"+ TextTools::toString(nMod));

    const auto mtm=dynamic_cast<const MixedTransitionModel*>(process_.getModel(nMod));

    std::vector<size_t> nsubmod(mtm?mtm->getNumberOfModels():0);
    std::iota (std::begin(nsubmod), std::end(nsubmod), 0);

    // assign model to branche id
    std::vector<uint> vId=process_.getNodesWithModel(nMod);
    ModelAssign mass(modelNode, nsubmod);
    
    for (auto id:vId)
      modelmap.emplace(id, mass);
  }

  ///////
  // tree node
  suff=spcm?("_"+TextTools::toString(spcm->getTreeNumber())):"";
  const ParametrizablePhyloTree& parTree = process_.getParametrizablePhyloTree();
  processNodes_.treeNode_ = makeProcessTree(context_, parTree, modelmap, paramList, suff);

  ///////////////////////////
  // rootFrequencies node

  auto root = process_.getRootFrequenciesSet();
  if (root)
  {
    suff=spcm?("_"+TextTools::toString(spcm->getRootFrequenciesNumber())):"";
    processNodes_.rootFreqsNode_ = ConfiguredParametrizable::createConfigured<FrequenciesSet, ConfiguredFrequenciesSet>(context_, *root, paramList, suff);
  }

  // get any modelNode from the map (only for StateMap)
  processNodes_.modelNode_ = modelmap.begin()->second.model_;
}


void LikelihoodCalculationSingleProcess::setNumericalDerivateConfiguration(double delta, const NumericalDerivativeType& config)
{
  auto deltaNode = NumericMutable<double>::create(context_, delta);

  if (processNodes_.ratesNode_)
  {
    processNodes_.ratesNode_->config.delta = deltaNode;
    processNodes_.ratesNode_->config.type = config;
  }

  /////////////////
  // model nodes

  vector<shared_ptr<ProcessEdge> > vpn=processNodes_.treeNode_->getAllEdges();

  for (auto& it: vpn)
  {
    auto mN=it->getModel();
    mN->config.delta = deltaNode;
    mN->config.type = config;
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

  auto rateNode = ConfiguredParameter::create(context_, pRate);

  auto rateRef= ValueFromConfiguredParameter::create(context_, {rateNode});

  /////////////////
  // brlen nodes

  vector<shared_ptr<ProcessEdge> > vpn=processNodes_.treeNode_->getAllEdges();

  for (auto& it: vpn)
  {
    auto cp=it->getBrLen();

    auto mulref = CWiseMul<double, std::tuple<double, double>>::create (context_, {cp->dependency(0), rateRef}, Dimension<double>());

    auto cp2=ConfiguredParameter::resetDependencies(context_, cp, {mulref});
    
    it->setBrLen(cp2);
  }

  // Remove all BrLen parameters
  auto parNames=getParameters().getParameterNames();
  
  for (auto& name:parNames)
  {
    if (name.substr(0,5)=="BrLen")
      deleteParameter_(name);
  }

  shareParameter_(rateNode);
}

const Eigen::RowVectorXd& LikelihoodCalculationSingleProcess::getSiteLikelihoodsForAClass(size_t nCat)
{
  return getSiteLikelihoodsTree(nCat)->getRoot()->getTargetValue();
}

const AllRatesSiteLikelihoods& LikelihoodCalculationSingleProcess::getSiteLikelihoodsForAllClasses()
{
  auto nbCat=vRateCatTrees_.size();
  auto allLk=std::make_shared<AllRatesSiteLikelihoods>(nbCat,getNumberOfDistinctSites());
  for (size_t nCat=0;nCat<nbCat;nCat++)
      allLk->row(nCat)=getSiteLikelihoodsForAClass(nCat);
  
  return *allLk;
}

std::shared_ptr<SiteLikelihoodsTree> LikelihoodCalculationSingleProcess::getSiteLikelihoodsTree(size_t nCat)
{
  if (nCat>=vRateCatTrees_.size())
    throw Exception("LikelihoodCalculationSingleProcess::getSiteLikelihoodsTree : Bad Class number " + TextTools::toString(nCat));

  if (!shrunkData_)
    throw Exception("LikelihoodCalculationSingleProcess::getSiteLikelihoodsTree : data not set.");
        
  if (!likelihood_)
    makeLikelihoodAtRoot_();

  if (vRateCatTrees_[nCat].lt==0)
    makeLikelihoodsAtNode_(getTreeNode()->getRootIndex());

  return vRateCatTrees_[nCat].lt;       
}

std::shared_ptr<ConditionalLikelihoodTree> LikelihoodCalculationSingleProcess::getConditionalLikelihoodTree(size_t nCat)
{
  if (nCat>=vRateCatTrees_.size())
    throw Exception("LikelihoodCalculationSingleProcess::getConditionalLikelihoodTree : Bad Class number " + TextTools::toString(nCat));
  
  if (shrunkData_ && !likelihood_)
    makeLikelihoodAtRoot_();
  
  if (vRateCatTrees_[nCat].clt==0)
    makeLikelihoodsAtNode_(getTreeNode()->getRootIndex());
  
  return vRateCatTrees_[nCat].clt;       
}


/****************************************
 * Construction methods
 ****************************************/

void LikelihoodCalculationSingleProcess::makeRootFreqs_()
{
// Set root frequencies 

  size_t nbState = getStateMap().getNumberOfModelStates();
  rFreqs_ = processNodes_.rootFreqsNode_?ConfiguredParametrizable::createVector<ConfiguredFrequenciesSet, FrequenciesFromFrequenciesSet> (
    context_, {processNodes_.rootFreqsNode_}, rowVectorDimension (Eigen::Index (nbState))):
    ConfiguredParametrizable::createVector<ConfiguredModel, EquilibriumFrequenciesFromModel> (
      context_, {processNodes_.modelNode_}, rowVectorDimension (Eigen::Index (nbState)));
}


void LikelihoodCalculationSingleProcess::makeForwardLikelihoodTree_()
{
  // Build conditional likelihoods up to root recursively.
  if (!processNodes_.treeNode_->isRooted ()) {
    throw Exception ("PhyloTree must be rooted");
  }
  
  if (processNodes_.ratesNode_)
  {
    auto brlenmap = createBrLenMap(context_, *processNodes_.treeNode_);
      
    uint nbCat=(uint)processNodes_.ratesNode_->getTargetValue()->getNumberOfCategories();

    vRateCatTrees_.resize(nbCat);

    for (uint nCat=0; nCat<nbCat; nCat++)
    {
      auto catRef = CategoryFromDiscreteDistribution::create(context_, {processNodes_.ratesNode_}, nCat);
  
      auto rateBrLen = multiplyBrLenMap(context_, *processNodes_.treeNode_, catRef);
      auto treeCat = std::shared_ptr<ProcessTree>(new ProcessTree(*processNodes_.treeNode_, std::move(rateBrLen)));
      vRateCatTrees_[nCat].phyloTree=treeCat;      

      auto flt=std::make_shared<ForwardLikelihoodTree>(context_, treeCat, getStateMap());
      flt->initialize(*getShrunkData());
      vRateCatTrees_[nCat].flt=flt;
    }
  }
  else
  {
    vRateCatTrees_.resize(1);
    vRateCatTrees_[0].phyloTree=processNodes_.treeNode_;
    
    auto flt=std::make_shared<ForwardLikelihoodTree >(context_, processNodes_.treeNode_, processNodes_.modelNode_->getTargetValue()->getStateMap());
    flt->initialize(*getShrunkData());
    vRateCatTrees_[0].flt=flt;
  }
}


void LikelihoodCalculationSingleProcess::makeLikelihoodAtRoot_()
{
  if (vRateCatTrees_.size()==0)
    makeForwardLikelihoodTree_();

  std::size_t nbSite = getShrunkData()->getNumberOfSites();    

  // Set root frequencies
  if (rFreqs_==0)
    makeRootFreqs_();
  
  if (processNodes_.ratesNode_)
  {
    std::vector<std::shared_ptr<Node>> vLogRoot;
      
    for (auto& rateCat: vRateCatTrees_)
    {
      vLogRoot.push_back(LikelihoodFromRootConditional::create (
                           context_, {rFreqs_, rateCat.flt->getForwardLikelihoodArrayAtRoot()}, rowVectorDimension (Eigen::Index (nbSite))));
    }

    auto catProb = ProbabilitiesFromDiscreteDistribution::create(context_, {processNodes_.ratesNode_});
    
    vLogRoot.push_back(catProb);
      
    siteLikelihoods_ = CWiseMean<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>, Eigen::RowVectorXd>::create(context_, std::move(vLogRoot), rowVectorDimension (Eigen::Index(nbSite)));
  }
  else
  {
    siteLikelihoods_ = LikelihoodFromRootConditional::create (
      context_, {rFreqs_, vRateCatTrees_[0].flt->getForwardLikelihoodArrayAtRoot()}, rowVectorDimension (Eigen::Index (nbSite)));
  }

  patternedSiteLikelihoods_= CWisePattern<Eigen::RowVectorXd>::create(context_,{siteLikelihoods_,rootPatternLinks_}, rowVectorDimension (Eigen::Index (getData()->getNumberOfSites())));

  auto totalLogLikelihood =
    SumOfLogarithms<Eigen::RowVectorXd>::create (context_, {siteLikelihoods_, rootWeights_}, rowVectorDimension (Eigen::Index (nbSite)));

// We want -log(likelihood)
  likelihood_ = CWiseNegate<double>::create (context_, {totalLogLikelihood}, Dimension<double> ());
}


void LikelihoodCalculationSingleProcess::makeLikelihoodsAtNode_(uint nodeId)
{
  if (vRateCatTrees_.size()==0)
    makeForwardLikelihoodTree_();
  
  if (rFreqs_==0)
    makeRootFreqs_();

  auto treeNode = processNodes_.treeNode_;
  
  const auto& stateMap = getStateMap();
  size_t nbSite = getShrunkData()->getNumberOfSites();
  size_t nbState = stateMap.getNumberOfModelStates();
  MatrixDimension likelihoodMatrixDim = conditionalLikelihoodDimension (nbState, nbSite);
  
  auto one=ConstantOne<Eigen::RowVectorXd>::create(context_, rowVectorDimension (Eigen::Index (nbState)));

  ValueRef<Eigen::RowVectorXd> siteLikelihoodsNode;

  if (processNodes_.ratesNode_)
  {
    std::vector<std::shared_ptr<Node>> vLogRoot;
      
    for (auto& rateCat: vRateCatTrees_)
    {
      if (!rateCat.blt)
        rateCat.blt=std::make_shared<BackwardLikelihoodTree>(context_, rateCat.flt, rateCat.phyloTree, rFreqs_, stateMap, nbSite);

      if (!rateCat.clt)
        rateCat.clt=std::make_shared<ConditionalLikelihoodTree>(treeNode->getGraph());

      if (!rateCat.lt)
        rateCat.lt=std::make_shared<SiteLikelihoodsTree>(treeNode->getGraph());

      auto condAbove = rateCat.blt->getBackwardLikelihoodArray(nodeId);
      auto condBelow = rateCat.flt->getForwardLikelihoodArray(nodeId);

      auto cond = BuildConditionalLikelihood::create (
        context_, {condAbove, condBelow}, likelihoodMatrixDim);
      
      rateCat.clt->associateNode(cond, treeNode->getNodeGraphid(treeNode->getNode(nodeId)));
      rateCat.clt->setNodeIndex(cond, nodeId);

      auto siteLikelihoodsCat = LikelihoodFromRootConditional::create (
        context_, {one, cond}, rowVectorDimension (Eigen::Index (nbSite)));

      auto patLikCat = CWisePattern<Eigen::RowVectorXd>::create(context_,{siteLikelihoodsCat, rootPatternLinks_}, rowVectorDimension (Eigen::Index(getData()->getNumberOfSites())));

      rateCat.lt->associateNode(patLikCat, treeNode->getNodeGraphid(treeNode->getNode(nodeId)));
      rateCat.lt->setNodeIndex(patLikCat, nodeId);

      vLogRoot.push_back(std::move(siteLikelihoodsCat));
    }

    auto catProb = ProbabilitiesFromDiscreteDistribution::create(context_, {processNodes_.ratesNode_});
    vLogRoot.push_back(catProb);
      
    siteLikelihoodsNode = CWiseMean<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>, Eigen::RowVectorXd>::create(context_, std::move(vLogRoot), rowVectorDimension (Eigen::Index(nbSite)));
  }
  else
  {
    auto& rateCat = vRateCatTrees_[0];
    
    if (!rateCat.blt)
      rateCat.blt=std::make_shared<BackwardLikelihoodTree>(context_, rateCat.flt, rateCat.phyloTree, rFreqs_, stateMap, nbSite);

    if (!rateCat.clt)
      rateCat.clt=std::make_shared<ConditionalLikelihoodTree>(treeNode->getGraph());
    if (!rateCat.lt)
      rateCat.lt=std::make_shared<SiteLikelihoodsTree>(treeNode->getGraph());
    
    auto condAbove = rateCat.blt->getBackwardLikelihoodArray(nodeId);
      
    auto condBelow = rateCat.flt->getForwardLikelihoodArray(nodeId);

    auto cond = BuildConditionalLikelihood::create (
      context_, {condAbove, condBelow}, likelihoodMatrixDim);
    
    rateCat.clt->associateNode(cond, treeNode->getNodeGraphid(treeNode->getNode(nodeId)));
    rateCat.clt->setNodeIndex(cond, nodeId);

    siteLikelihoodsNode = LikelihoodFromRootConditional::create (
      context_, {one, cond}, rowVectorDimension (Eigen::Index (nbSite)));

    auto patLikCat = CWisePattern<Eigen::RowVectorXd>::create(context_,{siteLikelihoodsNode, rootPatternLinks_}, rowVectorDimension (Eigen::Index(getData()->getNumberOfSites())));

    rateCat.lt->associateNode(patLikCat, treeNode->getNodeGraphid(treeNode->getNode(nodeId)));
    rateCat.lt->setNodeIndex(patLikCat, nodeId);
  }
//   auto totalLogLikelihood =
//     SumOfLogarithms<Eigen::RowVectorXd>::create (context_, {siteLikelihoodsNode}, rowVectorDimension (Eigen::Index (nbSite)));

// // We want -log(likelihood)
//   auto totalNegLogLikelihood =
//     CWiseNegate<double>::create (context_, {totalLogLikelihood}, Dimension<double> ());

//   return totalNegLogLikelihood;
}

