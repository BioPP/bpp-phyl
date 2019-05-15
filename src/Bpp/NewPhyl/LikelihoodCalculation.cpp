//
// File: LikelihoodCalculation.cpp
// Authors: François Gindraud, Laurent Guéguen (2018)
// Created: jeudi 28 février 2019, à 07h 22
//

#include "Bpp/NewPhyl/LikelihoodCalculation.h"
#include "Bpp/NewPhyl/ForwardLikelihoodTree.h"
#include "Bpp/NewPhyl/BackwardLikelihoodTree.h"

#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>

#include <unordered_map>
#include <list>

using namespace std;
using namespace bpp;
using namespace bpp::dataflow;

LikelihoodCalculation::LikelihoodCalculation(Context& context,
                                             const AlignedValuesContainer & sites,
                                             const SubstitutionProcess& process):
  AbstractParametrizable(""),
  context_(context), process_(process), psites_(&sites),
  rootPatternLinks_(), rootWeights_(), shrunkData_(),
  treeNode_(), modelNode_(), rootFreqsNode_(),  ratesNode_(), rFreqs_(),
  likelihood_(), siteLikelihoods_(), vRateCatTrees_()
{
  setPatterns_();
  makeProcessNodes_();
}

LikelihoodCalculation::LikelihoodCalculation(const LikelihoodCalculation& lik) :
  AbstractParametrizable(lik),
  context_(*std::shared_ptr<Context>().get()), process_(lik.process_), psites_(lik.psites_),
  rootPatternLinks_(lik.rootPatternLinks_), rootWeights_(), shrunkData_(lik.shrunkData_),
  treeNode_(), modelNode_(), rootFreqsNode_(), ratesNode_(), rFreqs_(),
  likelihood_(), siteLikelihoods_(), vRateCatTrees_()
{
  setPatterns_();
  makeProcessNodes_();
}

LikelihoodCalculation::LikelihoodCalculation(Context & context,
                                             const SubstitutionProcess& process):
  AbstractParametrizable(""),
  context_(context), process_(process), psites_(),
  rootPatternLinks_(), rootWeights_(), shrunkData_(),
  treeNode_(), modelNode_(), rootFreqsNode_(), ratesNode_(), rFreqs_(),
  likelihood_(), siteLikelihoods_(), vRateCatTrees_()
{
  makeProcessNodes_();
}


void LikelihoodCalculation::setPatterns_()
{
  SitePatterns patterns(psites_);
  shrunkData_       = patterns.getSites();
  rootPatternLinks_ = patterns.getIndices();
  size_t nbSites    = shrunkData_->getNumberOfSites();
  Eigen::RowVectorXi weights(nbSites);
  for (std::size_t i=0;i<nbSites;i++)
    weights(Eigen::Index(i))=patterns.getWeights()[i];
  rootWeights_ = SiteWeights::create(context_, std::move(weights));
}

void LikelihoodCalculation::makeProcessNodes_()
{
  // add Independent Parameters
  const auto& paramProc=process_.getIndependentParameters();
  for (size_t i=0;i<paramProc.size();i++)
    shareParameter_(ConfiguredParameter::create(context_, paramProc[i]));

  // Share dependencies with aliased parameters

  for (size_t i=0;i<paramProc.size();i++)
  {
    auto vs=process_.getAlias(paramProc[i].getName());
    auto dep=dynamic_cast<const ConfiguredParameter*>(&getParameters()[i])->dependency(0);
    for (const auto& s:vs)
    {
      auto newacp = ConfiguredParameter::create(context_, {dep}, process_.getParameter(s));
      shareParameter_(newacp);
    }
  }
  
  // rates nodes
  auto rates = process_.getRateDistribution();
  
  if (rates && dynamic_cast<const ConstantRateDistribution*>(rates)==nullptr)
  {
    std::unique_ptr<DiscreteDistribution> newRates(rates->clone());
    const auto& rateParams=newRates->getIndependentParameters();
    std::vector<NodeRef> deps;
    for (size_t i=0;i<rateParams.size();i++)
    {
      std::string name=rateParams[i].getName()+(hasParameter(rateParams[i].getName())?"":"_1");
      deps.push_back(ConfiguredParameter::create(context_, {dynamic_cast<ConfiguredParameter*>(getSharedParameter(name).get())->dependency(0)}, rateParams[i]));
    }
    
    ratesNode_ = ConfiguredParametrizable::createConfigured<DiscreteDistribution, ConfiguredDistribution>(
      context_,
      std::move(deps),
      std::move(newRates));
  }

  ///////
  // tree node
  const ParametrizablePhyloTree& parTree = process_.getParametrizablePhyloTree();
  std::vector<std::shared_ptr<PhyloBranchParam> > vB=parTree.getAllEdges();

  BrLenMap mapBr;

  for (auto& branch:vB)
  {
    const auto& bp=branch->getParameters()[0];
    
    std::string name=bp.getName()+(hasParameter(bp.getName())?"":"_1");
    mapBr.emplace(parTree.getEdgeIndex(branch),
                  ConfiguredParameter::create(context_, {dynamic_cast<ConfiguredParameter*>(getSharedParameter(name).get())->dependency(0)}, bp));
  }
  
  treeNode_ = std::shared_ptr<PhyloTree_BrRef>(new PhyloTree_BrRef(parTree, mapBr));

  /////////////////
  // model nodes
  
  std::vector<std::shared_ptr<ConfiguredModel>> vModelNodes;

  auto vnMod=process_.getModelNumbers();
        
  ModelMap modelmap;

  for (auto nMod:vnMod)
  {
    std::unique_ptr<TransitionModel> newModel(process_.getModel(nMod)->clone());
    
    const auto& modelParams=process_.getModel(nMod)->getParameters();

    std::vector<NodeRef> depModel;
    for (size_t np=0;np<modelParams.size();np++)
    {
      std::string name = modelParams[np].getName();
      
      if (hasParameter(name+"_"+ TextTools::toString(nMod)))
        depModel.push_back(ConfiguredParameter::create(context_, {dynamic_cast<ConfiguredParameter*>(getSharedParameter(name+"_"+ TextTools::toString(nMod)).get())->dependency(0)}, modelParams[np]));
      else
      {
        if (vnMod.size()==1)
          depModel.push_back(ConfiguredParameter::create(context_, {dynamic_cast<ConfiguredParameter*>(getSharedParameter(name).get())->dependency(0)}, modelParams[np]));
        else
          throw ParameterNotFoundException("LikelihoodCalculation::makeProcessNodes_", modelParams[np].getName());
      }
    }
    
    auto modelNode = ConfiguredParametrizable::createConfigured<TransitionModel, ConfiguredModel>(context_, std::move(depModel), std::move(newModel));
    
    // assign model to branche id
    std::vector<uint> vId=process_.getNodesWithModel(nMod);
    for (auto id:vId)
      modelmap.emplace(id, modelNode);
  }

  // assign models to branches
  treeNode_->setBranchModels(modelmap);

  ///////////////////////////
  // rootFrequencies node

  auto root = process_.getRootFrequenciesSet();

  if (root)
  {
    auto rootFreqs = std::unique_ptr<FrequenciesSet>(root->clone());
    const auto& rootParams= root->getParameters();
    
    std::vector<NodeRef> depRoot;
    for (size_t i=0;i<rootParams.size();i++)
      depRoot.push_back(ConfiguredParameter::create(context_, {dynamic_cast<ConfiguredParameter*>(getSharedParameter(rootParams[i].getName()).get())->dependency(0)}, rootParams[i]));
    
    rootFreqsNode_ = ConfiguredParametrizable::createConfigured<FrequenciesSet, ConfiguredFrequenciesSet>(context_, std::move(depRoot), std::move(rootFreqs));
  }

  // get a modelNode from the map
  modelNode_ = modelmap.begin()->second;
  
}

void LikelihoodCalculation::setNumericalDerivateConfiguration(double delta, const NumericalDerivativeType& config)
{
  auto deltaNode = NumericMutable<double>::create(context_, delta);

  if (ratesNode_)
  {
    ratesNode_->config.delta = deltaNode;
    ratesNode_->config.type = config;
  }

  /////////////////
  // model nodes

  vector<shared_ptr<BrRef> > vpn=treeNode_->getAllEdges();

  for (auto& it: vpn)
  {
    auto mN=it->getModel();
    mN->config.delta = deltaNode;
    mN->config.type = config;
  }

  ///////////////////////////
  // rootFrequencies node

  if (rootFreqsNode_)
  {
    rootFreqsNode_->config.delta = deltaNode;
    rootFreqsNode_->config.type = config;
  }
}

void LikelihoodCalculation::setClockLike(double rate)
{
  Parameter pRate("BrLen_rate", rate, Parameter::R_PLUS_STAR);

  auto rateNode = ConfiguredParameter::create(context_, pRate);

  auto rateRef= ValueFromConfiguredParameter::create(context_, {rateNode});

  /////////////////
  // brlen nodes

  vector<shared_ptr<BrRef> > vpn=treeNode_->getAllEdges();

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


/****************************************
 * Construction methods
 ****************************************/

void LikelihoodCalculation::makeRootFreqs_()
{
  
// Set root frequencies 

  size_t nbState = getStateMap().getNumberOfModelStates();
  rFreqs_ = rootFreqsNode_?ConfiguredParametrizable::createVector<ConfiguredFrequenciesSet, FrequenciesFromFrequenciesSet> (
    context_, {rootFreqsNode_}, rowVectorDimension (Eigen::Index (nbState))):
    ConfiguredParametrizable::createVector<ConfiguredModel, EquilibriumFrequenciesFromModel> (
      context_, {modelNode_}, rowVectorDimension (Eigen::Index (nbState)));
}


void LikelihoodCalculation::makeForwardLikelihoodTree_()
{
  // Build conditional likelihoods up to root recursively.
  if (!treeNode_->isRooted ()) {
    throw Exception ("PhyloTree must be rooted");
  }
        
  if (ratesNode_)
  {
    auto brlenmap = createBrLenMap(context_, *treeNode_);
      
    uint nbCat=(uint)ratesNode_->getTargetValue()->getNumberOfCategories();

    vRateCatTrees_.resize(nbCat);

    for (uint nCat=0; nCat<nbCat; nCat++)
    {
      auto catRef = CategoryFromDiscreteDistribution::create(context_, {ratesNode_}, nCat);
  
      auto rateBrLen = multiplyBrLenMap(context_, *treeNode_, catRef);
      auto treeCat = std::shared_ptr<PhyloTree_BrRef>(new PhyloTree_BrRef(*treeNode_, std::move(rateBrLen)));
      vRateCatTrees_[nCat].phyloTree=treeCat;      

      auto flt=std::make_shared<ForwardLikelihoodTree>(context_, treeCat, getStateMap());
      flt->initialize(*getShrunkData());
      vRateCatTrees_[nCat].flt=flt;
    }
  }
  else
  {
    vRateCatTrees_.resize(1);
    vRateCatTrees_[0].phyloTree=treeNode_;
    
    auto flt=std::make_shared<ForwardLikelihoodTree >(context_, treeNode_, modelNode_->getTargetValue()->getStateMap());
    flt->initialize(*getShrunkData());
    vRateCatTrees_[0].flt=flt;
  }
}


void LikelihoodCalculation::makeLikelihoodAtRoot_()
{
  if (vRateCatTrees_.size()==0)
    makeForwardLikelihoodTree_();

  std::size_t nbSite = getShrunkData()->getNumberOfSites();    
  auto rootIndex = treeNode_->getRootIndex();
  
  // Set root frequencies
  if (rFreqs_==0)
    makeRootFreqs_();
  
  if (ratesNode_)
  {
    std::vector<std::shared_ptr<Node>> vLogRoot;
      
    for (auto rateCat: vRateCatTrees_)
    {
      vLogRoot.push_back(LikelihoodFromRootConditional::create (
                           context_, {rFreqs_, rateCat.flt->getForwardLikelihoodArray(rootIndex)}, rowVectorDimension (Eigen::Index (nbSite))));
    }

    auto catProb = ProbabilitiesFromDiscreteDistribution::create(context_, {ratesNode_});
    
    vLogRoot.push_back(catProb);
      
    siteLikelihoods_ = CWiseMean<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>, Eigen::RowVectorXd>::create(context_, std::move(vLogRoot), rowVectorDimension (Eigen::Index(nbSite)));
  }
  else
  {
    siteLikelihoods_ = LikelihoodFromRootConditional::create (
      context_, {rFreqs_, vRateCatTrees_[0].flt->getForwardLikelihoodArray(rootIndex)}, rowVectorDimension (Eigen::Index (nbSite)));
  }

  auto totalLogLikelihood =
    SumOfLogarithms<Eigen::RowVectorXd>::create (context_, {siteLikelihoods_, rootWeights_}, rowVectorDimension (Eigen::Index (nbSite)));

// We want -log(likelihood)
  likelihood_ = CWiseNegate<double>::create (context_, {totalLogLikelihood}, Dimension<double> ());
}


void LikelihoodCalculation::makeLikelihoodsAtNode_(uint nodeId)
{
  if (vRateCatTrees_.size()==0)
    makeForwardLikelihoodTree_();
  
  if (rFreqs_==0)
    makeRootFreqs_();
  
  const auto& stateMap = getStateMap();
  size_t nbSite = getShrunkData()->getNumberOfSites();
  size_t nbState = stateMap.getNumberOfModelStates();
  MatrixDimension likelihoodMatrixDim = conditionalLikelihoodDimension (nbState, nbSite);
  
  auto one=ConstantOne<Eigen::RowVectorXd>::create(context_, rowVectorDimension (Eigen::Index (nbState)));

  ValueRef<Eigen::RowVectorXd> siteLikelihoodsNode;

  if (ratesNode_)
  {
    std::vector<std::shared_ptr<Node>> vLogRoot;
      
    for (auto rateCat: vRateCatTrees_)
    {
      if (!rateCat.blt)
        rateCat.blt=std::make_shared<BackwardLikelihoodTree>(context_, rateCat.flt, rateCat.phyloTree, rFreqs_, stateMap, nbSite);

      if (!rateCat.lt)
        rateCat.lt=std::make_shared<ConditionalLikelihoodTree>(treeNode_->getGraph());
      
      auto condAbove = rateCat.blt->getBackwardLikelihoodArray(nodeId);
      
      auto condBelow = rateCat.flt->getForwardLikelihoodArray(nodeId);

      auto cond = BuildConditionalLikelihood::create (
        context_, {condAbove, condBelow}, likelihoodMatrixDim);
      
      rateCat.lt->associateNode(cond, treeNode_->getNodeGraphid(treeNode_->getNode(nodeId)));
      rateCat.lt->setNodeIndex(cond, nodeId);

      auto siteLikelihoodsCat = LikelihoodFromRootConditional::create (
        context_, {one, cond}, rowVectorDimension (Eigen::Index (nbSite)));
      
      vLogRoot.push_back(std::move(siteLikelihoodsCat));
    }

    auto catProb = ProbabilitiesFromDiscreteDistribution::create(context_, {ratesNode_});
    vLogRoot.push_back(catProb);
      
    siteLikelihoodsNode = CWiseMean<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>, Eigen::RowVectorXd>::create(context_, std::move(vLogRoot), rowVectorDimension (Eigen::Index(nbSite)));
  }
  else
  {
    auto& rateCat = vRateCatTrees_[0];
    
    if (!rateCat.blt)
      rateCat.blt=std::make_shared<BackwardLikelihoodTree>(context_, rateCat.flt, rateCat.phyloTree, rFreqs_, stateMap, nbSite);

    if (!rateCat.lt)
      rateCat.lt=std::make_shared<ConditionalLikelihoodTree>(treeNode_->getGraph());
    auto condAbove = rateCat.blt->getBackwardLikelihoodArray(nodeId);
      
    auto condBelow = rateCat.flt->getForwardLikelihoodArray(nodeId);

    auto cond = BuildConditionalLikelihood::create (
      context_, {condAbove, condBelow}, likelihoodMatrixDim);
    
    rateCat.lt->associateNode(cond, treeNode_->getNodeGraphid(treeNode_->getNode(nodeId)));
    rateCat.lt->setNodeIndex(cond, nodeId);

    siteLikelihoodsNode = LikelihoodFromRootConditional::create (
      context_, {one, cond}, rowVectorDimension (Eigen::Index (nbSite)));
  }
//   auto totalLogLikelihood =
//     SumOfLogarithms<Eigen::RowVectorXd>::create (context_, {siteLikelihoodsNode}, rowVectorDimension (Eigen::Index (nbSite)));

// // We want -log(likelihood)
//   auto totalNegLogLikelihood =
//     CWiseNegate<double>::create (context_, {totalLogLikelihood}, Dimension<double> ());

//   return totalNegLogLikelihood;
}

