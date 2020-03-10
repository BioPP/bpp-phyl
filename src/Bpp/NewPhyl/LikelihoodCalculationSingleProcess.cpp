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


  ///////
  // tree node
  suff=spcm?("_"+TextTools::toString(spcm->getTreeNumber())):"";
  processNodes_.treeNode_ = makeProcessTree(context_, paramList, process_, suff);

  ///////////////////////////
  // rootFrequencies node

  auto root = process_.getRootFrequenciesSet();
  if (root)
  {
    suff=spcm?("_"+TextTools::toString(spcm->getRootFrequenciesNumber())):"";
    processNodes_.rootFreqsNode_ = ConfiguredParametrizable::createConfigured<FrequenciesSet, ConfiguredFrequenciesSet>(context_, *root, paramList, suff);
  }

  auto itE = processNodes_.treeNode_->allEdgesIterator();
  // get any modelNode from the map (only for StateMap)
  for (itE->start();!itE->end();itE->next())
  {
    if ((*(*itE))->getModel()!=0)
    {
      processNodes_.modelNode_ = (*(*itE))->getModel();
      break;
    }
  }
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

  auto rateNode = ConfiguredParameter::create(context_, pRate);

  auto rateRef= ValueFromConfiguredParameter::create(context_, {rateNode});

  /////////////////
  // brlen nodes

  vector<shared_ptr<ProcessEdge> > vpn=processNodes_.treeNode_->getAllEdges();

  for (auto& it: vpn)
  {
    auto cp=it->getBrLen();

    if (cp)
    {
      auto mulref = CWiseMul<double, std::tuple<double, double>>::create (context_, {cp->dependency(0), rateRef}, Dimension<double>());

      auto cp2=ConfiguredParameter::resetDependencies(context_, cp, {mulref});
  
      it->setBrLen(cp2);
    }
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

const Eigen::RowVectorXd& LikelihoodCalculationSingleProcess::getSiteLikelihoodsForAClass(size_t nCat, bool shrunk)
{
  if (shrunk)
    return getSiteLikelihoodsTree(nCat)->getRoot()->getTargetValue();
  else
   return expandVector(getSiteLikelihoodsTree(nCat)->getRoot())->getTargetValue();
}

const AllRatesSiteLikelihoods& LikelihoodCalculationSingleProcess::getSiteLikelihoodsForAllClasses(bool shrunk)
{
  auto nbCat=vRateCatTrees_.size();
  auto allLk=std::make_shared<AllRatesSiteLikelihoods>(nbCat,shrunk?getNumberOfDistinctSites():getNumberOfSites());
  
  for (size_t nCat=0;nCat<nbCat;nCat++)
    allLk->row(nCat)=getSiteLikelihoodsForAClass(nCat, shrunk);
  
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

// std::shared_ptr<ConditionalLikelihoodTree> LikelihoodCalculationSingleProcess::getConditionalLikelihoodTree(size_t nCat)
// {
//   if (nCat>=vRateCatTrees_.size())
//     throw Exception("LikelihoodCalculationSingleProcess::getConditionalLikelihoodTree : Bad Class number " + TextTools::toString(nCat));
  
//   if (shrunkData_ && !likelihood_)
//     makeLikelihoodAtRoot_();
  
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

    for (size_t nCat=0;nCat<vRateCatTrees_.size();nCat++)
      vLogRoot.push_back(ProbabilityFromDiscreteDistribution::create(context_, {processNodes_.ratesNode_},(uint)nCat));
    
    siteLikelihoods_ = CWiseMean<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>, ReductionOf<double>>::create(context_, std::move(vLogRoot), rowVectorDimension (Eigen::Index(nbSite)));
    
  }
  else
  {
    siteLikelihoods_ = LikelihoodFromRootConditional::create (
      context_, {rFreqs_, vRateCatTrees_[0].flt->getForwardLikelihoodArrayAtRoot()}, rowVectorDimension (Eigen::Index (nbSite)));
  }

  // likelihoods per site
  patternedSiteLikelihoods_= expandVector(siteLikelihoods_);

  auto totalLogLikelihood =
    SumOfLogarithms<Eigen::RowVectorXd>::create (context_, {siteLikelihoods_, rootWeights_}, rowVectorDimension (Eigen::Index (nbSite)));

// We want -log(likelihood)
  likelihood_ = CWiseNegate<double>::create (context_, {totalLogLikelihood}, Dimension<double> ());

  // using bpp::dataflow::DotOptions;
  // writeGraphToDot(
  //   "debug_lik.dot", {likelihood_.get()});//, DotOptions::DetailedNodeInfo | DotOp
}


ValueRef<double> LikelihoodCalculationSingleProcess::makeLikelihoodsAtNode_(uint speciesId)
{
  if (vRateCatTrees_.size()==0)
    makeForwardLikelihoodTree_();

  if (rFreqs_==0)
    makeRootFreqs_();
  
  auto processTree = processNodes_.treeNode_;

  const auto& stateMap = getStateMap();
  size_t nbSite = getShrunkData()->getNumberOfSites();
  size_t nbState = stateMap.getNumberOfModelStates();
  MatrixDimension likelihoodMatrixDim = conditionalLikelihoodDimension (nbState, nbSite);
  
  ValueRef<Eigen::RowVectorXd> siteLikelihoodsNode;

  std::shared_ptr<ConditionalLikelihood> cond(0);


  SiteLikelihoodsRef distinctSiteLikelihoodsNode;

  std::vector<std::shared_ptr<Node>> vLogRoot; // if several rates

  for (auto& rateCat: vRateCatTrees_)
  {
    if (!rateCat.blt)
      rateCat.blt=std::make_shared<BackwardLikelihoodTree>(context_, rateCat.flt, rateCat.phyloTree, rFreqs_, stateMap, nbSite);

    if (!rateCat.clt)
      rateCat.clt=std::make_shared<ConditionalLikelihoodTree>(rateCat.flt->getGraph());

    if (!rateCat.lt)
      rateCat.lt=std::make_shared<SiteLikelihoodsTree>(processTree->getGraph());

    auto& dagIndexes = rateCat.flt->getDAGNodesIndexes(speciesId);

    std::vector<std::shared_ptr<Node>> vCond;
    
    for (const auto& index : dagIndexes)
    {
      auto condAbove = rateCat.blt->getBackwardLikelihoodArray(index);
      auto condBelow = rateCat.flt->getForwardLikelihoodArray(index);
      
      cond = BuildConditionalLikelihood::create (
        context_, {condAbove, condBelow}, likelihoodMatrixDim);
      
      if (dagIndexes.size()>1) // for sum 
        vCond.push_back(cond);
      
      rateCat.clt->associateNode(cond, rateCat.flt->getNodeGraphid(rateCat.flt->getNode(index)));
      rateCat.clt->setNodeIndex(cond, index);
    }

    /*
     * If several DAG nodes related with this species node, sum the
     * likelihoods of all (already multiplied by their probability).
     *
     */

    if (dagIndexes.size()>1)
      cond = CWiseAdd<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>::create(context_, std::move(vCond), likelihoodMatrixDim);
    
    auto one=ConstantOne<Eigen::RowVectorXd>::create(context_, rowVectorDimension (Eigen::Index (nbState)));
    
    auto siteLikelihoodsCat = LikelihoodFromRootConditional::create (
      context_, {one, cond}, rowVectorDimension (Eigen::Index (nbSite)));
    
    rateCat.lt->associateNode(siteLikelihoodsCat, processTree->getNodeGraphid(processTree->getNode(speciesId)));
    rateCat.lt->setNodeIndex(siteLikelihoodsCat, speciesId);

    if (!processNodes_.ratesNode_)
    {
      distinctSiteLikelihoodsNode = siteLikelihoodsCat;
      break;
    }
    vLogRoot.push_back(siteLikelihoodsCat);
  }

  
  if (processNodes_.ratesNode_)
  {
    auto catProb = ProbabilitiesFromDiscreteDistribution::create(context_, {processNodes_.ratesNode_});
    vLogRoot.push_back(catProb);
    
    distinctSiteLikelihoodsNode = CWiseMean<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>, Eigen::RowVectorXd>::create(context_, std::move(vLogRoot), rowVectorDimension (Eigen::Index(nbSite)));
  }

  // And sum all ll
  auto totalLogLikelihood =
    SumOfLogarithms<Eigen::RowVectorXd>::create (context_, {distinctSiteLikelihoodsNode, rootWeights_}, rowVectorDimension (Eigen::Index (nbSite)));

  // We want -log(likelihood)
  auto totalNegLogLikelihood =
    CWiseNegate<double>::create (context_, {totalLogLikelihood}, Dimension<double> ());
   return totalNegLogLikelihood;
}

