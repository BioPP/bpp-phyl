//
// File: LikelihoodCalculationSingleProcess.cpp
// Authors: François Gindraud, Laurent Guéguen (2018)
// Created: jeudi 28 février 2019, à 07h 22
//

#include "Bpp/Phyl/NewLikelihood/DataFlow/LikelihoodCalculationSingleProcess.h"
#include "Bpp/Phyl/NewLikelihood/DataFlow/ForwardLikelihoodTree.h"
#include "Bpp/Phyl/NewLikelihood/DataFlow/BackwardLikelihoodTree.h"

#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/MixedTransitionModel.h>
#include "Bpp/Phyl/NewLikelihood/SubstitutionProcessCollectionMember.h"

#include <unordered_map>
#include <list>
#include <numeric>

using namespace std;
using namespace bpp;

LikelihoodCalculationSingleProcess::LikelihoodCalculationSingleProcess(Context& context,
                                                                       const AlignedValuesContainer & sites,
                                                                       const SubstitutionProcess& process):
  AlignedLikelihoodCalculation(context), process_(process), psites_(&sites),
  rootPatternLinks_(), rootWeights_(), shrunkData_(),
  processNodes_(), rFreqs_(),
  vRateCatTrees_(), condLikelihoodTree_(0)
{
  setPatterns_();
  makeProcessNodes_();

  // Default Derivate 
  setNumericalDerivateConfiguration(0.001, NumericalDerivativeType::ThreePoints);
}

LikelihoodCalculationSingleProcess::LikelihoodCalculationSingleProcess(Context & context,
                                                                       const SubstitutionProcess& process):
  AlignedLikelihoodCalculation(context),
  process_(process), psites_(),
  rootPatternLinks_(), rootWeights_(), shrunkData_(),
  processNodes_(), rFreqs_(),
  vRateCatTrees_(), condLikelihoodTree_(0)
{
  makeProcessNodes_();

  // Default Derivate 
  setNumericalDerivateConfiguration(0.001, NumericalDerivativeType::ThreePoints);
}

LikelihoodCalculationSingleProcess::LikelihoodCalculationSingleProcess(Context & context,
                                                                       const AlignedValuesContainer & sites,
                                                                       const SubstitutionProcess& process,
                                                                       ParameterList& paramList):
  AlignedLikelihoodCalculation(context),
  process_(process), psites_(&sites),
  rootPatternLinks_(), rootWeights_(), shrunkData_(),
  processNodes_(), rFreqs_(),
  vRateCatTrees_(), condLikelihoodTree_(0)
{
  setPatterns_();
  makeProcessNodes_(paramList);

  // Default Derivate 
  setNumericalDerivateConfiguration(0.001, NumericalDerivativeType::ThreePoints);
}


LikelihoodCalculationSingleProcess::LikelihoodCalculationSingleProcess(CollectionNodes& collection,
                                                                       const AlignedValuesContainer & sites,
                                                                       size_t nProcess):
  AlignedLikelihoodCalculation(collection.getContext()), process_(collection.getCollection().getSubstitutionProcess(nProcess)), psites_(&sites),
  rootPatternLinks_(), rootWeights_(), shrunkData_(),
  processNodes_(), rFreqs_(),
  vRateCatTrees_(), condLikelihoodTree_(0)
{
  setPatterns_();
  makeProcessNodes_(collection, nProcess);

  // Default Derivate 
  setNumericalDerivateConfiguration(0.001, NumericalDerivativeType::ThreePoints);
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
  setNumericalDerivateConfiguration(0.001, NumericalDerivativeType::ThreePoints);
}

void LikelihoodCalculationSingleProcess::setPatterns_()
{
  SitePatterns patterns(psites_);
  shrunkData_       = patterns.getSites();
  rootPatternLinks_ = NumericConstant<PatternType>::create(getContext_(), patterns.getIndices());
  size_t nbSites    = shrunkData_->getNumberOfSites();
  Eigen::RowVectorXi weights(nbSites);
  for (std::size_t i=0;i<nbSites;i++)
    weights(Eigen::Index(i))=patterns.getWeights()[i];
  rootWeights_ = SiteWeights::create(getContext_(), std::move(weights));
}

void LikelihoodCalculationSingleProcess::makeProcessNodes_()
{
  ParameterList paramList;
  
  // add Independent Parameters
  const auto& paramProc=process_.getIndependentParameters();
  
  for (size_t i=0;i<paramProc.size();i++)
    paramList.shareParameter(ConfiguredParameter::create(getContext_(), paramProc[i]));
  
  // Share dependencies with aliased parameters

  for (size_t i=0;i<paramProc.size();i++)
  {
    auto vs=process_.getAlias(paramProc[i].getName());
    auto dep=dynamic_cast<const ConfiguredParameter*>(&paramList[i])->dependency(0);
    for (const auto& s:vs)
    {
      auto newacp = ConfiguredParameter::create(getContext_(), {dep}, process_.getParameter(s));
      paramList.shareParameter(newacp);
    }
  }
  makeProcessNodes_(paramList);
}

void LikelihoodCalculationSingleProcess::makeProcessNodes_(ParameterList& paramList)
{
  const auto spcm=dynamic_cast<const SubstitutionProcessCollectionMember*>(&process_);

  // share process_ parameters with those of the paramList
  const auto& paramProc=process_.getParameters();
  
  for (size_t i=0;i<paramProc.size();i++)
  {
    auto name=paramProc[i].getName();
    if (!paramList.hasParameter(name))
      throw Exception("LikelihoodCalculationSingleProcess::makeProcessNodes_ : paramList does not have parameter " + name);
    auto* confPar=dynamic_cast<ConfiguredParameter*>(&paramList.getParameter(name));
    if (!confPar)
      throw Exception("LikelihoodCalculationSingleProcess::makeProcessNodes_ : parameter " + name + "is not a ConfiguredParameter.");
      
    shareParameter_(paramList.getSharedParameter(name));
  }
  
  // // Share dependencies with aliased parameters

  // for (size_t i=0;i<paramProc.size();i++)
  // {
  //   auto vs=process_.getAlias(paramProc[i].getName());
  //   auto dep=dynamic_cast<const ConfiguredParameter*>(&paramList.getParameters()[i])->dependency(0);
  //   for (const auto& s:vs)
  //   {
  //     auto newacp = ConfiguredParameter::create(getContext_(), {dep}, process_.getParameter(s));
  //     paramList.shareParameter_(newacp);
  //   }
  // }

  // rates node
  std::string suff=spcm?("_"+TextTools::toString(spcm->getRateDistributionNumber())):"";
  
  auto rates = process_.getRateDistribution(); 
  if (rates && dynamic_cast<const ConstantRateDistribution*>(rates)==nullptr)
    processNodes_.ratesNode_ = ConfiguredParametrizable::createConfigured<DiscreteDistribution, ConfiguredDistribution>(getContext_(), *rates, paramList, suff);


  ///////
  // tree node
  suff=spcm?("_"+TextTools::toString(spcm->getTreeNumber())):"";
  processNodes_.treeNode_ = ProcessTree::makeProcessTree(getContext_(), process_, paramList, suff);

  ///////////////////////////
  // rootFrequencies node

  auto root = process_.getRootFrequenciesSet();
  if (root)
  {
    suff=spcm?("_"+TextTools::toString(spcm->getRootFrequenciesNumber())):"";
    processNodes_.rootFreqsNode_ = ConfiguredParametrizable::createConfigured<FrequenciesSet, ConfiguredFrequenciesSet>(getContext_(), *root, paramList, suff);
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

void LikelihoodCalculationSingleProcess::makeProcessNodes_(CollectionNodes& collection, size_t nProc)
{
  auto& spcm = collection.getCollection().getSubstitutionProcess(nProc);

  // share process parameters with those of the collection
  const auto& paramProc=spcm.getParameters();
  
  for (size_t i=0;i<paramProc.size();i++)
  {
    auto name=paramProc[i].getName();
    if (!collection.hasParameter(name))
      throw Exception("LikelihoodCalculationSingleProcess::makeProcessNodes_ : CollectionNodes does not have parameter " + name);
    
    shareParameter_(collection.getSharedParameter(name));
  }

  // rates node
  
  const DiscreteDistribution* rates = spcm.getRateDistribution();
  
  if (dynamic_cast<const ConstantRateDistribution*>(rates)==nullptr)
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
  auto deltaNode = NumericMutable<double>::create(getContext_(), delta);

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

  auto rateNode = ConfiguredParameter::create(getContext_(), pRate);

  auto rateRef= ValueFromConfiguredParameter::create(getContext_(), {rateNode});

  /////////////////
  // brlen nodes

  vector<shared_ptr<ProcessEdge> > vpn=processNodes_.treeNode_->getAllEdges();

  for (auto& it: vpn)
  {
    auto cp=it->getBrLen();

    if (cp)
    {
      auto mulref = CWiseMul<double, std::tuple<double, double>>::create (getContext_(), {cp->dependency(0), rateRef}, Dimension<double>());

      auto cp2=ConfiguredParameter::resetDependencies(getContext_(), cp, {mulref});
  
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

Eigen::RowVectorXd LikelihoodCalculationSingleProcess::getSiteLikelihoodsForAClass(size_t nCat, bool shrunk)
{  
  if (shrunk)
    return getSiteLikelihoodsTree_(nCat)->getRoot()->getTargetValue();
  else
    return expandVector(getSiteLikelihoodsTree_(nCat)->getRoot())->getTargetValue();
}

AllRatesSiteLikelihoods LikelihoodCalculationSingleProcess::getSiteLikelihoodsForAllClasses(bool shrunk)
{
  auto nbCat=vRateCatTrees_.size();
  auto allLk=std::make_shared<AllRatesSiteLikelihoods>(nbCat,shrunk?getNumberOfDistinctSites():getNumberOfSites());

  for (size_t nCat=0;nCat<nbCat;nCat++)
    allLk->row(nCat)=getSiteLikelihoodsForAClass(nCat, shrunk);

  return *allLk;
}

// std::shared_ptr<ConditionalLikelihoodTree> LikelihoodCalculationSingleProcess::getConditionalLikelihoodTree(size_t nCat)
// {
//   if (nCat>=vRateCatTrees_.size())
//     throw Exception("LikelihoodCalculationSingleProcess::getConditionalLikelihoodTree : Bad Class number " + TextTools::toString(nCat));
  
//   if (shrunkData_ && !likelihood_)
//     makeLikelihoodsAtRoot_();
  
//   if (vRateCatTrees_[nCat].clt==0)
//     makeLikelihoodsAtNode_(getTreeNode_()->getRootIndex());
  
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
    getContext_(), {processNodes_.rootFreqsNode_}, rowVectorDimension (Eigen::Index (nbState))):
    ConfiguredParametrizable::createVector<ConfiguredModel, EquilibriumFrequenciesFromModel> (
      getContext_(), {processNodes_.modelNode_}, rowVectorDimension (Eigen::Index (nbState)));
}


void LikelihoodCalculationSingleProcess::makeForwardLikelihoodTree_()
{
  // Build conditional likelihoods up to root recursively.
  if (!processNodes_.treeNode_->isRooted ()) {
    throw Exception ("LikelihoodCalculationSingleProcess;;makeForwardLikelihoodTree_ : PhyloTree must be rooted");
  }
  
  if (processNodes_.ratesNode_)
  {
    uint nbCat=(uint)processNodes_.ratesNode_->getTargetValue()->getNumberOfCategories();

    vRateCatTrees_.resize(nbCat);

    for (uint nCat=0; nCat<nbCat; nCat++)
    {
      ValueRef<double> catRef = CategoryFromDiscreteDistribution::create(getContext_(), {processNodes_.ratesNode_}, nCat);

      auto treeCat = std::make_shared<ProcessTree>(*processNodes_.treeNode_, catRef);
      
      vRateCatTrees_[nCat].phyloTree=treeCat;      

      auto flt=std::make_shared<ForwardLikelihoodTree>(getContext_(), treeCat, getStateMap());
      flt->initialize(*getShrunkData());
      vRateCatTrees_[nCat].flt=flt;
    }
  }
  else
  {
    vRateCatTrees_.resize(1);
    vRateCatTrees_[0].phyloTree=processNodes_.treeNode_;

    auto flt=std::make_shared<ForwardLikelihoodTree >(getContext_(), processNodes_.treeNode_, processNodes_.modelNode_->getTargetValue()->getStateMap());
    flt->initialize(*getShrunkData());
    vRateCatTrees_[0].flt=flt;
  }
}


void LikelihoodCalculationSingleProcess::makeLikelihoodsAtRoot_()
{
  if (vRateCatTrees_.size()==0)
    makeForwardLikelihoodTree_();

  std::size_t nbSite = getShrunkData()->getNumberOfSites();    

  // Set root frequencies
  if (rFreqs_==0)
    makeRootFreqs_();

  if (processNodes_.ratesNode_)
  {
    std::vector<std::shared_ptr<Node_DF>> vLogRoot;
      
    for (auto& rateCat: vRateCatTrees_)
    {
      vLogRoot.push_back(LikelihoodFromRootConditional::create (
                           getContext_(), {rFreqs_, rateCat.flt->getForwardLikelihoodArrayAtRoot()}, rowVectorDimension (Eigen::Index (nbSite))));
    }

    auto catProb = ProbabilitiesFromDiscreteDistribution::create(getContext_(), {processNodes_.ratesNode_});

    for (size_t nCat=0;nCat<vRateCatTrees_.size();nCat++)
      vLogRoot.push_back(ProbabilityFromDiscreteDistribution::create(getContext_(), {processNodes_.ratesNode_},(uint)nCat));
    
    auto sL = CWiseMean<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>, ReductionOf<double>>::create(getContext_(), std::move(vLogRoot), rowVectorDimension (Eigen::Index(nbSite)));

    setSiteLikelihoods(sL, true);
  }
  else
  {
    auto sL = LikelihoodFromRootConditional::create (
      getContext_(), {rFreqs_, vRateCatTrees_[0].flt->getForwardLikelihoodArrayAtRoot()}, rowVectorDimension (Eigen::Index (nbSite)));

    setSiteLikelihoods(sL, true);
  }

  // likelihoods per site

  setSiteLikelihoods(expandVector(patternedSiteLikelihoods_), false);

  likelihood_ =
    SumOfLogarithms<Eigen::RowVectorXd>::create (getContext_(), {patternedSiteLikelihoods_, rootWeights_}, rowVectorDimension (Eigen::Index (nbSite)));

  // using bpp::DotOptions;
  // writeGraphToDot(
  //   "debug_lik.dot", {likelihood_.get()});//, DotOptions::DetailedNodeInfo | DotOp
}


void LikelihoodCalculationSingleProcess::makeLikelihoodsAtNode_(uint speciesId)
{
  // Already built
  if (condLikelihoodTree_ && condLikelihoodTree_->hasNode(speciesId))
    return;
  
  if (vRateCatTrees_.size()==0)
    makeForwardLikelihoodTree_();

  if (rFreqs_==0)
    makeRootFreqs_();
  
  const auto& stateMap = getStateMap();
  size_t nbSite = getShrunkData()->getNumberOfSites();
  size_t nbState = stateMap.getNumberOfModelStates();
  MatrixDimension likelihoodMatrixDim = conditionalLikelihoodDimension (nbState, nbSite);

  const auto& phylotree = process_.getParametrizablePhyloTree();
  
  ValueRef<Eigen::RowVectorXd> siteLikelihoodsNode;

  std::shared_ptr<ConditionalLikelihood> cond(0);

  std::vector<NodeRef> vCondRate;
  
  SiteLikelihoodsRef distinctSiteLikelihoodsNode;
  ConditionalLikelihoodRef conditionalLikelihoodsNode;
  
  std::vector<std::shared_ptr<Node_DF>> vRoot; // if several rates

  if (!condLikelihoodTree_)
    condLikelihoodTree_ = std::make_shared<ConditionalLikelihoodTree>(phylotree.getGraph());
  
  auto one=ConstantOne<Eigen::RowVectorXd>::create(getContext_(), rowVectorDimension (Eigen::Index (nbState)));
    
  for (auto& rateCat: vRateCatTrees_)
  {
    if (!rateCat.blt)
      rateCat.blt=std::make_shared<BackwardLikelihoodTree>(getContext_(), rateCat.flt, rateCat.phyloTree, rFreqs_, stateMap, nbSite);

    if (!rateCat.clt)
      rateCat.clt=std::make_shared<ConditionalLikelihoodDAG>(rateCat.flt->getGraph());

    if (!rateCat.lt)
      rateCat.lt=std::make_shared<SiteLikelihoodsTree>(phylotree.getGraph());

    auto& dagIndexes = rateCat.flt->getDAGNodesIndexes(speciesId);

    std::vector<std::shared_ptr<Node_DF>> vCond;

    for (const auto& index : dagIndexes)
    {
      if (rateCat.clt->hasNode(index))
      {
        cond = rateCat.clt->getNode(index);
        if (dagIndexes.size()>1) // for sum 
          vCond.push_back(cond);
        continue;
      }

      auto condAbove = rateCat.blt->getBackwardLikelihoodArray(index);
      auto condBelow = rateCat.flt->getForwardLikelihoodArray(index);

      cond = BuildConditionalLikelihood::create (
        getContext_(), {condAbove, condBelow}, likelihoodMatrixDim);

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
      cond = CWiseAdd<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>::create(getContext_(), std::move(vCond), likelihoodMatrixDim);

    // for Lik at Node
    auto siteLikelihoodsCat = LikelihoodFromRootConditional::create (
      getContext_(), {one, cond}, rowVectorDimension (Eigen::Index (nbSite)));

    if (!rateCat.lt->hasNode(speciesId))
    {
      rateCat.lt->associateNode(siteLikelihoodsCat, phylotree.getNodeGraphid(phylotree.getNode(speciesId)));
      rateCat.lt->setNodeIndex(siteLikelihoodsCat, speciesId);
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
    vRoot.push_back(catProb);
    vCondRate.push_back(catProb);

    distinctSiteLikelihoodsNode = CWiseMean<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>, Eigen::RowVectorXd>::create(getContext_(), std::move(vRoot), rowVectorDimension (Eigen::Index(nbSite)));

    conditionalLikelihoodsNode = CWiseMean<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>, Eigen::RowVectorXd>::create(getContext_(), std::move(vCondRate), MatrixDimension (Eigen::Index(nbState), Eigen::Index(nbSite)));
  }

  condLikelihoodTree_->associateNode(conditionalLikelihoodsNode, phylotree.getNodeGraphid(phylotree.getNode(speciesId)));
  condLikelihoodTree_->setNodeIndex(conditionalLikelihoodsNode, speciesId);
    
  // // And sum all ll
  // auto totalLogLikelihood =
  //   SumOfLogarithms<Eigen::RowVectorXd>::create (getContext_(), {distinctSiteLikelihoodsNode, rootWeights_}, rowVectorDimension (Eigen::Index (nbSite)));

  // return totalLogLikelihood;

  // We want -log(likelihood)
  // auto totalNegLogLikelihood =
  //   CWiseNegate<double>::create (getContext_(), {totalLogLikelihood}, Dimension<double> ());
  //  return totalNegLogLikelihood;

}


std::shared_ptr<SiteLikelihoodsTree> LikelihoodCalculationSingleProcess::getSiteLikelihoodsTree_(size_t nCat)
{
  if (nCat>=vRateCatTrees_.size())
    throw Exception("LikelihoodCalculationSingleProcess::getSiteLikelihoodsTree : Bad Class number " + TextTools::toString(nCat));

  if (!shrunkData_)
    throw Exception("LikelihoodCalculationSingleProcess::getSiteLikelihoodsTree : data not set.");
        
  if (!likelihood_)
    makeLikelihoodsAtRoot_();
  
  if (vRateCatTrees_[nCat].lt==0)
    makeLikelihoodsAtNode_(getTreeNode_()->getRootIndex());
  
  return vRateCatTrees_[nCat].lt;
}

