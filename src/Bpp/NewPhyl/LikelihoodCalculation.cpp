//
// File: LikelihoodCalculation.cpp
// Authors: François Gindraud, Laurent Guéguen (2018)
// Created: jeudi 28 février 2019, à 07h 22
//

#include "Bpp/NewPhyl/LikelihoodCalculation.h"
#include "Bpp/NewPhyl/ForwardLikelihoodTree.h"

#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>

#include <unordered_map>
#include <list>

using namespace std;
using namespace bpp;
using namespace bpp::dataflow;

LikelihoodCalculation::LikelihoodCalculation(dataflow::Context & context,
                                             const AlignedValuesContainer & sites,
                                             const SubstitutionProcess& process):
  AbstractParametrizable(""),
  likelihood_(), context_(context), process_(process), psites_(&sites),
  treeNode_(), modelNode_(), rootFreqsNode_(),
  ratesNode_()
{
  makeProcessNodes_();
}

LikelihoodCalculation::LikelihoodCalculation(const LikelihoodCalculation& lik) :
  AbstractParametrizable(lik),
  likelihood_(), context_(), process_(lik.process_), psites_(lik.psites_),
  treeNode_(), modelNode_(), rootFreqsNode_(),
  ratesNode_()
{
  makeProcessNodes_();
}

LikelihoodCalculation::LikelihoodCalculation(dataflow::Context & context,
                                             const SubstitutionProcess& process):
  AbstractParametrizable(""),
  likelihood_(), context_(context), process_(process), psites_(),
  treeNode_(), modelNode_(), rootFreqsNode_(),
  ratesNode_()
{
  makeProcessNodes_();
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
    std::vector<bpp::dataflow::NodeRef> deps;
    for (size_t i=0;i<rateParams.size();i++)
      deps.push_back(bpp::dataflow::ConfiguredParameter::create(context_, {dynamic_cast<ConfiguredParameter*>(getSharedParameter(rateParams[i].getName()).get())->dependency(0)}, rateParams[i]));
    
    ratesNode_ = bpp::dataflow::ConfiguredParametrizable::createConfigured<bpp::DiscreteDistribution, bpp::dataflow::ConfiguredDistribution>(
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
    
    mapBr.emplace(parTree.getEdgeIndex(branch),
                  bpp::dataflow::ConfiguredParameter::create(context_, {dynamic_cast<ConfiguredParameter*>(getSharedParameter(bp.getName()).get())->dependency(0)}, bp));
  }
  
  treeNode_ = std::shared_ptr<bpp::dataflow::PhyloTree_BrRef>(new bpp::dataflow::PhyloTree_BrRef(parTree, mapBr));

  /////////////////
  // model nodes
  
  std::vector<std::shared_ptr<ConfiguredModel>> vModelNodes;

  size_t nMod=process_.getNumberOfModels();
        
  bpp::dataflow::ModelMap modelmap;

  for (size_t nm=0; nm<nMod;nm++)
  {
    std::unique_ptr<TransitionModel> newModel(process_.getModel(nm)->clone());
    const auto& modelParams=process_.getModel(nm)->getParameters();
    
    std::vector<bpp::dataflow::NodeRef> depModel;
    for (size_t np=0;np<modelParams.size();np++)
    {
      if (hasParameter(modelParams[np].getName()+"_"+ TextTools::toString(nm+1)))
        depModel.push_back(bpp::dataflow::ConfiguredParameter::create(context_, {dynamic_cast<ConfiguredParameter*>(getSharedParameter(modelParams[np].getName()+"_"+ TextTools::toString(nm+1)).get())->dependency(0)}, modelParams[np]));
      else
      {
        if (nMod==1)
          depModel.push_back(bpp::dataflow::ConfiguredParameter::create(context_, {dynamic_cast<ConfiguredParameter*>(getSharedParameter(modelParams[np].getName()).get())->dependency(0)}, modelParams[np]));
        else
          throw ParameterNotFoundException("LikelihoodCalculation::makeProcessNodes_", modelParams[np].getName());
      }
    }  
    auto modelNode = bpp::dataflow::ConfiguredParametrizable::createConfigured<bpp::TransitionModel, bpp::dataflow::ConfiguredModel>(context_, std::move(depModel), std::move(newModel));
    
    // assign model to branche id
    std::vector<uint> vId=process_.getNodesWithModel(nm);
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
    
    std::vector<bpp::dataflow::NodeRef> depRoot;
    for (size_t i=0;i<rootParams.size();i++)
      depRoot.push_back(bpp::dataflow::ConfiguredParameter::create(context_, {dynamic_cast<ConfiguredParameter*>(getSharedParameter(rootParams[i].getName()).get())->dependency(0)}, rootParams[i]));
    
    rootFreqsNode_ = bpp::dataflow::ConfiguredParametrizable::createConfigured<bpp::FrequenciesSet, bpp::dataflow::ConfiguredFrequenciesSet>(context_, std::move(depRoot), std::move(rootFreqs));
  }

  // get a modelNode from the map
  modelNode_ = modelmap.begin()->second;
  
}

void LikelihoodCalculation::setNumericalDerivateConfiguration(double delta, const NumericalDerivativeType& config)
{
  auto deltaNode = bpp::dataflow::NumericMutable<double>::create(context_, delta);

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

  auto rateNode = bpp::dataflow::ConfiguredParameter::create(context_, pRate);

  auto rateRef= bpp::dataflow::ValueFromConfiguredParameter::create(context_, {rateNode});

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


ValueRef<double> LikelihoodCalculation::makeLikelihoodNodes_()
{
  std::size_t nbSite = psites_->getNumberOfSites();
    
  const auto nbState = modelNode_->getTargetValue()->getNumberOfStates (); // Number of stored state values !
  // Build conditional likelihoods up to root recursively.
  if (!treeNode_->isRooted ()) {
    throw Exception ("PhyloTree must be rooted");
  }
        
  // Set root frequencies 
  auto rFreqs = rootFreqsNode_?dataflow::ConfiguredParametrizable::createVector<dataflow::ConfiguredFrequenciesSet, dataflow::FrequenciesFromFrequenciesSet> (
    context_, {rootFreqsNode_}, rowVectorDimension (Eigen::Index (nbState))):
    dataflow::ConfiguredParametrizable::createVector<dataflow::ConfiguredModel, dataflow::EquilibriumFrequenciesFromModel> (
      context_, {modelNode_}, rowVectorDimension (Eigen::Index (nbState)));

  dataflow::ValueRef<Eigen::RowVectorXd> siteLikelihoods;

  if (ratesNode_)
  {
    auto brlenmap = bpp::dataflow::createBrLenMap(context_, *treeNode_);
      
    uint nbCat=(uint)ratesNode_->getTargetValue()->getNumberOfCategories();

    std::vector<std::shared_ptr<bpp::dataflow::Node>> vLogRoot;
      
    for (uint nCat=0; nCat<nbCat; nCat++)
    {
      auto catRef = bpp::dataflow::CategoryFromDiscreteDistribution::create(context_, {ratesNode_}, nCat);
  
      auto rateBrLen = bpp::dataflow::multiplyBrLenMap(context_, *treeNode_, catRef);

      auto treeCat = std::shared_ptr<bpp::dataflow::PhyloTree_BrRef>(new bpp::dataflow::PhyloTree_BrRef(*treeNode_, std::move(rateBrLen)));

      bpp::dataflow::ForwardLikelihoodTree flt(context_, treeCat, modelNode_->getTargetValue()->getStateMap());
      flt.initialize(*psites_);
      auto rootConditionalLikelihoods=flt.getForwardLikelihoodArray(treeCat->getRootIndex ());

      auto siteLikelihoodsCat = dataflow::LikelihoodFromRootConditional::create (
        context_, {rFreqs, rootConditionalLikelihoods}, rowVectorDimension (Eigen::Index (nbSite)));

      vLogRoot.push_back(std::move(siteLikelihoodsCat));
    }

    auto catProb = bpp::dataflow::ProbabilitiesFromDiscreteDistribution::create(context_, {ratesNode_});

    vLogRoot.push_back(catProb);
      
    siteLikelihoods = bpp::dataflow::CWiseMean<Eigen::RowVectorXd, bpp::dataflow::ReductionOf<Eigen::RowVectorXd>, Eigen::RowVectorXd>::create(context_, std::move(vLogRoot), rowVectorDimension (Eigen::Index(nbSite)));
  }
  else
  {
    bpp::dataflow::ForwardLikelihoodTree flt(context_, treeNode_, modelNode_->getTargetValue()->getStateMap());
    flt.initialize(*psites_);

    auto rootConditionalLikelihoods=flt.getForwardLikelihoodArray(treeNode_->getRootIndex ());

    siteLikelihoods = dataflow::LikelihoodFromRootConditional::create (
      context_, {rFreqs, rootConditionalLikelihoods}, rowVectorDimension (Eigen::Index (nbSite)));
  }
    
  auto totalLogLikelihood =
    dataflow::SumOfLogarithms<Eigen::RowVectorXd>::create (context_, {siteLikelihoods}, rowVectorDimension (Eigen::Index (nbSite)));

// We want -log(likelihood)
  auto totalNegLogLikelihood =
    dataflow::CWiseNegate<double>::create (context_, {totalLogLikelihood}, Dimension<double> ());

  return totalNegLogLikelihood;
}

