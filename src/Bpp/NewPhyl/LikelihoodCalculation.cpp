//
// File: LikelihoodCalculation.cpp
// Authors: François Gindraud, Laurent Guéguen (2018)
// Created: jeudi 28 février 2019, à 07h 22
//

#include "Bpp/NewPhyl/LikelihoodCalculation.h"

#include "DataFlowWrappers.h"
#include "Bpp/NewPhyl/ForwardLikelihoodTree.h"

#include <unordered_map>
#include <list>

using namespace bpp;
using namespace bpp::dataflow;

LikelihoodCalculation::LikelihoodCalculation(dataflow::Context & context,
                                             const AlignedValuesContainer & sites,
                                             const SubstitutionProcess& process):
  likelihood_(), context_(context), process_(process), psites_(&sites),
  parameters_(), deltaNode_(), treeNode_(), modelNode_(), rootFreqsNode_(),
  ratesNode_()
{
  setUpParameters_();
  likelihood_=makeLikelihoodNodes_(context, *psites_, treeNode_, modelNode_, rootFreqsNode_, ratesNode_);
}

LikelihoodCalculation::LikelihoodCalculation(dataflow::Context & context,
                                             const SubstitutionProcess& process):
  likelihood_(), context_(context), process_(process), psites_(&sites),
  parameters_(), deltaNode_(), treeNode_(), modelNode_(), rootFreqsNode_(),
  ratesNode_()
{
  setUpParameters_();
}


void LikelihoodCalculation::setUpParameters_()
{
  double delta_=0.001;
  
  auto deltaNode_ = bpp::dataflow::NumericMutable<double>::create(context, delta_);

  // add Independent Parameters
  const auto& paramProc=process_.getIndependentParameters();

  for (size_t i=0;i<paramProc.size();i++)
    parameters_.shareParameter(ConfiguredParameter::create(context, paramProc[i]));

  // Share dependencies with aliased parameters

  for (size_t i=0;i<paramProc.size();i++)
  {
    auto vs=process_.getAlias(paramProc[i].getName());
    auto dep=dynamic_cast<ConfiguredParameter*>(&parameters_[i])->dependency(0);
    for (const auto& s:vs)
    {
      auto newacp = ConfiguredParameter::create(context, {dep}, process_.getParameter(s));
      parameters_.shareParameter(newacp);
    }
  }
  
  // rates nodes
  auto rates = process_.getRateDistribution();
  
  if (rates)
  {
    std::unique_ptr<DiscreteDistribution> newRates(rates->clone());
    const auto& rateParams=newRates->getParameters();
    std::vector<bpp::dataflow::NodeRef> deps;
    for (size_t i=0;i<rateParams.size();i++)
      deps.push_back(bpp::dataflow::ConfiguredParameter::create(context, {dynamic_cast<ConfiguredParameter*>(parameters_.getSharedParameter(rateParams[i].getName()).get())->dependency(0)}, rateParams[i]));
    
    ratesNode_ = bpp::dataflow::ConfiguredParametrizable::createConfigured<bpp::DiscreteDistribution, bpp::dataflow::ConfiguredDistribution>(
      context,
      std::move(deps),
      std::move(newRates));
    
    auto deltaRate = bpp::dataflow::NumericMutable<double>::create(context, 0.001);
    ratesNode_->config.delta = deltaNode;
    ratesNode_->config.type = bpp::dataflow::NumericalDerivativeType::ThreePoints;
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
                  bpp::dataflow::ConfiguredParameter::create(context, {dynamic_cast<ConfiguredParameter*>(parameters_.getSharedParameter(bp.getName()).get())->dependency(0)}, bp));
  }
  
  treeNode_ = std::shared_ptr<bpp::dataflow::PhyloTree_BrRef>(new bpp::dataflow::PhyloTree_BrRef(parTree, mapBr));

  /////////////////
  // model nodes
  
  std::vector<std::shared_ptr<ConfiguredModel>> vModelNodes;

  size_t nMod=process.getNumberOfModels();
        
  bpp::dataflow::ModelMap modelmap;

  for (size_t nm=0; nm<nMod;nm++)
  {
    std::unique_ptr<TransitionModel> newModel(process_.getModel(nm)->clone());
    const auto& modelParams=process_.getModel(nm)->getParameters();
    
    std::vector<bpp::dataflow::NodeRef> depModel;
    for (size_t np=0;np<modelParams.size();np++)
      depModel.push_back(bpp::dataflow::ConfiguredParameter::create(context, {dynamic_cast<ConfiguredParameter*>(parameters_.getSharedParameter(modelParams[np].getName()+"_"+ TextTools::toString(nm+1)).get())->dependency(0)}, modelParams[np]));

    auto modelNode = bpp::dataflow::ConfiguredParametrizable::createConfigured<bpp::TransitionModel, bpp::dataflow::ConfiguredModel>(context, std::move(depModel), std::move(newModel));

    auto delta = bpp::dataflow::NumericMutable<double>::create(context, 0.001);
    modelNode->config.delta = deltaNode;
    modelNode->config.type = bpp::dataflow::NumericalDerivativeType::ThreePoints;

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
      depRoot.push_back(bpp::dataflow::ConfiguredParameter::create(context, {dynamic_cast<ConfiguredParameter*>(parameters_.getSharedParameter(rootParams[i].getName()).get())->dependency(0)}, rootParams[i]));
    
    rootFreqsNode_ = bpp::dataflow::ConfiguredParametrizable::createConfigured<bpp::FrequenciesSet, bpp::dataflow::ConfiguredFrequenciesSet>(context, std::move(depRoot), std::move(rootFreqs));
  
    rootFreqsNode_->config.delta = deltaNode;
    rootFreqsNode_->config.type = bpp::dataflow::NumericalDerivativeType::ThreePoints;
  }

  // get a modelNode from the map
  modelNode_ = modelmap.begin()->second;
  
}


ValueRef<double> LikelihoodCalculation::makeLikelihoodNodes_(dataflow::Context & c,
                                                                       const AlignedValuesContainer & sites,
                                                                       std::shared_ptr<PhyloTree_BrRef> tree,
                                                                       std::shared_ptr<ConfiguredModel> model,
                                                                       std::shared_ptr<ConfiguredFrequenciesSet> rootFreqs,
                                                                       std::shared_ptr<ConfiguredDistribution> rates) {

  std::size_t nbSite = sites.getNumberOfSites();
    
  const auto nbState = model->getTargetValue()->getNumberOfStates (); // Number of stored state values !
  // Build conditional likelihoods up to root recursively.
  if (!tree->isRooted ()) {
    throw Exception ("PhyloTree must be rooted");
  }
        
  // Set root frequencies 
  auto rFreqs = rootFreqs?dataflow::ConfiguredParametrizable::createVector<dataflow::ConfiguredFrequenciesSet, dataflow::FrequenciesFromFrequenciesSet> (
    c, {rootFreqs}, rowVectorDimension (Eigen::Index (nbState))):
    dataflow::ConfiguredParametrizable::createVector<dataflow::ConfiguredModel, dataflow::EquilibriumFrequenciesFromModel> (
      c, {model}, rowVectorDimension (Eigen::Index (nbState)));

  dataflow::ValueRef<Eigen::RowVectorXd> siteLikelihoods;
    
  if (rates)
  {
    auto brlenmap = bpp::dataflow::createBrLenMap(c, *tree);
      
    uint nbCat=(uint)rates->getTargetValue()->getNumberOfCategories();

    std::vector<std::shared_ptr<bpp::dataflow::Node>> vLogRoot;
      
    for (uint nCat=0; nCat<nbCat; nCat++)
    {
      auto catRef = bpp::dataflow::CategoryFromDiscreteDistribution::create(c, {rates}, nCat);
  
      auto rateBrLen = bpp::dataflow::multiplyBrLenMap(c, *tree, std::move(catRef));

      auto treeCat = std::shared_ptr<bpp::dataflow::PhyloTree_BrRef>(new bpp::dataflow::PhyloTree_BrRef(*tree, std::move(rateBrLen)));

      bpp::dataflow::ForwardLikelihoodTree flt(c, treeCat, model->getTargetValue()->getStateMap());
      flt.initialize(sites);
      auto rootConditionalLikelihoods=flt.getForwardLikelihoodArray(treeCat->getRootIndex ());

      auto siteLikelihoodsCat = dataflow::LikelihoodFromRootConditional::create (
        c, {rFreqs, rootConditionalLikelihoods}, rowVectorDimension (Eigen::Index (nbSite)));

      vLogRoot.push_back(std::move(siteLikelihoodsCat));
    }

    auto catProb = bpp::dataflow::ProbabilitiesFromDiscreteDistribution::create(c, {rates});

    vLogRoot.push_back(catProb);
      
    siteLikelihoods = bpp::dataflow::CWiseMean<Eigen::RowVectorXd, bpp::dataflow::ReductionOf<Eigen::RowVectorXd>, Eigen::RowVectorXd>::create(c, std::move(vLogRoot), rowVectorDimension (Eigen::Index(nbSite)));
  }
  else
  {
    bpp::dataflow::ForwardLikelihoodTree flt(c, tree, model->getTargetValue()->getStateMap());
    flt.initialize(sites);

    auto rootConditionalLikelihoods=flt.getForwardLikelihoodArray(tree->getRootIndex ());

    siteLikelihoods = dataflow::LikelihoodFromRootConditional::create (
      c, {rFreqs, rootConditionalLikelihoods}, rowVectorDimension (Eigen::Index (nbSite)));
  }
    
  auto totalLogLikelihood =
    dataflow::SumOfLogarithms<Eigen::RowVectorXd>::create (c, {siteLikelihoods}, rowVectorDimension (Eigen::Index (nbSite)));

// We want -log(likelihood)
  auto totalNegLogLikelihood =
    dataflow::CWiseNegate<double>::create (c, {totalLogLikelihood}, Dimension<double> ());

  return totalNegLogLikelihood;
}

