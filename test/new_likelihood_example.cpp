//
// File: new_phyl_dataflow.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-19 00:00:00
// Last modified: 2017-06-08
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

  This software is a computer program whose purpose is to provide classes
  for numerical calculus. This file is part of the Bio++ project.

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

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

//#define ENABLE_OLD
#define ENABLE_NEW
#define ENABLE_DF

// Common stuff
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Function/SimpleNewtonMultiDimensions.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/FrequenciesSet/NucleotideFrequenciesSet.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <chrono>

// Old likelihood
#ifdef ENABLE_OLD
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#endif
// Newlik
#ifdef ENABLE_NEW
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>
#endif
// DF
#ifdef ENABLE_DF
#include <Bpp/NewPhyl/FrequenciesSet.h>
#include <Bpp/NewPhyl/Model.h>
#include <Bpp/NewPhyl/PhyloTree_BrRef.h>
#include <Bpp/NewPhyl/Parametrizable.h>
#include <Bpp/NewPhyl/Likelihood.h>
#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/HomogeneousLikelihoodExample.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Text/TextTools.h>
#endif

static bool enableDotOutput = false;
using namespace std;

static void dotOutput(const std::string& testName, const std::vector<const bpp::dataflow::Node*>& nodes)
{
  if (enableDotOutput)
  {
    using bpp::dataflow::DotOptions;
    writeGraphToDot(
      "debug_" + testName + ".dot", nodes, DotOptions::DetailedNodeInfo | DotOptions::ShowDependencyIndex);
  }
}

namespace
{
  using TimePoint = typename std::chrono::high_resolution_clock::time_point;
  TimePoint timingStart(void) { return std::chrono::high_resolution_clock::now(); }
  void timingEnd(TimePoint start, const std::string& prefix)
  {
    auto end = timingStart(); // ill named, just to get now()
    std::cout << "[time-ns] " << prefix << " "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << "\n";
  }
  void printLik(double logLik, const std::string& prefix)
  {
    std::cout << "[log-lik] " << prefix << " " << logLik << "\n";
  }

  template<typename Func>
  void do_func_multiple_times(const std::string& timePrefix, Func f)
  {
    constexpr std::size_t updatesNbIterations = 1000;
    auto ts = timingStart();
    for (std::size_t i = 0; i < updatesNbIterations; ++i)
      f();
    timingEnd(ts, timePrefix);
  }
  void do_param_changes_multiple_times(bpp::DerivableSecondOrder& llh,
                                       const std::string& timePrefix,
                                       const bpp::ParameterList& p1,
                                       const bpp::ParameterList& p2)
  {

    llh.matchParametersValues(p1);
    printLik(llh.getValue(), timePrefix);
    llh.matchParametersValues(p2);
    printLik(llh.getValue(), timePrefix);

    do_func_multiple_times(timePrefix, [&]() {
      llh.matchParametersValues(p1);
      llh.getValue();
      llh.matchParametersValues(p2);
      llh.getValue();
    });
  }

  void optimize_for_params(bpp::DerivableSecondOrder& llh,
                           const std::string& prefix,
                           const bpp::ParameterList& params)
  {
    auto ts = timingStart();
    bpp::ConjugateGradientMultiDimensions optimizer(&llh);
    // bpp::SimpleNewtonMultiDimensions optimizer(&llh);
    optimizer.setVerbose(0);
    optimizer.setProfiler(nullptr);       // bpp::ApplicationTools::message);
    optimizer.setMessageHandler(nullptr); // bpp::ApplicationTools::message);
    optimizer.setMaximumNumberOfEvaluations(100);
    optimizer.getStopCondition()->setTolerance(0.000001);
    optimizer.setConstraintPolicy(bpp::AutoParameter::CONSTRAINTS_AUTO);
    optimizer.init(params);
    optimizer.optimize();
    timingEnd(ts, prefix);
    printLik(llh.getValue(), prefix);
  }

  struct CommonStuff
  {
    const bpp::NucleicAlphabet& alphabet;
    bpp::VectorSiteContainer sites;
    const char* treeStr;
    bpp::ParameterList paramModel1;
    bpp::ParameterList paramModel2;
    bpp::ParameterList paramBrLen1;
    bpp::ParameterList paramBrLen2;
    bpp::ParameterList paramRoot1;
    bpp::ParameterList paramRoot2;

    CommonStuff()
      : alphabet(bpp::AlphabetTools::DNA_ALPHABET)
      , sites(&alphabet)
      , treeStr("(((A:0.01, B:0.02):0.03,C:0.01):0.3,D:0.1);")
    {
      // Init sequences
      sites.addSequence(
        bpp::BasicSequence("A", "ATCCAGACATGCCGGGACTTTGCAGAGAAGGAGTTGTTTCCCATTGCAGCCCAGGTGGATAAGGAACAGC", &alphabet));
      sites.addSequence(
        bpp::BasicSequence("B", "CGTCAGACATGCCGTGACTTTGCCGAGAAGGAGTTGGTCCCCATTGCGGCCCAGCTGGACAGGGAGCATC", &alphabet));
      sites.addSequence(
        bpp::BasicSequence("C", "GGTCAGACATGCCGGGAATTTGCTGAAAAGGAGCTGGTTCCCATTGCAGCCCAGGTAGACAAGGAGCATC", &alphabet));
      sites.addSequence(
        bpp::BasicSequence("D", "TTCCAGACATGCCGGGACTTTACCGAGAAGGAGTTGTTTTCCATTGCAGCCCAGGTGGATAAGGAACATC", &alphabet));

      // Set of parameters to apply to tree + model
      paramModel1.addParameter(bpp::Parameter("T92.kappa", 0.1));
      paramModel2.addParameter(bpp::Parameter("T92.kappa", 0.2));
      paramRoot1.addParameter(bpp::Parameter("GC.theta", 0.1));
      paramRoot2.addParameter(bpp::Parameter("GC.theta", 0.2));
      paramBrLen1.addParameter(bpp::Parameter("BrLen1", 0.1));
      paramBrLen2.addParameter(bpp::Parameter("BrLen1", 0.2));
    }
  };
}

#ifdef ENABLE_OLD
TEST_CASE("old")
{
  const CommonStuff c;
  auto ts = timingStart();
  auto model = new bpp::T92(&c.alphabet, 3.);
  auto distribution = new bpp::ConstantDistribution(1.0);
  auto tree = std::unique_ptr<bpp::TreeTemplate<bpp::Node>>(bpp::TreeTemplateTools::parenthesisToTree(c.treeStr));
  bpp::RHomogeneousTreeLikelihood llh(*tree, c.sites, model, distribution, false, false);
  timingEnd(ts, "old_setup");
  ts = timingStart();
  llh.initialize();
  auto logLik = llh.getValue();
  timingEnd(ts, "old_init_value");
  printLik(logLik, "old_init_value");

  std::cout << "[dbrlen1] " << llh.getFirstOrderDerivative("BrLen1") << "\n";
  do_param_changes_multiple_times(llh, "old_param_model_change", c.paramModel1, c.paramModel2);
  do_param_changes_multiple_times(llh, "old_param_brlen_change", c.paramBrLen1, c.paramBrLen2);
  optimize_for_params(llh, "old_brlens_opt", llh.getBranchLengthsParameters());
}
#endif

#ifdef ENABLE_NEW
TEST_CASE("new")
{
  const CommonStuff c;
  auto ts = timingStart();
  auto model = new bpp::T92(&c.alphabet, 3.);
  auto rootFreqs = new bpp::GCFrequenciesSet(&c.alphabet, 0.1);
  auto distribution = new bpp::GammaDiscreteRateDistribution(3, 0.5);
  
  bpp::Newick reader;
  auto phyloTree = std::unique_ptr<bpp::PhyloTree>(reader.parenthesisToPhyloTree(c.treeStr, false, "", false, false));
  auto paramPhyloTree = new bpp::ParametrizablePhyloTree(*phyloTree);
  std::vector<std::string> globalParameterNames;

  auto process =
    std::unique_ptr<bpp::NonHomogeneousSubstitutionProcess>(bpp::NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(model, distribution, rootFreqs, paramPhyloTree, globalParameterNames));
  
  auto likelihoodCompStruct = std::unique_ptr<bpp::RecursiveLikelihoodTreeCalculation>(
    new bpp::RecursiveLikelihoodTreeCalculation(c.sites, process.get(), false, true));
  bpp::SingleProcessPhyloLikelihood llh(process.get(), likelihoodCompStruct.release());
  timingEnd(ts, "new_setup");
  ts = timingStart();
  llh.computeLikelihood();
  auto logLik = llh.getValue();
  timingEnd(ts, "new_init_value");
  printLik(logLik, "new_init_value");

  std::cout << "[dbrlen1] " << llh.getFirstOrderDerivative("BrLen1") << "\n";
  // do_param_changes_multiple_times(llh, "new_param_model_change", c.paramModel1, c.paramModel2);
  // do_param_changes_multiple_times(llh, "new_param_root_change", c.paramRoot1, c.paramRoot2);
  // do_param_changes_multiple_times(llh, "new_param_brlen_change", c.paramBrLen1, c.paramBrLen2);
  optimize_for_params(llh, "new_brlens_opt", llh.getBranchLengthParameters());
}
#endif

TEST_CASE("df")
{
  const CommonStuff c;
  bpp::dataflow::Context context;

  auto ts = timingStart();

  // Rate

//  auto rates = std::unique_ptr<bpp::GammaDiscreteRateDistribution>(new bpp::GammaDiscreteRateDistribution(3, 0.5));
  auto rates = std::unique_ptr<bpp::ConstantRateDistribution>(new bpp::ConstantRateDistribution());

  auto ratesParameters = bpp::dataflow::createParameterMap(context, *rates);

  auto depvecRates=bpp::dataflow::createDependencyVector(
    *rates, [&ratesParameters](const std::string& paramName) { return ratesParameters[paramName]; });

  auto ratesNode = bpp::dataflow::ConfiguredParametrizable::createConfigured<bpp::DiscreteDistribution, bpp::dataflow::ConfiguredDistribution>(
    context,
    std::move(depvecRates),
    std::move(rates));

  auto deltaRate = bpp::dataflow::NumericMutable<double>::create(context, 0.001);
  ratesNode->config.delta = deltaRate;
  ratesNode->config.type = bpp::dataflow::NumericalDerivativeType::ThreePoints;

  
  // Read tree structure
  bpp::Newick reader;
  
  auto phyloTree = std::unique_ptr<bpp::PhyloTree>(reader.parenthesisToPhyloTree(c.treeStr, false, "", false, false));

  auto treeBrLen = bpp::dataflow::createBrLenMap(context, *phyloTree);

  auto treeNode = std::shared_ptr<bpp::dataflow::PhyloTree_BrRef>(new bpp::dataflow::PhyloTree_BrRef(*phyloTree, treeBrLen));

  auto parTree = std::unique_ptr<bpp::ParametrizablePhyloTree>(new bpp::ParametrizablePhyloTree(*phyloTree));
  
  // Model: create simple leaf nodes as model parameters
  auto model = std::unique_ptr<bpp::T92>(new bpp::T92(&c.alphabet, 3.));
  auto modelParameters = bpp::dataflow::createParameterMap(context, *model);
  
  auto depvecModel=bpp::dataflow::createDependencyVector(
    *model, [&modelParameters](const std::string& paramName) { return modelParameters[paramName]; });

  auto modelNode = bpp::dataflow::ConfiguredParametrizable::createConfigured<bpp::TransitionModel, bpp::dataflow::ConfiguredModel>(
    context,
    std::move(depvecModel),
    std::move(model));

  auto delta = bpp::dataflow::NumericMutable<double>::create(context, 0.001);
  modelNode->config.delta = delta;
  modelNode->config.type = bpp::dataflow::NumericalDerivativeType::ThreePoints;


  // auto model2 = std::unique_ptr<bpp::T92>(new bpp::T92(&c.alphabet, 2.));
  // auto model2Parameters = bpp::dataflow::createParameterMap(context, *model2);
  
  // auto depvecModel2=bpp::dataflow::createDependencyVector(
  //   *model2, [&model2Parameters](const std::string& paramName) { return model2Parameters[paramName]; });

  // auto model2Node = bpp::dataflow::ConfiguredParametrizable::createConfigured<bpp::TransitionModel, bpp::dataflow::ConfiguredModel>(
  //   context,
  //   std::move(depvecModel2),
  //   std::move(model2));

  // auto delta2 = bpp::dataflow::NumericMutable<double>::create(context, 0.001);
  // model2Node->config.delta = delta2;
  // model2Node->config.type = bpp::dataflow::NumericalDerivativeType::ThreePoints;

  bpp::dataflow::ModelMap modelmap;
  auto vId=treeNode->getAllEdgesIndexes();

  for (const auto& id: vId)
  {
    // if (id%2==0)
    //   modelmap.emplace(id, model2Node);
    // else
      modelmap.emplace(id, modelNode);
  }

  treeNode->setBranchModels(modelmap);
  
  // Root Frequencies
  auto rootFreqs = std::unique_ptr<bpp::GCFrequenciesSet>(new bpp::GCFrequenciesSet(&c.alphabet, 0.1));
  auto rootFreqsParameters = bpp::dataflow::createParameterMap(context, *rootFreqs);
  auto depvecRootFreqs=bpp::dataflow::createDependencyVector(
    *rootFreqs, [&rootFreqsParameters](const std::string& paramName) { return rootFreqsParameters[paramName]; });

  auto rootFreqsNode = bpp::dataflow::ConfiguredParametrizable::createConfigured<bpp::FrequenciesSet, bpp::dataflow::ConfiguredFrequenciesSet>(
    context,
    std::move(depvecRootFreqs),
    std::move(rootFreqs));

  rootFreqsNode->config.delta = delta;
  rootFreqsNode->config.type = bpp::dataflow::NumericalDerivativeType::ThreePoints;

  // Build likelihood value node
  auto l = makeHomogeneousLikelihoodNodes(context, c.sites, treeNode, modelNode, rootFreqsNode, ratesNode);

  
  // Assemble a bpp::Optimizer compatible interface to HomogeneousLikelihoodNodes.
  bpp::ParameterList brlenOnlyParameters;


  for (const auto& id: vId)
  {
    auto param = bpp::dataflow::DataFlowParameter(parTree->getParameter("BrLen"+bpp::TextTools::toString(id)), dynamic_pointer_cast<bpp::dataflow::NumericMutable<double>>(treeNode->getEdge(id)->getBrLen()));
    brlenOnlyParameters.addParameter(std::move(param));
  }
  
  bpp::ParameterList allParameters;
  allParameters.addParameters(brlenOnlyParameters);
  
  // for (const auto& p : modelParameters)
  // {
  //   auto modelb=modelNode->getValue();
    
  //   auto param = bpp::dataflow::DataFlowParameter(*p.second->getValue(), Name(),()modelb->getParameter(modelb->getParameterNameWithoutNamespace(p.first)), p.second);
  //   param = std::shared_ptr<Parameter>(p.second->getValue()->clone());
  //   allParameters.addParameter(std::move(param));
  // }

  // for (const auto& p : model2Parameters)
  // {
  //   auto modelb=model2Node->getValue();
    
  //   auto param = bpp::dataflow::DataFlowParameter(modelb->getParameter(modelb->getParameterNameWithoutNamespace(p.first)), p.second);
  //   param.setName(param.getName()+"_2");
    
  //   allParameters.addParameter(std::move(param));
  // }

  // for (const auto& p : rootFreqsParameters)
  // {
  //   auto root2=rootFreqsNode->getValue();
    
  //   auto param = bpp::dataflow::DataFlowParameter(root2->getParameter(root2->getParameterNameWithoutNamespace(p.first)), p.second);
  //   allParameters.addParameter(param);
  // }

  
  bpp::dataflow::DataFlowFunction llh(context, l, allParameters);
  timingEnd(ts, "df_setup");

  ts = timingStart();
  auto logLik = llh.getValue();
  timingEnd(ts, "df_init_value");
  printLik(logLik, "df_init_value");
  dotOutput("likelihood_example_value", {l.get()});

  // Manual access to dbrlen1
  auto dlogLik_dbrlen1 = l->deriveAsValue(context, *treeNode->getEdge(1)->getBrLen());
  std::cout << "[dbrlen1] " << dlogLik_dbrlen1->getTargetValue() << "\n";
  dotOutput("likelihood_example_dbrlen1", {dlogLik_dbrlen1.get()});

  auto dlogLik_dkappa = l->deriveAsValue(context,  *modelParameters["T92.kappa"]);
  
  std::cout << "[dkappa] " << dlogLik_dkappa->getTargetValue() << "\n";
  dotOutput("likelihood_example_dkappa", {dlogLik_dkappa.get()});

  // do_param_changes_multiple_times(llh, "df_param_model_change", c.paramModel1, c.paramModel2);
  // do_param_changes_multiple_times(llh, "df_param_root_change", c.paramRoot1, c.paramRoot2);
  // do_param_changes_multiple_times(llh, "df_param_brlen_change", c.paramBrLen1, c.paramBrLen2);
  
  optimize_for_params(llh, "df_brlens_opt", brlenOnlyParameters);
  optimize_for_params(llh, "df_all_opt", allParameters);
  // TODO test optimization with model params
}

int main(int argc, char** argv)
{
  const std::string keyword = "dot_output";
  for (int i = 1; i < argc; ++i)
  {
    if (argv[i] == keyword)
    {
      enableDotOutput = true;
    }
  }
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  return context.run();
}
