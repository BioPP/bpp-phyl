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

// Common stuff
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Function/SimpleNewtonMultiDimensions.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/MixtureOfASubstitutionModel.h>
#include <Bpp/Phyl/Model/MultinomialFromTransitionModel.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <chrono>


#include <Bpp/Phyl/Likelihood/DataFlow/BackwardLikelihoodTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/OneProcessSequencePhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include "Bpp/Phyl/Likelihood/SubstitutionProcess.h"
#include <Bpp/Text/TextTools.h>

static bool enableDotOutput = true;
using namespace std;
using namespace bpp;

static void dotOutput(const std::string& testName, const std::vector<const Node_DF*>& nodes)
{
  if (enableDotOutput)
  {
    using bpp::DotOptions;
    writeGraphToDot(
      "debug_" + testName + ".dot", nodes);//, DotOptions::DetailedNodeInfo | DotOptions::ShowDependencyIndex);
  }
}

using namespace std;

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

  // template<typename Func>
  // void do_func_multiple_times(const std::string& timePrefix, Func f)
  // {
  //   constexpr std::size_t updatesNbIterations = 1000;
  //   auto ts = timingStart();
  //   for (std::size_t i = 0; i < updatesNbIterations; ++i)
  //     f();
  //   timingEnd(ts, timePrefix);
  // }
  // void do_param_changes_multiple_times(DerivableSecondOrder& llh,
  //                                      const std::string& timePrefix,
  //                                      const ParameterList& p1,
  //                                      const ParameterList& p2)
  // {

  //   llh.matchParametersValues(p1);
  //   printLik(llh.getValue(), timePrefix);
  //   llh.matchParametersValues(p2);
  //   printLik(llh.getValue(), timePrefix);

  //   do_func_multiple_times(timePrefix, [&]() {
  //     llh.matchParametersValues(p1);
  //     llh.getValue();
  //     llh.matchParametersValues(p2);
  //     llh.getValue();
  //   });
  // }

  void optimize_for_params(DerivableSecondOrder& llh,
                           const std::string& prefix,
                           const ParameterList& params)
  {
    auto ts = timingStart();
    ConjugateGradientMultiDimensions optimizer(&llh);
    // SimpleNewtonMultiDimensions optimizer(&llh);
    optimizer.setVerbose(0);
    optimizer.setProfiler(nullptr);       // ApplicationTools::message);
    optimizer.setMessageHandler(nullptr); // ApplicationTools::message);
    optimizer.setMaximumNumberOfEvaluations(10000);
    optimizer.getStopCondition()->setTolerance(0.000001);
    optimizer.setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);

    optimizer.init(params);
    optimizer.optimize();
    timingEnd(ts, prefix);
    printLik(llh.getValue(), prefix);
  }

  struct CommonStuff
  {
    const NucleicAlphabet& alphabet;
    VectorSiteContainer sites;
    const char* treeStr;
    ParameterList paramModel1;
    ParameterList paramModel2;
    ParameterList paramBrLen1;
    ParameterList paramBrLen2;
    ParameterList paramRoot1;
    ParameterList paramRoot2;

    CommonStuff()
      : alphabet(AlphabetTools::DNA_ALPHABET)
      , sites(&alphabet)
      , treeStr("(A:0.01, B:0.02);")
//      , treeStr("((A:0.01, B:0.02):0.03,C:0.01);")
    {
      // Init sequences
      sites.addSequence(
        BasicSequence("A", "ATCCAGACATGCCGGGACTTTGCAGAGAAGGAGTTGTTTCCCATTGCAGCCCAGGTGGATAAGGAACAGC", &alphabet));
      sites.addSequence(
        BasicSequence("B", "CGTCAGACATGCCGTGACTTTGCCGAGAAGGAGTTGGTCCCCATTGCGGCCCAGCTGGACAGGGAGCATC", &alphabet));
      sites.addSequence(
        BasicSequence("C", "GGTCAGACATGCCGGGAATTTGCTGAAAAGGAGCTGGTTCCCATTGCAGCCCAGGTAGACAAGGAGCATC", &alphabet));
      sites.addSequence(
        BasicSequence("D", "TTCCAGACATGCCGGGACTTTACCGAGAAGGAGTTGTTTTCCATTGCAGCCCAGGTGGATAAGGAACATC", &alphabet));

      // Set of parameters to apply to tree + model
      paramModel1.addParameter(Parameter("T92.kappa", 0.1));
      paramModel2.addParameter(Parameter("T92.kappa", 0.2));
      paramRoot1.addParameter(Parameter("GC.theta", 0.1));
      paramRoot2.addParameter(Parameter("GC.theta", 0.2));
      paramBrLen1.addParameter(Parameter("BrLen1", 0.1));
      paramBrLen2.addParameter(Parameter("BrLen1", 0.2));
    }
  };
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

  const CommonStuff c;

  Context context;

  auto ts = timingStart();
  
  auto t92 = std::make_shared<T92>(&c.alphabet, 3., 0.7);
  auto t922 = std::make_shared<T92>(&c.alphabet, 3., 0.7);

  std::map<std::string, DiscreteDistribution*> mapParam1;

  mapParam1["kappa"]=new GammaDiscreteDistribution(2, 1);

//  auto mt92 = std::make_shared<MixtureOfASubstitutionModel>(&c.alphabet, t92.get(), mapParam1);

  auto k80 = std::make_shared<K80>(&c.alphabet, 2.);
  // std::map<std::string, DiscreteDistribution*> mapParam2;
  // auto sdm=std::map<double, double>({{0.0001,0.3},{200.,0.7}});
  
  // mapParam2["kappa"]=new SimpleDiscreteDistribution(sdm);
  // auto mk80 = std::make_shared<MixtureOfASubstitutionModel>(&c.alphabet, k80.get(), mapParam2);

  /* scenario
     auto scenario = std::make_shared<ModelScenario>();

     auto mp1=make_shared<ModelPath>();
     //mp1->setModel(mt92,Vuint({0}));
     mp1->setModel(mk80,Vuint({0}));
     scenario->addModelPath(mp1);

     mp1=make_shared<ModelPath>();
     //mp1->setModel(mt92,Vuint({1}));
     mp1->setModel(mk80,Vuint({1}));
     scenario->addModelPath(mp1);

  */

  auto rootFreqs = std::make_shared<GCFrequencySet>(&c.alphabet, 0.1);

  auto distribution = std::make_shared<ConstantRateDistribution>();
  //auto distribution = new GammaDiscreteRateDistribution(3, 1);

  // Read tree structure
  Newick reader;
  auto phyloTree = std::shared_ptr<PhyloTree>(reader.parenthesisToPhyloTree(c.treeStr, false, "", false, false));

//  std::vector<std::string> globalParameterNames({"T92.kappa"});

  // auto process =
  //   std::unique_ptr<NonHomogeneousSubstitutionProcess>(NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(model.get(), distribution, paramPhyloTree, rootFreqs, globalParameterNames));


//  process->setModelScenario(scenario);

  
  // auto process  = std::make_shared<NonHomogeneousSubstitutionProcess>(distribution, paramPhyloTree, rootFreqs);

  // process->addModel(k80, Vuint({0,1,3}));

  // process->addModel(t92, Vuint({2}));
    
  auto process  = NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(k80, distribution, phyloTree, rootFreqs);//, scenario));

  process->getParameters().printParameters(cerr);
  
  // Build likelihood value node
  auto l = std::make_shared<LikelihoodCalculationSingleProcess>(context, c.sites, *process);

  
  l->setNumericalDerivateConfiguration(0.001, NumericalDerivativeType::ThreePoints);
//  l->setClockLike();

  OneProcessSequenceEvolution ope(*process);

  OneProcessSequencePhyloLikelihood llh(context, c.sites,  ope);

  timingEnd(ts, "df_setup");
  auto lik = llh.getLikelihoodCalculation();
  dotOutput("likelihood_example_value", {lik->getLikelihoodNode().get()});

  ts = timingStart();
  
  auto logLik = llh.getValue();
  timingEnd(ts, "df_init_value");
  printLik(logLik, "df_init_value");

  // Manual access to dbrlen
  auto br= dynamic_cast<ConfiguredParameter*>(lik->hasParameter("BrLen1")?lik->getSharedParameter("BrLen1").get():lik->getSharedParameter("BrLen_rate").get());

  auto dlogLik_dbrlen1 = lik->getLikelihoodNode()->deriveAsValue(context, *br->dependency(0));

  dotOutput("likelihood_example_dbrlen1", {dlogLik_dbrlen1.get()});
  std::cout << "[dbrlen1] " << dlogLik_dbrlen1->getTargetValue() << "\n";
  std::cout << "[dbrlen1] " << llh.getFirstOrderDerivative("BrLen1") << std::endl;
  std::cout << "[dbrlen1] " << llh.getSecondOrderDerivative("BrLen1") << std::endl;

  // // Manual access to dkappa
  
  auto kappa= dynamic_cast<ConfiguredParameter*>(llh.getLikelihoodCalculation()->getSharedParameter("K80.kappa").get());
  auto dlogLik_dkappa = lik->getLikelihoodNode()->deriveAsValue(context, *kappa->dependency(0));
  std::cout << "[dkappa] " << dlogLik_dkappa->getTargetValue() << "\n";
  dotOutput("likelihood_example_dkappa", {dlogLik_dkappa.get()});
  
  // auto d2logLik_dkappa2 = dlogLik_dkappa->deriveAsValue(context, *kappa->dependency(0));
  // std::cout << "[d2kappa] " << d2logLik_dkappa2->getTargetValue() << "\n";
  // dotOutput("likelihood_example_dkappa2", {d2logLik_dkappa2.get()});

  // Manual access to dalpha

  // auto alpha= dynamic_cast<ConfiguredParameter*>(llh.getLikelihoodCalculation()->getSharedParameter("Gamma.alpha").get());
  // auto dlogLik_dalpha = lik->getLikelihood()->deriveAsValue(context, *alpha->dependency(0));
  // std::cout << "[dalpha] " << dlogLik_dalpha->getTargetValue() << "\n";
  // dotOutput("likelihood_example_dalpha", {dlogLik_dalpha.get()});
  // */
  
  // for (size_t pos=0;pos<c.sites.getNumberOfSites();pos++)
  //   std::cout << pos << " : " << lik->getLikelihoodForASite(pos) << std::endl;


  // Test on nodes
  // auto lik2 = lik->getLogLikelihoodAtNode(2);
  // std::cout << "[lik2] " << lik2->getTargetValue() << "\n";
  // dotOutput("likelihood_2", {lik2.get()});
  
  l->getParameters().printParameters(cerr);
  
  optimize_for_params(llh, "df_all_opt", l->getParameters());
  dotOutput("likelihood_optim_value", {lik->getLikelihoodNode().get()});
  llh.getParameters().printParameters(std::cerr);  

}
