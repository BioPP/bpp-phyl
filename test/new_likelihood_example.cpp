// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/OnABranchPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationOnABranch.h>
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
        "debug_" + testName + ".dot", nodes); // , DotOptions::DetailedNodeInfo | DotOptions::ShowDependencyIndex);
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

void optimize_for_params(shared_ptr<SecondOrderDerivable> llh,
    const std::string& prefix,
    const ParameterList& params)
{
  auto ts = timingStart();
  ConjugateGradientMultiDimensions optimizer(llh);
  // SimpleNewtonMultiDimensions optimizer(llh);
  optimizer.setVerbose(0);
  optimizer.setProfiler(nullptr);       // ApplicationTools::message);
  optimizer.setMessageHandler(nullptr); // ApplicationTools::message);
  optimizer.setMaximumNumberOfEvaluations(10000);
  optimizer.getStopCondition()->setTolerance(0.000001);
  optimizer.setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);

  optimizer.init(params);
  optimizer.optimize();
  timingEnd(ts, prefix);
  printLik(llh->getValue(), prefix);
}

struct CommonStuff
{
  shared_ptr<const Alphabet> alphabet;
  shared_ptr<const NucleicAlphabet> nucAlphabet;
  shared_ptr<VectorSiteContainer> sites;
  const char* treeStr;
  ParameterList paramModel1;
  ParameterList paramModel2;
  ParameterList paramBrLen1;
  ParameterList paramBrLen2;
  ParameterList paramRoot1;
  ParameterList paramRoot2;

  CommonStuff()
    : alphabet(AlphabetTools::DNA_ALPHABET)
    , nucAlphabet(AlphabetTools::DNA_ALPHABET)
    , sites(new VectorSiteContainer(alphabet))
//      , treeStr("(A:0.01, B:0.02);")
    , treeStr("((A:0.01, B:0.02):0.03,C:0.01);")
    , paramModel1()
    , paramModel2()
    , paramBrLen1()
    , paramBrLen2()
    , paramRoot1()
    , paramRoot2()
  {
    // Init sequences
    auto seqA = make_unique<Sequence>("A", "ATCCAGACATGCCGGGACTTTGCAGAGAAGGAGTTGTTTCCCATTGCAGCCCAGGTGGATAAGGAACAGC", alphabet);
    sites->addSequence("A", seqA);
    auto seqB = make_unique<Sequence>("B", "CGTCAGACATGCCGTGACTTTGCCGAGAAGGAGTTGGTCCCCATTGCGGCCCAGCTGGACAGGGAGCATC", alphabet);
    sites->addSequence("B", seqB);
    auto seqC = make_unique<Sequence>("C", "GGTCAGACATGCCGGGAATTTGCTGAAAAGGAGCTGGTTCCCATTGCAGCCCAGGTAGACAAGGAGCATC", alphabet);
    sites->addSequence("C", seqC);
    auto seqD = make_unique<Sequence>("D", "TTCCAGACATGCCGGGACTTTACCGAGAAGGAGTTGTTTTCCATTGCAGCCCAGGTGGATAAGGAACATC", alphabet);
    sites->addSequence("D", seqD);

    // Set of parameters to apply to tree + model
    paramModel1.addParameter(Parameter("T92.kappa", 0.1));
    paramModel2.addParameter(Parameter("T92.kappa", 0.2));
    paramRoot1.addParameter(Parameter("GC.theta", 0.1));
    paramRoot2.addParameter(Parameter("GC.theta", 0.2));
    paramBrLen1.addParameter(Parameter("BrLen1", 0.1));
    paramBrLen2.addParameter(Parameter("BrLen1", 0.2));
  }

  CommonStuff(const CommonStuff& cs)
    : alphabet(cs.alphabet)
    , nucAlphabet(cs.nucAlphabet)
    , sites(cs.sites)
    , treeStr(cs.treeStr)
    , paramModel1(cs.paramModel1)
    , paramModel2(cs.paramModel2)
    , paramBrLen1(cs.paramBrLen1)
    , paramBrLen2(cs.paramBrLen2)
    , paramRoot1(cs.paramRoot1)
    , paramRoot2(cs.paramRoot2)
  {}

  CommonStuff& operator=(const CommonStuff& cs)
  {
    alphabet = cs.alphabet;
    nucAlphabet = cs.nucAlphabet;
    sites = cs.sites;
    treeStr = cs.treeStr;
    paramModel1 = cs.paramModel1;
    paramModel2 = cs.paramModel2;
    paramBrLen1 = cs.paramBrLen1;
    paramBrLen2 = cs.paramBrLen2;
    paramRoot1 = cs.paramRoot1;
    paramRoot2 = cs.paramRoot2;
    return(*this);
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

  // Read tree structure
  Newick reader;
  auto phyloTree = std::shared_ptr<PhyloTree>(reader.parenthesisToPhyloTree(c.treeStr, false, "", false, false));


  // Models
  auto t92 = std::make_unique<T92>(c.nucAlphabet, 2., 0.5);


  auto t922 = std::make_unique<T92>(c.nucAlphabet, 3., 0.7);

  std::map<std::string, unique_ptr<DiscreteDistributionInterface>> mapParam1;

  auto sdm = std::map<double, double>({{0.1, 0.3}, {0.6, 0.7}});
  mapParam1["theta"].reset(new SimpleDiscreteDistribution(sdm));
  auto mt92 = std::make_shared<MixtureOfASubstitutionModel>(c.alphabet, std::move(t922), mapParam1);


  auto scenario = std::make_shared<ModelScenario>();

  auto mp1 = make_shared<ModelPath>();
  mp1->setModel(mt92, Vuint({0}));
  scenario->addModelPath(mp1);

  mp1 = make_shared<ModelPath>();
  mp1->setModel(mt92, Vuint({1}));
  scenario->addModelPath(mp1);


  auto rootFreqs = std::make_shared<GCFrequencySet>(c.nucAlphabet, 0.1);

  auto distribution = std::make_shared<ConstantRateDistribution>();

  std::vector<std::string> globalParameterNames({"T92.kappa"});


  shared_ptr<SubstitutionProcessInterface> process = NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(mt92, distribution, phyloTree, rootFreqs, scenario);

  process->getParameters().printParameters(cerr);

  // Build likelihoodx value node
  auto ope = make_shared<OneProcessSequenceEvolution>(process);

  auto llh = make_shared<OneProcessSequencePhyloLikelihood>(context, c.sites,  ope);

  timingEnd(ts, "df_setup");
  auto lik = llh->getLikelihoodCalculation();
  dotOutput("likelihood_example_value", {lik->getLikelihoodNode().get()});

  ts = timingStart();

  auto logLik = llh->getValue();
  timingEnd(ts, "df_init_value");
  printLik(logLik, "df_init_value");

  // Likelihoods after changing a model on a branch
  OnABranchPhyloLikelihood obp(context, llh->getLikelihoodCalculationSingleProcess(), 2);

  auto ct92 = ConfiguredParametrizable::createConfigured<BranchModelInterface, ConfiguredModel>(context, std::move(t92));

  obp.setModel(ct92);

  dotOutput("likelihood_example_onabranch", {obp.likelihoodCalculation().getLikelihoodNode().get()});

  auto logLikob = obp.getValue();

  printLik(logLikob, "branch 2 init_value");

  // Manual access to dbrlen
  auto br = dynamic_cast<ConfiguredParameter*>(lik->hasParameter("BrLen1") ? lik->getParameter("BrLen1").get() : lik->getParameter("BrLen_rate").get());

  auto dlogLik_dbrlen1 = lik->getLikelihoodNode()->deriveAsValue(context, *br->dependency(0));

  dotOutput("likelihood_example_dbrlen1", {dlogLik_dbrlen1.get()});
  std::cout << "[dbrlen1] " << dlogLik_dbrlen1->targetValue() << "\n";
  std::cout << "[dbrlen1] " << llh->getFirstOrderDerivative("BrLen1") << std::endl;
  std::cout << "[dbrlen1] " << llh->getSecondOrderDerivative("BrLen1") << std::endl;

  // // Manual access to dkappa

  auto kappa = dynamic_cast<ConfiguredParameter*>(llh->likelihoodCalculation().getParameter("T92.kappa").get());
  auto dlogLik_dkappa = lik->getLikelihoodNode()->deriveAsValue(context, *kappa->dependency(0));
  std::cout << "[dkappa] " << dlogLik_dkappa->targetValue() << "\n";
  dotOutput("likelihood_example_dkappa", {dlogLik_dkappa.get()});

  auto d2logLik_dkappa2 = dlogLik_dkappa->deriveAsValue(context, *kappa->dependency(0));
  std::cout << "[d2kappa] " << d2logLik_dkappa2->targetValue() << "\n";
  dotOutput("likelihood_example_dkappa2", {d2logLik_dkappa2.get()});

  // Manual access to dalpha

  // auto alpha= dynamic_cast<ConfiguredParameter*>(llh.getLikelihoodCalculation()->getParameter("Gamma.alpha").get());
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

//  l->getParameters().printParameters(cerr);

  optimize_for_params(llh, "df_all_opt", llh->getParameters());
  dotOutput("likelihood_optim_value", {lik->getLikelihoodNode().get()});
  llh->getParameters().printParameters(std::cerr);
}
