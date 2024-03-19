// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Legacy/Likelihood/RHomogeneousTreeLikelihood.h>

#include <Bpp/Phyl/Legacy/OptimizationTools.h>

#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/SimpleSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>

#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

#include <iostream>


using namespace bpp;
using namespace std;

void fitModelHSR(std::shared_ptr<SubstitutionModelInterface> model,
    std::shared_ptr<DiscreteDistributionInterface> rdist,
    const Tree& tree,
    std::shared_ptr<ParametrizablePhyloTree> partree,
    std::shared_ptr<const SiteContainerInterface> sites,
    double initialValue, double finalValue)
{
  ApplicationTools::startTimer();

  RHomogeneousTreeLikelihood tl(tree, *sites, model, rdist, false, false);
  tl.initialize();

  ApplicationTools::displayResult("Test model", model->getName());
  cout << "OldTL: " << setprecision(20) << tl.getValue() << endl;
  cout << "OldTL D1: " << setprecision(20) << tl.getFirstOrderDerivative("BrLen2") << endl;
  cout << "OldTL D2: " << setprecision(20) << tl.getSecondOrderDerivative("BrLen2") << endl;
  ApplicationTools::displayResult("* initial likelihood", tl.getValue());
  if (abs(tl.getValue() - initialValue) > 0.001)
    throw Exception("Incorrect initial value.");
  cout << endl;

  Parameter p1("T92.kappa", 0.1);
  Parameter p2("T92.kappa", 3);

  ParameterList pl1; pl1.addParameter(p1);
  ParameterList pl2; pl2.addParameter(p2);

  Parameter p3("BrLen1", 0.1);
  Parameter p4("BrLen1", 0.2);

  ParameterList pl3; pl3.addParameter(p3);
  ParameterList pl4; pl4.addParameter(p4);

  unsigned int n = 100000;

  ApplicationTools::startTimer();
  for (size_t i = 0; i < n; ++i)
  {
    ApplicationTools::displayGauge(i, n - 1);
    tl.matchParametersValues(pl1);
    tl.getValue();
    tl.matchParametersValues(pl2);
    tl.getValue();
  }
  cout << endl;
  ApplicationTools::displayTime("Old Likelihood: model upgrade");

  ApplicationTools::startTimer();
  for (size_t i = 0; i < n; ++i)
  {
    ApplicationTools::displayGauge(i, n - 1);
    tl.matchParametersValues(pl3);
    tl.getValue();
    tl.matchParametersValues(pl4);
    tl.getValue();
  }
  cout << endl;

  ApplicationTools::displayTime("Old Likelihood: brlen upgrade");

  cout << "=============================" << endl;

  auto process = std::make_shared<RateAcrossSitesSubstitutionProcess>(model, rdist, partree);

  Context context;
  auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, sites, process);
  SingleProcessPhyloLikelihood llh(context, lik);

  llh.getFirstOrderDerivative("BrLen0");
  llh.getFirstOrderDerivative("BrLen1");
  llh.getFirstOrderDerivative("BrLen2");
  llh.getFirstOrderDerivative("BrLen3");
  llh.getFirstOrderDerivative("BrLen4");

  cout << "NewTL: " << setprecision(20) << llh.getValue() << endl;
  cout << "NewTL D1: " << setprecision(20) << llh.getFirstOrderDerivative("BrLen2") << endl;
  cout << "NewTL D2: " << setprecision(20) << llh.getSecondOrderDerivative("BrLen2") << endl;

  ApplicationTools::displayResult("* initial likelihood", llh.getValue());
  if (abs(llh.getValue() - initialValue) > 0.001)
  {
    cerr << "Incorrect initial value." << endl;
    throw Exception("Incorrect initial value.");
  }
  cout << endl;
  for (size_t i = 0; i < n; ++i)
  {
    ApplicationTools::displayGauge(i, n - 1);
    llh.matchParametersValues(pl1);
    llh.getValue();
    llh.matchParametersValues(pl2);
    llh.getValue();
  }

  cout << endl;
  ApplicationTools::displayTime("New Likelihood: model upgrade");

  ApplicationTools::startTimer();
  for (size_t i = 0; i < n; ++i)
  {
    ApplicationTools::displayGauge(i, n - 1);
    llh.matchParametersValues(pl3);
    llh.getValue();
    llh.matchParametersValues(pl4);
    llh.getValue();
  }
  cout << endl;

  ApplicationTools::displayTime("New Likelihood: brlen upgrade");

  cout << endl;

  cout << "==========================================" << endl;
  cout << "==========================================" << endl;
  cout << endl;

  cout << "Optimization : " << endl;
  cout << endl;

  uint nboptim = 1000;

  ///////////////

  auto tlop = make_shared<RHomogeneousTreeLikelihood>(
        tree,
        *sites,
        shared_ptr<SubstitutionModelInterface>(model->clone()),
        shared_ptr<DiscreteDistributionInterface>(rdist->clone()),
        false,
        false);
  tlop->initialize();

  LegacyOptimizationTools::optimizeNumericalParameters2(tlop, tlop->getParameters(), 0, 0.000001, nboptim, 0, 0);
  cout << setprecision(20) << tlop->getValue() << endl;
  ApplicationTools::displayResult("* lnL after full optimization (old)", tlop->getValue());
  if (abs(tlop->getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value.");
  tlop->getParameters().printParameters(cout);


  process = std::make_shared<RateAcrossSitesSubstitutionProcess>(
        shared_ptr<SubstitutionModelInterface>(model->clone()),
        shared_ptr<DiscreteDistributionInterface>(rdist->clone()),
        partree);

  Context context2;

  lik.reset(new LikelihoodCalculationSingleProcess(context2, sites, process));

  auto llh2 = make_shared<SingleProcessPhyloLikelihood>(context2, lik);

  ParameterList opln1 = process->getBranchLengthParameters(true);

  OptimizationTools::optimizeNumericalParameters2(llh2, llh2->getParameters(), 0, 0.000001, nboptim, 0, 0);
  cout << setprecision(20) << llh2->getValue() << endl;
  ApplicationTools::displayResult("* lnL after full optimization (new)", llh2->getValue());
  if (abs(llh2->getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value.");
  llh2->getParameters().printParameters(cout);
}


int main()
{
  unique_ptr<TreeTemplate<Node>> tree(TreeTemplateTools::parenthesisToTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);"));
  vector<string> seqNames = tree->getLeavesNames();
  vector<int> ids = tree->getNodesId();

  Newick reader;
  auto pTree = unique_ptr<PhyloTree>(reader.parenthesisToPhyloTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);", false, "", false, false));
  auto paramphyloTree = make_shared<ParametrizablePhyloTree>(*pTree);

  // -------------

  shared_ptr<const NucleicAlphabet> nucAlphabet = AlphabetTools::DNA_ALPHABET;
  shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET;

  auto sites = make_shared<VectorSiteContainer>(alphabet);
  auto seqA = make_unique<Sequence>("A", "ATCCAGACATGCCGGGACTTTGCAGAGAAGGAGTTGTTTCCCATTGCAGCCCAGGTGGATAAGGAACAGC", alphabet);
  sites->addSequence("A", seqA);
  auto seqB = make_unique<Sequence>("B", "CGTCAGACATGCCGTGACTTTGCCGAGAAGGAGTTGGTCCCCATTGCGGCCCAGCTGGACAGGGAGCATC", alphabet);
  sites->addSequence("B", seqB);
  auto seqC = make_unique<Sequence>("C", "GGTCAGACATGCCGGGAATTTGCTGAAAAGGAGCTGGTTCCCATTGCAGCCCAGGTAGACAAGGAGCATC", alphabet);
  sites->addSequence("C", seqC);
  auto seqD = make_unique<Sequence>("D", "TTCCAGACATGCCGGGACTTTACCGAGAAGGAGTTGTTTTCCATTGCAGCCCAGGTGGATAAGGAACATC", alphabet);
  sites->addSequence("D", seqD);

  auto model = make_shared<T92>(nucAlphabet, 3.);
  auto rdist = make_shared<GammaDiscreteRateDistribution>(4, 1.0);
  try
  {
    cout << "Testing Single Tree Traversal likelihood class..." << endl;
    fitModelHSR(model, rdist, *tree, paramphyloTree, sites, 228.6333642493463, 198.47216106233);
  }
  catch (exception& ex)
  {
    cerr << "ERROR!!!" << endl;
    cerr << ex.what() << endl;
    return 1;
  }

  return 0;
}
