// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Io/Pasta.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Legacy/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Legacy/OptimizationTools.h>

#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/SimpleSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>

#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

#include <iostream>


using namespace bpp;
using namespace std;

void fitModelHSR(
    shared_ptr<SubstitutionModelInterface> model,
    shared_ptr<DiscreteDistributionInterface> rdist,
    const Tree& tree, const ParametrizablePhyloTree&  new_tree,
    shared_ptr<const ProbabilisticSiteContainerInterface> sites,
    double initialValue, double finalValue)
{
  RHomogeneousTreeLikelihood tl(tree, *sites, model, rdist, false, false, false);
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

  ApplicationTools::startTimer();

  auto process = make_shared<RateAcrossSitesSubstitutionProcess>(
        shared_ptr<SubstitutionModelInterface>(model->clone()),
        shared_ptr<DiscreteDistributionInterface>(rdist->clone()),
        shared_ptr<ParametrizablePhyloTree>(new_tree.clone()));

  Context context;
  auto lik = make_shared<LikelihoodCalculationSingleProcess>(context, sites, process);
  auto newTl = make_shared<SingleProcessPhyloLikelihood>(context, lik);


  cout << "NewTL:    " << setprecision(20) << newTl->getValue() << endl;
  cout << "NewTL D1: " << setprecision(20) << newTl->getFirstOrderDerivative("BrLen2") << endl;
  cout << "NewTL D2: " << setprecision(20) << newTl->getSecondOrderDerivative("BrLen2") << endl;
  ApplicationTools::displayResult("* initial likelihood", newTl->getValue());
  if (abs(newTl->getValue() - initialValue) > 0.001)
    throw Exception("Incorrect initial value.");
  cout << endl;

  for (size_t i = 0; i < n; ++i)
  {
    ApplicationTools::displayGauge(i, n - 1);
    newTl->matchParametersValues(pl1);
    newTl->getValue();
    newTl->matchParametersValues(pl2);
    newTl->getValue();
  }

  cout << endl;
  ApplicationTools::displayTime("New Likelihood: model upgrade");

  ApplicationTools::startTimer();
  for (size_t i = 0; i < n; ++i)
  {
    ApplicationTools::displayGauge(i, n - 1);
    newTl->matchParametersValues(pl3);
    newTl->getValue();
    newTl->matchParametersValues(pl4);
    newTl->getValue();
  }
  cout << endl;

  ApplicationTools::displayTime("New Likelihood: brlen upgrade");

  cout << endl;

  cout << "==========================================" << endl;
  cout << "==========================================" << endl;
  cout << endl;

  cout << "Optimization : " << endl;
  cout << endl;

  unsigned int nboptim = 1000;

  auto tlop = make_shared<RHomogeneousTreeLikelihood>(tree, *sites,
        shared_ptr<SubstitutionModelInterface>(model->clone()),
        shared_ptr<DiscreteDistributionInterface>(rdist->clone()),
        false, false);
  tlop->initialize();

  LegacyOptimizationTools::optimizeNumericalParameters2(
      tlop, tlop->getParameters(),
      0, 0.000001, nboptim, 0, 0);
  cout << setprecision(20) << tlop->getValue() << endl;
  ApplicationTools::displayResult("* lnL after full optimization (old)", tlop->getValue());
  if (abs(tlop->getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value.");
  tlop->getParameters().printParameters(cout);


  process.reset(new RateAcrossSitesSubstitutionProcess(model, std::shared_ptr<DiscreteDistributionInterface>(rdist->clone()), std::shared_ptr<ParametrizablePhyloTree>(new_tree.clone())));
  lik.reset(new LikelihoodCalculationSingleProcess(context, sites, process));
  newTl.reset(new SingleProcessPhyloLikelihood(context, lik));

  ParameterList opln1 = process->getBranchLengthParameters(true);

  OptimizationTools::OptimizationOptions optopt;

  optopt.parameters = newTl->getParameters();
  optopt.nbEvalMax = nboptim;
  optopt.verbose = 0;

  OptimizationTools::optimizeNumericalParameters2(newTl, optopt);
  
  cout << setprecision(20) << newTl->getValue() << endl;
  ApplicationTools::displayResult("* lnL after full optimization (new)", newTl->getValue());
  if (abs(newTl->getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value.");
  newTl->getParameters().printParameters(cout);
}


int main()
{
  unique_ptr<TreeTemplate<Node>> tree(TreeTemplateTools::parenthesisToTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);"));
  vector<string> seqNames = tree->getLeavesNames();
  vector<int> ids = tree->getNodesId();

  Newick reader;
  unique_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);", false, "", false, false));

  // -------------

  shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET;
  shared_ptr<const NucleicAlphabet> nucAlphabet = AlphabetTools::DNA_ALPHABET;

  Pasta pasta;

  auto sites = make_shared<ProbabilisticVectorSiteContainer>(alphabet);
  pasta.readAlignment("example1.pa", *sites);

  auto model = make_shared<T92>(nucAlphabet, 3.);
  auto rdist = make_shared<ConstantRateDistribution>();
  try
  {
    cout << "Testing Single Tree Traversal likelihood class..." << endl;
    fitModelHSR(model, rdist, *tree, *pTree, sites, 222.26297478, 215.7976882);
  }
  catch (Exception& ex)
  {
    cerr << ex.what() << endl;
    return 1;
  }

  return 0;
}
