// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Io/Pasta.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/MultinomialFromTransitionModel.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/OptimizationTools.h>

#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>

#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

#include <iostream>


using namespace bpp;
using namespace std;

static bool enableDotOutput = true;
static void dotOutput(const std::string& testName, const std::vector<const Node_DF*>& nodes)
{
  if (enableDotOutput)
  {
    using bpp::DotOptions;
    writeGraphToDot(
        "debug_" + testName + ".dot", nodes); // , DotOptions::DetailedNodeInfo | DotOptions::ShowDependencyIndex);
  }
}


int main()
{
  Newick reader;
  shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((A:0.01, B:0.02):0.03,C:0.01):0.01,D:0.1);", false, "", false, false));

  // -------------

  shared_ptr<const NucleicAlphabet> nucAlphabet = AlphabetTools::DNA_ALPHABET;
  shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET;

  Pasta pasta;
  auto sites = make_shared<ProbabilisticVectorSiteContainer>(alphabet);
  pasta.readAlignment("counts.pa", *sites);


  // model
  auto t92 = make_shared<T92>(nucAlphabet, 3.);

  auto multimodel = make_shared<MultinomialFromTransitionModel>(t92); // t92 is copied there


  auto rdist = make_shared<ConstantRateDistribution>();
  auto process = make_shared<NonHomogeneousSubstitutionProcess>(rdist, pTree);


  // internal leaves
  process->addModel(t92, {2, 4}); // internal branches
  process->addModel(multimodel, {0, 1, 3, 5}); // leaves

  if (!process->isFullySetUp())
    throw Exception("test_likelihood_multinomial: process not fully set up.");

  process->getParameters().printParameters(cerr);

  process->aliasParameters("MultinomialFrom.T92.kappa_2", "T92.kappa_1");
  process->aliasParameters("MultinomialFrom.T92.theta_2", "T92.theta_1");

  process->getIndependentParameters().printParameters(cerr);

  cerr << endl;

  Context context;

  auto lik = make_shared<LikelihoodCalculationSingleProcess>(context, sites, process);
  auto newtl = make_shared<SingleProcessPhyloLikelihood>(context, lik);

  cerr << "StartLik: " << setprecision(20) << newtl->getValue() << endl;
  newtl->getParameters().printParameters(cout);
  dotOutput("lik_multinomial_value", {lik->getLikelihoodNode().get()});

  OptimizationTools::OptimizationOptions optopt;

  optopt.parameters = newtl->getParameters();
  optopt.nbEvalMax = 1000;
  optopt.verbose = 0;

  
  OptimizationTools::optimizeNumericalParameters2(newtl, optopt);
  cout << "NewLik: " << newtl->getValue() << endl;
  newtl->getParameters().printParameters(cout);

  return 0;
}
