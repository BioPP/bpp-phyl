// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <iostream>

#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

using namespace bpp;
using namespace std;

void fitModelH(shared_ptr<SubstitutionModelInterface> model,
    shared_ptr<DiscreteDistributionInterface> rdist,
    shared_ptr<PhyloTree> tree,
    shared_ptr<const SiteContainerInterface> sites,
    double initialValue, double finalValue)
{
  auto process = make_shared<RateAcrossSitesSubstitutionProcess>(model, rdist, tree);

  Context context;

  auto lik = make_shared<LikelihoodCalculationSingleProcess>(context, sites, process);

  auto llh = make_shared<SingleProcessPhyloLikelihood>(context, lik);

  ApplicationTools::displayResult("Test model", model->getName());
  double initValue = llh->getValue();

  cout << setprecision(20) << initValue << endl;
  ApplicationTools::displayResult("* initial likelihood", initValue);
  if (abs(initValue - initialValue) > 0.001)
    throw Exception("Incorrect initial value:" + TextTools::toString(initValue) + "<>" + TextTools::toString(initialValue));
  shared_ptr<OutputStream> messenger(new StlOutputStream(make_unique<ofstream>("messages.txt", ios::out)));
  shared_ptr<OutputStream> profiler(new StlOutputStream(make_unique<ofstream>("profile.txt", ios::out)));
  profiler->setPrecision(20);


  OptimizationTools::optimizeNumericalParameters2(llh, llh->getParameters(), 0, 0.000001, 10000, messenger, profiler, false, true, 2, OptimizationTools::OPTIMIZATION_NEWTON);
  cout << setprecision(20) << llh->getValue() << endl;
  ApplicationTools::displayResult("* likelihood after full optimization", llh->getValue());
  llh->getParameters().printParameters(cout);
  if (abs(llh->getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value:" + TextTools::toString(llh->getValue()) + "<>" + TextTools::toString(finalValue));
}

void fitModelHClock(shared_ptr<SubstitutionModelInterface> model,
    shared_ptr<DiscreteDistributionInterface> rdist,
    shared_ptr<PhyloTree> tree,
    shared_ptr<const SiteContainerInterface> sites,
    double initialValue, double finalValue)
{
  auto process = make_shared<RateAcrossSitesSubstitutionProcess>(model, rdist, tree);

  Context context;

  auto lik = make_shared<LikelihoodCalculationSingleProcess>(context, sites, process);
  lik->setClockLike();

  auto llh = make_shared<SingleProcessPhyloLikelihood>(context, lik);

  ApplicationTools::displayResult("Test model", model->getName());
  double initValue = llh->getValue();

  cout << setprecision(20) << initValue << endl;
  ApplicationTools::displayResult("* initial likelihood", initValue);
  if (abs(initValue - initialValue) > 0.001)
    throw Exception("Incorrect initial value:" + TextTools::toString(initValue) + "<>" + TextTools::toString(initialValue));
  shared_ptr<OutputStream> messenger(new StlOutputStream(make_unique<ofstream>("messages.txt", ios::out)));
  shared_ptr<OutputStream> profiler(new StlOutputStream(make_unique<ofstream>("profile.txt", ios::out)));
  profiler->setPrecision(20);

  OptimizationTools::optimizeNumericalParameters2(llh, llh->getParameters(), 0, 0.000001, 10000, messenger, profiler);
  cout << setprecision(20) << llh->getValue() << endl;
  ApplicationTools::displayResult("* likelihood after full optimization", llh->getValue());
  llh->getParameters().printParameters(cout);
  if (abs(llh->getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value:" + TextTools::toString(llh->getValue()) + "<>" + TextTools::toString(finalValue));
}

int main()
{
  bpp::Newick reader;
  auto phyloTree = std::shared_ptr<bpp::PhyloTree>(reader.parenthesisToPhyloTree("(((A:0.01, B:0.01):0.02,C:0.03):0.01,D:0.04);", false, "", false, false));

  shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET;
  shared_ptr<const NucleicAlphabet> nucAlphabet = AlphabetTools::DNA_ALPHABET;
  auto model = make_shared<T92>(nucAlphabet, 3.);
  auto rdist = make_shared<ConstantRateDistribution>();

  auto sites = make_shared<VectorSiteContainer>(alphabet);
  auto seqA = make_unique<Sequence>("A", "AAATGGCTGTGCACGTC", alphabet);
  sites->addSequence("A", seqA);
  auto seqB = make_unique<Sequence>("B", "AACTGGATCTGCATGTC", alphabet);
  sites->addSequence("B", seqB);
  auto seqC = make_unique<Sequence>("C", "ATCTGGACGTGCACGTG", alphabet);
  sites->addSequence("C", seqC);
  auto seqD = make_unique<Sequence>("D", "CAACGGGAGTGCGCCTA", alphabet);
  sites->addSequence("D", seqD);

  try
  {
    fitModelH(
        std::shared_ptr<SubstitutionModelInterface>(model->clone()),
        std::shared_ptr<DiscreteDistributionInterface>(rdist->clone()),
        std::shared_ptr<PhyloTree>(phyloTree->clone()),
        sites,
        94.3957, 71.0564);
  }
  catch (Exception& ex)
  {
    cerr << ex.what() << endl;
    return 1;
  }

  cout << endl << endl;

  try
  {
    fitModelHClock(
        model,
        std::shared_ptr<DiscreteDistributionInterface>(rdist->clone()),
        std::shared_ptr<PhyloTree>(phyloTree->clone()),
        sites,
        94.395699, 72.7196);
  }
  catch (Exception& ex)
  {
    cerr << ex.what() << endl;
    return 1;
  }

  // -------------

  return 0;
}
