// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/POMO.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequencySet/FrequencySet.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Phyl/Legacy/OptimizationTools.h>

#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/SimpleSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>

#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

#include <iostream>


using namespace bpp;
using namespace std;

int main()
{
  Newick reader;
  std::shared_ptr<PhyloTree> phyloTree = reader.readPhyloTree("lysozymeLarge.dnd");

  // -------------

  string nameSeq = "lysozymeLarge.fasta";
  Fasta fasta;

  auto alphaall = make_shared<AllelicAlphabet>(AlphabetTools::DNA_ALPHABET, 3);
  auto alpha = alphaall->getStateAlphabet();

  std::shared_ptr<SiteContainerInterface> sites = fasta.readAlignment(nameSeq, alpha);

  SiteContainerTools::changeGapsToUnknownCharacters(*sites);

  vector<unique_ptr<ProbabilisticSequence>> vseq;
  for (size_t ns = 0; ns < sites->getNumberOfSequences(); ++ns)
  {
    vseq.push_back(alphaall->convertFromStateAlphabet(sites->sequence(ns)));
  }

  auto sites2 = make_shared<ProbabilisticVectorSiteContainer>(alphaall);

  for (auto& seq : vseq)
  {
    sites2->addSequence(seq->getName(), seq);
  }

  auto t92 = make_unique<T92>(AlphabetTools::DNA_ALPHABET);

  auto fitness = make_unique<FullNucleotideFrequencySet>(AlphabetTools::DNA_ALPHABET);

  auto model = make_shared<POMO>(alphaall, std::move(t92), std::move(fitness));

  auto rootFreqs = make_shared<FullFrequencySet>(model->getStateMap());

  auto distribution = make_shared<ConstantRateDistribution>();

  shared_ptr<SubstitutionProcessInterface> process = NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model, distribution, phyloTree, rootFreqs);

//  auto process = make_shared<SimpleSubstitutionProcess>(model, phyloTree);

  Context context;

  auto lik = make_shared<LikelihoodCalculationSingleProcess>(context, sites2, process);

  auto llh = make_shared<SingleProcessPhyloLikelihood>(context, lik);


  cout << "NewTL: " << setprecision(20) << llh->getValue() << endl;

  // Set up optimization options
  OptimizationTools::OptimizationOptions optopt;

  optopt.messenger = std::shared_ptr<OutputStream>(new StlOutputStream(make_unique<ofstream>("messages.txt", ios::out)));
  optopt.profiler = std::shared_ptr<OutputStream>(new StlOutputStream(make_unique<ofstream>("profile.txt", ios::out)));
  optopt.profiler->setPrecision(20);

  optopt.parameters = llh->getParameters();
  optopt.useClock = true;
  optopt.verbose = 2;
  
  OptimizationTools::optimizeNumericalParameters2(llh, optopt);

  llh->getParameters().printParameters(cerr);

  cout << setprecision(20) << llh->getValue() << endl;

  return 0;
}
