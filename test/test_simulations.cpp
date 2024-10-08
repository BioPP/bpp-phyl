// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

#include <Bpp/Seq/Io/BppOAlignmentWriterFormat.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSiteSimulator.h>
#include <Bpp/Phyl/Simulation/GivenDataSubstitutionProcessSequenceSimulator.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <iostream>

using namespace bpp;
using namespace std;

int main()
{
  Newick reader;
  auto phyloTree = std::shared_ptr<PhyloTree>(reader.parenthesisToPhyloTree("(((A:0.1, B:0.2):0.3,C:0.4):0.1,D:0.5);", false, "", false, false));

  vector<string> seqNames = phyloTree->getAllLeavesNames();
  // -------------

  shared_ptr<const NucleicAlphabet> nucAlphabet = AlphabetTools::DNA_ALPHABET;
  auto model = std::make_shared<T92>(nucAlphabet, 3.);
  auto rdist = std::make_shared<ConstantRateDistribution>();
  auto rootFreqs = std::make_shared<GCFrequencySet>(nucAlphabet);
  std::vector<std::string> globalParameterNames({"T92.kappa"});

  shared_ptr<SubstitutionProcessInterface> process_sim = NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(model, rdist, phyloTree, rootFreqs, globalParameterNames);
  shared_ptr<SubstitutionProcessInterface> process = NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(model, rdist, phyloTree, rootFreqs, globalParameterNames);
  shared_ptr<SubstitutionProcessInterface> process2 = NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(model, rdist, phyloTree, rootFreqs, globalParameterNames);

  vector<double> thetas;
  for (unsigned int i = 0; i < process_sim->getNumberOfModels(); ++i)
  {
    double theta = RandomTools::giveRandomNumberBetweenZeroAndEntry(0.8) + 0.1;
    cout << "Theta" << i + 1 << " set to " << theta << endl;
    process_sim->setParameterValue("T92.theta_" + TextTools::toString(i + 1), theta);
    thetas.push_back(theta);
  }

  double theta= RandomTools::giveRandomNumberBetweenZeroAndEntry(0.8) + 0.1;
  cout << "Theta root set to " << theta << endl;
  process_sim->setParameterValue("GC.theta", theta);

  SimpleSubstitutionProcessSiteSimulator simulator(process_sim);

  unsigned int n = 100000;
  auto profiler  = make_shared<StlOutputStream>(make_unique<ofstream>("profile.txt", ios::out));
  auto messenger = make_shared<StlOutputStream>(make_unique<ofstream>("messages.txt", ios::out));

  // Check fast simulation first:

  cout << "Fast check:" << endl;

  // Generate data set:
  auto sites = make_shared<VectorSiteContainer>(seqNames, nucAlphabet);
  for (unsigned int i = 0; i < n; ++i)
  {
    auto simSite = simulator.simulateSite();
    simSite->setCoordinate(static_cast<int>(i));
    sites->addSite(simSite, false);
  }

  cout << "fit model" << endl;

  // Now fit model:
  Context context;
  auto l = std::make_shared<LikelihoodCalculationSingleProcess>(context, sites, process);
  auto llh = make_shared<SingleProcessPhyloLikelihood>(context, l);

  OptimizationTools::OptimizationOptions optopt;

  optopt.parameters = llh->getParameters();
  optopt.verbose = 0;
  optopt.messenger = messenger;
  optopt.profiler = profiler;
  
  OptimizationTools::optimizeNumericalParameters2(llh, optopt);
  
  process->matchParametersValues(llh->getParameters());

  // Now compare estimated values to real ones:
  for (size_t i = 0; i < thetas.size(); ++i)
  {
    cout << thetas[i] << "\t" << process->model(i + 1).parameter("theta").getValue() << endl;
    double diff = abs(thetas[i] - process->model(i + 1).parameter("theta").getValue());
    if (diff > 0.1)
    {
      cout << "difference too large" << endl;
      return 1;
    }
  }

  // Now try detailed simulations:

  cout << "Detailed check:" << endl;

  // Generate data set:
  auto sites2 = make_shared<VectorSiteContainer>(seqNames, nucAlphabet);
  for (unsigned int i = 0; i < n; ++i)
  {
    auto result = simulator.dSimulateSite();
    auto simSite = result->getSite();
    unique_ptr<Site> site(dynamic_cast<Site*>(simSite.release()));
    site->setCoordinate(static_cast<int>(i));
    sites2->addSite(site, false);
  }

  // Now fit model:

  auto l2 = make_shared<LikelihoodCalculationSingleProcess>(context, sites2, process2);
  auto llh2 = make_shared<SingleProcessPhyloLikelihood>(context, l2);

  optopt.parameters = llh2->getParameters();
 
  OptimizationTools::optimizeNumericalParameters2(llh2, optopt);
 
  process2->matchParametersValues(llh2->getParameters());

  // Now compare estimated values to real ones:
  for (size_t i = 0; i < thetas.size(); ++i)
  {
    cout << thetas[i] << "\t" << process2->model(i + 1).parameter("theta").getValue() << endl;
    double diff = abs(thetas[i] - process2->model(i + 1).parameter("theta").getValue());
    if (diff > 0.1)
    {
      cout << "difference too large" << endl;
      return 1;
    }
  }

  cout << "Estimates fine." << endl;

  // -------------


  GivenDataSubstitutionProcessSequenceSimulator gdps(llh2->getLikelihoodCalculationSingleProcess());

  auto vec2 = gdps.simulate();

  BppOAlignmentWriterFormat bppoWriter(1);
  auto oAln = bppoWriter.read("Fasta");

  oAln->writeAlignment("seq1.fasta", *sites2, true);
  oAln->writeAlignment("seq2.fasta", *vec2, true);
  // compare

  for (const auto& name : seqNames)
  {
    const auto& seq1 = sites2->sequence(name);
    const auto& seq2 = vec2->sequence(name);

    auto sim = SiteContainerTools::computeSimilarity(seq1, seq2);
    cout << name << ":" << sim << endl;
    if (sim!=1)
    {
      cout << "Bad similarity in sequence " << name << endl;
      return 1;
    }
  }

  return 0;
}
