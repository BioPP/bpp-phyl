// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/SequenceTools.h>

#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Io/Newick.h>

#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Model/FrequencySet/FrequencySet.h>

#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSequenceSimulator.h>

#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>

#include <iostream>

using namespace bpp;
using namespace std;

double testBowker(const SimpleSubstitutionProcessSequenceSimulator& sim, size_t seqlen)
{
  auto sites(sim.simulate(seqlen));
  unique_ptr<BowkerTest> bTest(SequenceTools::bowkerTest(sites->sequence(0), sites->sequence(1)));
  return bTest->getPValue();
}

int main()
{
  Newick reader;
  auto phyloTree = std::shared_ptr<PhyloTree>(reader.parenthesisToPhyloTree("(A:0.02, B:0.02);", false, "", false, false));

  // First test stationnary model:
  cout << "..:: Testing with stationary model ::.." << endl;
  shared_ptr<NucleicAlphabet> alphabet(new DNA());
  auto model1 = make_shared<T92>(alphabet, 3., 0.65);
  shared_ptr<DiscreteDistributionInterface> rdist(new ConstantRateDistribution());

  shared_ptr<SubstitutionProcessInterface> process1 = NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model1, rdist, phyloTree);

  SimpleSubstitutionProcessSequenceSimulator simulatorS(process1);

  unsigned int nsim = 1000;
  unsigned int seqlen = 2000;
  unsigned int count05 = 0;
  unsigned int count01 = 0;

  for (unsigned int i = 0; i < nsim; ++i)
  {
    ApplicationTools::displayGauge(i, nsim - 1);
    double pvalue = testBowker(simulatorS, seqlen);
    if (pvalue < 0.05)
      count05++;
    if (pvalue < 0.01)
      count01++;
  }

  double p05 = (static_cast<double>(count05) / static_cast<double>(nsim));
  double p01 = (static_cast<double>(count01) / static_cast<double>(nsim));
  cout << "\n" << p05 << "\t" << p01 << endl;

  if (abs(p05 - 0.05) > 0.05)
    return 1;
  if (abs(p01 - 0.01) > 0.01)
    return 1;


  ///////////////////////////////////////////////////////
  // Then test homogeneous, non-stationary model:
  cout << "..:: Testing with homogeneous, non-stationary model ::.." << endl;

  auto model2 = make_shared<T92>(alphabet, 3., 0.65);
  auto rootFreqs2 = make_shared<GCFrequencySet>(alphabet, 0.4);

  shared_ptr<SubstitutionProcessInterface> process2 = NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model2, rdist, phyloTree, rootFreqs2);

  SimpleSubstitutionProcessSequenceSimulator simulatorNS(process2);

  count05 = 0;
  count01 = 0;
  for (unsigned int i = 0; i < nsim; ++i)
  {
    ApplicationTools::displayGauge(i, nsim - 1);
    double pvalue = testBowker(simulatorNS, seqlen);
    if (pvalue < 0.05)
      count05++;
    if (pvalue < 0.01)
      count01++;
  }

  p05 = (static_cast<double>(count05) / static_cast<double>(nsim));
  p01 = (static_cast<double>(count01) / static_cast<double>(nsim));
  cout << "\n" << p05 << "\t" << p01 << endl;
  if (abs(p05 - 0.05) > 0.05)
    return 1;
  if (abs(p01 - 0.01) > 0.01)
    return 1;


  ///////////////////////////////////////////////////////////////
  // Now test non-homogeneous model, with distinct GC content:
  cout << "..:: Testing with non-homogeneous, non-stationary model ::.." << endl;

  auto model3 = make_shared<T92>(alphabet, 3., 0.65);
  auto rootFreqs3 = make_shared<GCFrequencySet>(alphabet, 0.65);

  std::vector<string> globalParameterNames = {"T92.kappa"};

  shared_ptr<SubstitutionProcessInterface> process3 = NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(model3, rdist, phyloTree, rootFreqs3, globalParameterNames);

  process3->setParameterValue("T92.theta_1", 0.3);
  process3->setParameterValue("T92.theta_2", 0.8);

  SimpleSubstitutionProcessSequenceSimulator simulatorNHGC(process3);

  count05 = 0;
  count01 = 0;
  for (unsigned int i = 0; i < nsim; ++i)
  {
    ApplicationTools::displayGauge(i, nsim - 1);
    double pvalue = testBowker(simulatorNHGC, seqlen);
    if (pvalue < 0.05)
      count05++;
    if (pvalue < 0.01)
      count01++;
  }
  p05 = (static_cast<double>(count05) / static_cast<double>(nsim));
  p01 = (static_cast<double>(count01) / static_cast<double>(nsim));

  cout << "\n" << p05 << "\t" << p01 << endl;
  if (p05 < 0.7)
    return 1;
  if (p01 < 0.7)
    return 1;


  ////////////////////////////////////////////////////////
  // Now test non-homogeneous model, with distinct ts/tv:
  cout << "..:: Testing with non-homogeneous, stationary model ::.." << endl;

  auto model4 = make_shared<T92>(alphabet, 3., 0.5);
  auto rootFreqs4 = make_shared<GCFrequencySet>(alphabet, 0.5);

  std::vector<string> globalParameterNames2 = {"T92.theta"};

  shared_ptr<SubstitutionProcessInterface> process4 = NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(model4, rdist, phyloTree, rootFreqs4, globalParameterNames2);

  process4->setParameterValue("T92.kappa_1", 2);
  process4->setParameterValue("T92.kappa_2", 7);

  SimpleSubstitutionProcessSequenceSimulator simulatorNHTsTv(process4);

  count05 = 0;
  count01 = 0;
  for (unsigned int i = 0; i < nsim; ++i)
  {
    ApplicationTools::displayGauge(i, nsim - 1);
    double pvalue = testBowker(simulatorNHTsTv, seqlen);
    if (pvalue < 0.05)
      count05++;
    if (pvalue < 0.01)
      count01++;
  }

  p05 = (static_cast<double>(count05) / static_cast<double>(nsim));
  p01 = (static_cast<double>(count01) / static_cast<double>(nsim));
  cout << "\n" << p05 << "\t" << p01 << endl;
  if (abs(p05 - 0.05) > 0.05)
    return 1;
  if (abs(p01 - 0.01) > 0.01)
    return 1;

  // -------------

  return 0;
}
