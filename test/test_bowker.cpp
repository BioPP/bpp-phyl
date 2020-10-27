//
// File: test_bowker.cpp
// Created by: Julien Dutheil
// Created on: Wed Apr 27 14:19 2011
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/SequenceTools.h>

#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Io/Newick.h>

#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Model/FrequencySet/FrequencySet.h>

#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSequenceSimulator.h>

#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>

#include <iostream>

using namespace bpp;
using namespace std;

double testBowker(const SimpleSubstitutionProcessSequenceSimulator& sim, size_t seqlen) {
  auto sites(sim.simulate(seqlen));
  unique_ptr<BowkerTest> bTest(SequenceTools::bowkerTest(sites->getSequence(0), sites->getSequence(1)));
  return bTest->getPValue();
}

int main() {

  Newick reader;
  auto phyloTree = std::unique_ptr<PhyloTree>(reader.parenthesisToPhyloTree("(A:0.02, B:0.02);", false, "", false, false));
  auto paramPhyloTree = new ParametrizablePhyloTree(*phyloTree);

  //First test stationnary model:
  cout << "..:: Testing with stationary model ::.." << endl;
  unique_ptr<NucleicAlphabet> alphabet(new DNA());
  std::shared_ptr<T92> model1(new T92(alphabet.get(), 3., 0.65));
  shared_ptr<DiscreteDistribution> rdist(new ConstantRateDistribution());
  
  auto process1 = NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model1, rdist.get(), paramPhyloTree);

  SimpleSubstitutionProcessSequenceSimulator simulatorS(*process1);

  unsigned int nsim = 1000;
  unsigned int seqlen = 2000;
  unsigned int count05 = 0;
  unsigned int count01 = 0;

  for (unsigned int i = 0; i < nsim; ++i) {
    ApplicationTools::displayGauge(i, nsim - 1);
    double pvalue = testBowker(simulatorS, seqlen);
    if (pvalue < 0.05) count05++;
    if (pvalue < 0.01) count01++;
  }

  double p05 = (static_cast<double>(count05) / static_cast<double>(nsim));
  double p01 = (static_cast<double>(count01) / static_cast<double>(nsim));
  cout << "\n" << p05 << "\t" << p01 << endl;

  if (abs(p05 - 0.05) > 0.05) return 1;
  if (abs(p01 - 0.01) > 0.01) return 1;

  
  ///////////////////////////////////////////////////////
  //Then test homogeneous, non-stationary model:
  cout << "..:: Testing with homogeneous, non-stationary model ::.." << endl;
  
  std::shared_ptr<T92> model2(new T92(alphabet.get(), 3., 0.65));
  auto rootFreqs2 = new GCFrequencySet(alphabet.get(), 0.4);

  auto process2 = NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model2, rdist.get(), paramPhyloTree, rootFreqs2);

  SimpleSubstitutionProcessSequenceSimulator simulatorNS(*process2);

  count05 = 0;
  count01 = 0;
  for (unsigned int i = 0; i < nsim; ++i) {
    ApplicationTools::displayGauge(i, nsim - 1);
    double pvalue = testBowker(simulatorNS, seqlen);
    if (pvalue < 0.05) count05++;
    if (pvalue < 0.01) count01++;
  }
  
  p05 = (static_cast<double>(count05) / static_cast<double>(nsim));
  p01 = (static_cast<double>(count01) / static_cast<double>(nsim)); 
  cout << "\n" << p05 << "\t" << p01 << endl;
  if (abs(p05 - 0.05) > 0.05) return 1;
  if (abs(p01 - 0.01) > 0.01) return 1;

  
  ///////////////////////////////////////////////////////////////
  //Now test non-homogeneous model, with distinct GC content:
  cout << "..:: Testing with non-homogeneous, non-stationary model ::.." << endl;

  std::shared_ptr<T92> model3(new T92(alphabet.get(), 3., 0.65));
  auto rootFreqs3 = new GCFrequencySet(alphabet.get(), 0.65);

  std::vector<string> globalParameterNames={"T92.kappa"};
  
  auto process3 = NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(model3, rdist.get(), paramPhyloTree, rootFreqs3, globalParameterNames);
  
  process3->setParameterValue("T92.theta_1", 0.3);
  process3->setParameterValue("T92.theta_2", 0.8);

  SimpleSubstitutionProcessSequenceSimulator simulatorNHGC(*process3);

  count05 = 0;
  count01 = 0;
  for (unsigned int i = 0; i < nsim; ++i) {
    ApplicationTools::displayGauge(i, nsim - 1);
    double pvalue = testBowker(simulatorNHGC, seqlen);
    if (pvalue < 0.05) count05++;
    if (pvalue < 0.01) count01++;
  }
  p05 = (static_cast<double>(count05) / static_cast<double>(nsim));
  p01 = (static_cast<double>(count01) / static_cast<double>(nsim));

  cout << "\n" << p05 << "\t" << p01 << endl;
  if (p05 < 0.7) return 1;
  if (p01 < 0.7) return 1;

  
  ////////////////////////////////////////////////////////
  //Now test non-homogeneous model, with distinct ts/tv:
  cout << "..:: Testing with non-homogeneous, stationary model ::.." << endl;

  std::shared_ptr<T92> model4(new T92(alphabet.get(), 3., 0.5));
  auto rootFreqs4 = new GCFrequencySet(alphabet.get(), 0.5);
  
  std::vector<string> globalParameterNames2={"T92.theta"};
  
  auto process4 = NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(model4, rdist.get(), paramPhyloTree, rootFreqs4, globalParameterNames2);
  
  process4->setParameterValue("T92.kappa_1", 2);
  process4->setParameterValue("T92.kappa_2", 7);
  
  SimpleSubstitutionProcessSequenceSimulator simulatorNHTsTv(*process4);

  count05 = 0;
  count01 = 0;
  for (unsigned int i = 0; i < nsim; ++i) {
    ApplicationTools::displayGauge(i, nsim - 1);
    double pvalue = testBowker(simulatorNHTsTv, seqlen);
    if (pvalue < 0.05) count05++;
    if (pvalue < 0.01) count01++;
  }

  p05 = (static_cast<double>(count05) / static_cast<double>(nsim));
  p01 = (static_cast<double>(count01) / static_cast<double>(nsim));
  cout << "\n" << p05 << "\t" << p01 << endl;
  if (abs(p05 - 0.05) > 0.05) return 1;
  if (abs(p01 - 0.01) > 0.01) return 1;

  //-------------

  return 0;
}
