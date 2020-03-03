//
// File: test_simulations.cpp
// Created by: Julien Dutheil
// Created on: Fri Jan 21 17:21 2011
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

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/FrequenciesSet/NucleotideFrequenciesSet.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Simulation/SubstitutionProcessSequenceSimulator.h>
#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/NewPhyl/SingleProcessPhyloLikelihood_DF.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {

  Newick reader;
  auto phyloTree = std::unique_ptr<bpp::PhyloTree>(reader.parenthesisToPhyloTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);", false, "", false, false));
  auto partree = new bpp::ParametrizablePhyloTree(*phyloTree);

  vector<string> seqNames= phyloTree->getAllLeavesNames();
  //-------------

  NucleicAlphabet* alphabet = new DNA();
  auto model = std::make_shared<T92>(alphabet, 3.);
  DiscreteDistribution* rdist = new ConstantRateDistribution();
  FrequenciesSet* rootFreqs = new GCFrequenciesSet(alphabet);
  std::vector<std::string> globalParameterNames({"T92.kappa"});

  auto process=NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(model, rdist, partree, rootFreqs, globalParameterNames);

  vector<double> thetas;
  for (unsigned int i = 0; i < process->getNumberOfModels(); ++i) {
    double theta = RandomTools::giveRandomNumberBetweenZeroAndEntry(0.99) + 0.005;
    cout << "Theta" << i << " set to " << theta << endl; 
    process->setParameterValue("T92.theta_" + TextTools::toString(i + 1), theta);
    thetas.push_back(theta);
  }

  SimpleSubstitutionProcessSequenceSimulator simulator(*process);

  unsigned int n = 100000;
  OutputStream* profiler  = new StlOutputStream(new ofstream("profile.txt", ios::out));
  OutputStream* messenger = new StlOutputStream(new ofstream("messages.txt", ios::out));

  //Check fast simulation first:
 
  cout << "Fast check:" << endl;
 
  //Generate data set:
  VectorSiteContainer sites(seqNames, alphabet);
  for (unsigned int i = 0; i < n; ++i) {
    unique_ptr<Site> site(simulator.simulateSite());
    site->setPosition(static_cast<int>(i));
    sites.addSite(*site, false);
  }

  cout << "fit model" << endl;
  
  //Now fit model:
  bpp::dataflow::Context context;
  auto l = std::make_shared<bpp::dataflow::LikelihoodCalculationSingleProcess>(context, sites, *process);
  bpp::dataflow::SingleProcessPhyloLikelihood_DF llh(context, l, l->getParameters());

  OptimizationTools::optimizeNumericalParameters2(
      llh, llh.getParameters(), 0,
      0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

  //Now compare estimated values to real ones:
  for (size_t i = 0; i < thetas.size(); ++i) {
    cout << thetas[i] << "\t" << process->getModel(i)->getParameter("theta").getValue() << endl;
    double diff = abs(thetas[i] - process->getModel(i)->getParameter("theta").getValue());
    if (diff > 0.1)
      return 1;
  }

  //Now try detailed simulations:

  cout << "Detailed check:" << endl;
  
  //Generate data set:
  VectorSiteContainer sites2(seqNames, alphabet);
  for (unsigned int i = 0; i < n; ++i) {
    auto result = simulator.dSimulateSite();
    unique_ptr<Site> site(result->getSite(dynamic_cast<const TransitionModel&>(*simulator.getSubstitutionProcess()->getModel(0))));
    site->setPosition(static_cast<int>(i));
    sites2.addSite(*site, false);
    delete result;
  }

  //Now fit model:
  auto process2 = std::shared_ptr<SubstitutionProcess>(process->clone());
  auto l2 = std::make_shared<bpp::dataflow::LikelihoodCalculationSingleProcess>(context, sites, *process2);
  bpp::dataflow::SingleProcessPhyloLikelihood_DF llh2(context, l2, l2->getParameters());

  OptimizationTools::optimizeNumericalParameters2(
    llh2, llh2.getParameters(), 0,
    0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

  //Now compare estimated values to real ones:
  for (size_t i = 0; i < thetas.size(); ++i) {
    cout << thetas[i] << "\t" << process2->getModel(i)->getParameter("theta").getValue() << endl;
    double diff = abs(thetas[i] - process2->getModel(i)->getParameter("theta").getValue());
    if (diff > 0.1)
    {
      cout << "difference too large" << endl;
      return 1;
    }
  }

  //-------------
  delete partree;
  delete alphabet;
  delete rdist;

  return 0;
}
