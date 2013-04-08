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
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/FrequenciesSet/NucleotideFrequenciesSet.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);");
  vector<string> seqNames= tree->getLeavesNames();
  vector<int> ids = tree->getNodesId();
  //-------------

  NucleicAlphabet* alphabet = new DNA();
  SubstitutionModel* model = new T92(alphabet, 3.);
  FrequenciesSet* rootFreqs = new GCFrequenciesSet(alphabet);
  std::vector<std::string> globalParameterNames;
  globalParameterNames.push_back("T92.kappa");
  SubstitutionModelSet* modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameterNames);
  DiscreteDistribution* rdist = new ConstantRateDistribution();
  vector<double> thetas;
  for (unsigned int i = 0; i < modelSet->getNumberOfModels(); ++i) {
    double theta = RandomTools::giveRandomNumberBetweenZeroAndEntry(0.99) + 0.005;
    cout << "Theta" << i << " set to " << theta << endl; 
    modelSet->setParameterValue("T92.theta_" + TextTools::toString(i + 1), theta);
    thetas.push_back(theta);
  }
  NonHomogeneousSequenceSimulator simulator(modelSet, rdist, tree);

  unsigned int n = 100000;
  OutputStream* profiler  = new StlOutputStream(new ofstream("profile.txt", ios::out));
  OutputStream* messenger = new StlOutputStream(new ofstream("messages.txt", ios::out));

  //Check fast simulation first:
  
  //Generate data set:
  VectorSiteContainer sites(seqNames, alphabet);
  for (unsigned int i = 0; i < n; ++i) {
    auto_ptr<Site> site(simulator.simulate());
    site->setPosition(i);
    sites.addSite(*site, false);
  }

  //Now fit model:
  SubstitutionModelSet* modelSet2 = modelSet->clone();
  RNonHomogeneousTreeLikelihood tl(*tree, sites, modelSet2, rdist);
  tl.initialize();

  OptimizationTools::optimizeNumericalParameters2(
      &tl, tl.getParameters(), 0,
      0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

  //Now compare estimated values to real ones:
  for (size_t i = 0; i < thetas.size(); ++i) {
    cout << thetas[i] << "\t" << modelSet2->getModel(i)->getParameter("theta").getValue() << endl;
    double diff = abs(thetas[i] - modelSet2->getModel(i)->getParameter("theta").getValue());
    if (diff > 0.1)
      return 1;
  }
  delete modelSet2;

  //Now try detailed simulations:

  //Generate data set:
  VectorSiteContainer sites2(seqNames, alphabet);
  for (unsigned int i = 0; i < n; ++i) {
    RASiteSimulationResult* result = simulator.dSimulate();
    auto_ptr<Site> site(result->getSite());
    site->setPosition(i);
    sites2.addSite(*site, false);
    delete result;
  }

  //Now fit model:
  SubstitutionModelSet* modelSet3 = modelSet->clone();
  RNonHomogeneousTreeLikelihood tl2(*tree, sites2, modelSet3, rdist);
  tl2.initialize();

  OptimizationTools::optimizeNumericalParameters2(
      &tl2, tl2.getParameters(), 0,
      0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

  //Now compare estimated values to real ones:
  for (size_t i = 0; i < thetas.size(); ++i) {
    cout << thetas[i] << "\t" << modelSet3->getModel(i)->getParameter("theta").getValue() << endl;
    double diff = abs(thetas[i] - modelSet3->getModel(i)->getParameter("theta").getValue());
    if (diff > 0.1)
      return 1;
  }
  delete modelSet3;

  //-------------
  delete tree;
  delete alphabet;
  delete modelSet;
  delete rdist;

  return 0;
}
