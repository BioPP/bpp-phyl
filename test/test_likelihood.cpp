//
// File: test_likelihood.cpp
// Created by: Julien Dutheil
// Created on: Mon APr 04 10:18 2011
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

#include <Bpp/Numeric/Prob.all>
#include <Bpp/Numeric/Matrix.all>
#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Model.all>
#include <Bpp/Phyl/Simulation.all>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/OptimizationTools.h>
#include <iostream>

using namespace bpp;
using namespace std;

void fitModelH(SubstitutionModel* model, DiscreteDistribution* rdist, const Tree& tree, const SiteContainer& sites,
    double initialValue, double finalValue) {
  RHomogeneousTreeLikelihood tl(tree, sites, model, rdist);
  tl.initialize();
  ApplicationTools::displayResult("Test model", model->getName());
  //cout << setprecision(20) << tl.getValue() << endl;
  ApplicationTools::displayResult("* initial likelihood", tl.getValue());
  if (abs(tl.getValue() - initialValue) > 0.00001)
    throw Exception("Incorrect initial value.");
  OptimizationTools::optimizeTreeScale(&tl);
  ApplicationTools::displayResult("* likelihood after tree scale", tl.getValue());
  OptimizationTools::optimizeNumericalParameters2(&tl, tl.getParameters(), 0, 0.000001, 10000, 0, 0);
  //cout << setprecision(20) << tl.getValue() << endl;
  ApplicationTools::displayResult("* likelihood after full optimization", tl.getValue());
  if (abs(tl.getValue() - finalValue) > 0.00001)
    throw Exception("Incorrect final value.");
}

int main() {
  TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);");
  vector<string> seqNames= tree->getLeavesNames();
  vector<int> ids = tree->getNodesId();
  //-------------

  const NucleicAlphabet* alphabet = &AlphabetTools::DNA_ALPHABET;
  SubstitutionModel* model = new T92(alphabet, 3.);
  DiscreteDistribution* rdist = new GammaDiscreteDistribution(1.0, 4);
  rdist->aliasParameters("alpha", "beta");

  VectorSiteContainer sites(alphabet);
  sites.addSequence(BasicSequence("A", "AAATGGCTGTGCACGTC", alphabet));
  sites.addSequence(BasicSequence("B", "GACTGGATCTGCACGTC", alphabet));
  sites.addSequence(BasicSequence("C", "CTCTGGATGTGCACGTG", alphabet));
  sites.addSequence(BasicSequence("D", "AAATGGCGGTGCGCCTA", alphabet));

  try {
    fitModelH(model, rdist, *tree, sites, 75.031104151696752069, 65.03473753351640596065);
  } catch (Exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }
  /*
  FrequenciesSet* rootFreqs = new GCFrequenciesSet(alphabet);
  std::vector<std::string> globalParameterNames;
  globalParameterNames.push_back("T92.kappa");
  SubstitutionModelSet* modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameterNames);
  DiscreteDistribution* rdist = new ConstantDistribution(1.0, true);
  vector<double> thetas;
  for (unsigned int i = 0; i < modelSet->getNumberOfModels(); ++i) {
    double theta = RandomTools::giveRandomNumberBetweenZeroAndEntry(0.99) + 0.005;
    cout << "Theta" << i << " set to " << theta << endl; 
    modelSet->setParameterValue("T92.theta_" + TextTools::toString(i + 1), theta);
    thetas.push_back(theta);
  }
  NonHomogeneousSequenceSimulator simulator(modelSet, rdist, tree);

  unsigned int n = 100000;
  OutputStream* profiler  = new StlOutputStream(auto_ptr<ostream>(new ofstream("profile.txt", ios::out)));
  OutputStream* messenger = new StlOutputStream(auto_ptr<ostream>(new ofstream("messages.txt", ios::out)));

  //Now fit model:
  SubstitutionModelSet* modelSet2 = modelSet->clone();
  RNonHomogeneousTreeLikelihood tl(*tree, sites, modelSet2, rdist);
  tl.initialize();

  OptimizationTools::optimizeNumericalParameters2(
      &tl, tl.getParameters(), 0,
      0.0001, 10000, messenger, profiler, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

  //Now compare estimated values to real ones:
  for (size_t i = 0; i < thetas.size(); ++i) {
    cout << thetas[i] << "\t" << modelSet2->getModel(i)->getParameter("theta").getValue() << endl;
    double diff = abs(thetas[i] - modelSet2->getModel(i)->getParameter("theta").getValue());
    if (diff > 0.1)
      return 1;
  }
  delete modelSet2;

  */

  //-------------
  delete tree;
  delete model;
  delete rdist;

  return 0;
}
