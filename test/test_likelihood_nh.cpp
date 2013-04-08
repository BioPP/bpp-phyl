//
// File: test_likelihood_nh.cpp
// Created by: Julien Dutheil
// Created on: Thu Jul 14 11:04 2011
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

#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/FrequenciesSet/NucleotideFrequenciesSet.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Simulation/NonHomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <iostream>

using namespace bpp;
using namespace std;

void fitModelNH(SubstitutionModelSet* model, DiscreteDistribution* rdist, const Tree& tree, const SiteContainer& sites,
    double initialValue, double finalValue, bool reparam) {
  DRNonHomogeneousTreeLikelihood tl(tree, sites, model, rdist, false, reparam);
  tl.initialize();
  cout << setprecision(20) << tl.getValue() << endl;
  ApplicationTools::displayResult("* initial likelihood", tl.getValue());
  if (abs(tl.getValue() - initialValue) > 0.0001)
    throw Exception("Incorrect initial value.");
  OptimizationTools::optimizeTreeScale(&tl);
  ApplicationTools::displayResult("* likelihood after tree scale", tl.getValue());
  OptimizationTools::optimizeNumericalParameters2(&tl, tl.getParameters(), 0, 0.000001, 10000, 0, 0);
  cout << setprecision(20) << tl.getValue() << endl;
  ApplicationTools::displayResult("* likelihood after full optimization", tl.getValue());
  if (abs(tl.getValue() - finalValue) > 0.0001)
    throw Exception("Incorrect final value.");
}

int main() {
  TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree("(((A:0.1, B:0.2):0.3,C:0.1):0.2,(D:0.3,(E:0.2,F:0.05):0.1):0.1);");
  vector<string> seqNames= tree->getLeavesNames();
  vector<int> ids = tree->getNodesId();
  //-------------

  const NucleicAlphabet* alphabet = &AlphabetTools::DNA_ALPHABET;
  FrequenciesSet* rootFreqs = new GCFrequenciesSet(alphabet);
  SubstitutionModel* model = new T92(alphabet, 3.);
  std::vector<std::string> globalParameterNames;
  globalParameterNames.push_back("T92.kappa");
  SubstitutionModelSet* modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameterNames);
  //DiscreteDistribution* rdist = new ConstantDistribution(1.0, true);
  //Very difficult to optimize on small datasets:
  DiscreteDistribution* rdist = new GammaDiscreteDistribution(4, 1.0);
  rdist->aliasParameters("alpha", "beta");

  size_t nsites = 1000;
  unsigned int nrep = 20;
  size_t nmodels = modelSet->getNumberOfModels();
  vector<double> thetas(nmodels);
  vector<double> thetasEst1(nmodels);
  vector<double> thetasEst2(nmodels);

  for (size_t i = 0; i < nmodels; ++i) {
    double theta = RandomTools::giveRandomNumberBetweenZeroAndEntry(0.99) + 0.005;
    cout << "Theta" << i << " set to " << theta << endl; 
    modelSet->setParameterValue("T92.theta_" + TextTools::toString(i + 1), theta);
    thetas[i] = theta;
  }
  NonHomogeneousSequenceSimulator simulator(modelSet, rdist, tree);
 
  for (unsigned int j = 0; j < nrep; j++) {

    OutputStream* profiler  = new StlOutputStream(new ofstream("profile.txt", ios::out));
    OutputStream* messenger = new StlOutputStream(new ofstream("messages.txt", ios::out));

    //Simulate data:
    auto_ptr<SiteContainer> sites(simulator.simulate(nsites));
    //Now fit model:
    auto_ptr<SubstitutionModelSet> modelSet2(modelSet->clone());
    auto_ptr<SubstitutionModelSet> modelSet3(modelSet->clone());
    RNonHomogeneousTreeLikelihood tl(*tree, *sites.get(), modelSet2.get(), rdist, true, true, false);
    tl.initialize();
    RNonHomogeneousTreeLikelihood tl2(*tree, *sites.get(), modelSet3.get(), rdist, true, true, true);
    tl2.initialize();
   
    unsigned int c1 = OptimizationTools::optimizeNumericalParameters2(
        &tl, tl.getParameters(), 0,
        0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

    unsigned int c2 = OptimizationTools::optimizeNumericalParameters2(
        &tl2, tl2.getParameters(), 0,
        0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

    cout << c1 << ": " << tl.getValue() << "\t" << c2 << ": " << tl2.getValue() << endl;
      
    for (size_t i = 0; i < nmodels; ++i) {
      cout << modelSet2->getModel(i)->getParameter("theta").getValue() << "\t" << modelSet3->getModel(i)->getParameter("theta").getValue() << endl;
      //if (abs(modelSet2->getModel(i)->getParameter("theta").getValue() - modelSet3->getModel(i)->getParameter("theta").getValue()) > 0.1)
      //  return 1;
      thetasEst1[i] +=  modelSet2->getModel(i)->getParameter("theta").getValue();
      thetasEst2[i] +=  modelSet3->getModel(i)->getParameter("theta").getValue();
    }
  }
  thetasEst1 /= static_cast<double>(nrep);
  thetasEst2 /= static_cast<double>(nrep);

  //Now compare estimated values to real ones:
  for (size_t i = 0; i < thetas.size(); ++i) {
     cout << thetas[i] << "\t" << thetasEst1[i] << "\t" << thetasEst2[i] << endl;
     double diff1 = abs(thetas[i] - thetasEst1[i]);
     double diff2 = abs(thetas[i] - thetasEst2[i]);
     if (diff1 > 0.2 || diff2 > 0.2)
        return 1;
  }

  //-------------
  delete tree;
  delete modelSet;
  delete rdist;

  return 0;
}
