//
// File: test_mapping.cpp
// Created by: Julien Dutheil
// Created on: Thu Nov 25 16:02 2010
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

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Prob.all>
#include <Bpp/Numeric/Matrix.all>
#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Model.all>
#include <Bpp/Phyl/Simulation.all>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/Mapping.all>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree("((A:0.001, B:0.002):0.008,C:0.01,D:0.1);");
  vector<int> ids = tree->getNodesId();
  ids.pop_back(); //Ignore root

  //-------------

  NucleicAlphabet* alphabet = new DNA();
  SubstitutionModel* model = new GTR(alphabet, 1, 0.2, 0.3, 0.4, 0.4, 0.1, 0.35, 0.35, 0.2);
  //DiscreteDistribution* rdist = new GammaDiscreteDistribution(4, 0.4, 0.4);
  DiscreteDistribution* rdist = new ConstantDistribution(1.0);
  HomogeneousSequenceSimulator simulator(model, rdist, tree);

  unsigned int n = 20000;
  vector< vector<double> > realMap(n);
  VectorSiteContainer sites(tree->getLeavesNames(), alphabet);
  for (unsigned int i = 0; i < n; ++i) {
    ApplicationTools::displayGauge(i, n-1, '=');
    RASiteSimulationResult* result = simulator.dSimulate();
    realMap[i].resize(ids.size());
    for (size_t j = 0; j < ids.size(); ++j) {
      realMap[i][j] = result->getSubstitutionCount(ids[j]);
    }
    auto_ptr<Site> site(result->getSite());
    site->setPosition(i);
    sites.addSite(*site);
    delete result;
  }
  ApplicationTools::displayTaskDone();
  
  //-------------
  //Now build the substitution vectors with the true model:
  //Fasta fasta;
  //fasta.write("Simulations.fasta", sites);
  DRHomogeneousTreeLikelihood drhtl(*tree, sites, model, rdist);
  drhtl.initialize();
  cout << drhtl.getValue() << endl;
  
  SubstitutionCount* sCount = new SimpleSubstitutionCount(alphabet);
  ProbabilisticSubstitutionMapping* probMap = 
    SubstitutionMappingTools::computeSubstitutionVectors(drhtl, *sCount, 0);

  //Check per branch:
  for (unsigned int j = 0; j < ids.size(); ++j) {
    double totalReal = 0;
    double totalObs  = 0;
    for (unsigned int i = 0; i < n; ++i) {
      totalReal += realMap[i][j];
      totalObs  += probMap->getNumberOfSubstitutions(ids[j], i);
    }
    if (tree->isLeaf(ids[j])) cout << tree->getNodeName(ids[j]) << "\t";
    cout << tree->getDistanceToFather(ids[j]) << "\t" << totalReal << "\t" << totalObs << endl;
    if (abs(totalReal - totalObs) / totalReal > 0.1) return 1;
  }
  //Check per site:
  for (unsigned int i = 0; i < n; ++i) {
  }

  //-------------
  delete tree;
  delete alphabet;
  delete model;
  delete rdist;
  delete sCount;
  delete probMap;

  //return (abs(obs - 0.001) < 0.001 ? 0 : 1);
  return 0;
}
