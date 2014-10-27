//
// File: test_likelihood.cpp
// Created by: Julien Dutheil
// Created on: Mon Apr 04 10:18 2011
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
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <iostream>

using namespace bpp;
using namespace std;

void fitModelHSR(SubstitutionModel* model, DiscreteDistribution* rdist, const Tree& tree, const SiteContainer& sites,
    double initialValue, double finalValue) {
  RHomogeneousTreeLikelihood tl(tree, sites, model, rdist);
  tl.initialize();
  ApplicationTools::displayResult("Test model", model->getName());
  cout << setprecision(20) << tl.getValue() << endl;
  ApplicationTools::displayResult("* initial likelihood", tl.getValue());
  if (abs(tl.getValue() - initialValue) > 0.001)
    throw Exception("Incorrect initial value.");
  OptimizationTools::optimizeTreeScale(&tl);
  ApplicationTools::displayResult("* likelihood after tree scale", tl.getValue());
  OptimizationTools::optimizeNumericalParameters2(&tl, tl.getParameters(), 0, 0.000001, 10000, 0, 0);
  cout << setprecision(20) << tl.getValue() << endl;
  ApplicationTools::displayResult("* likelihood after full optimization", tl.getValue());
  if (abs(tl.getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value.");
}

void fitModelHDR(SubstitutionModel* model, DiscreteDistribution* rdist, const Tree& tree, const SiteContainer& sites,
    double initialValue, double finalValue) {
  DRHomogeneousTreeLikelihood tl(tree, sites, model, rdist);
  tl.initialize();
  ApplicationTools::displayResult("Test model", model->getName());
  cout << setprecision(20) << tl.getValue() << endl;
  ApplicationTools::displayResult("* initial likelihood", tl.getValue());
  if (abs(tl.getValue() - initialValue) > 0.001)
    throw Exception("Incorrect initial value.");
  OptimizationTools::optimizeTreeScale(&tl);
  ApplicationTools::displayResult("* likelihood after tree scale", tl.getValue());
  OptimizationTools::optimizeNumericalParameters2(&tl, tl.getParameters(), 0, 0.000001, 10000, 0, 0);
  cout << setprecision(20) << tl.getValue() << endl;
  ApplicationTools::displayResult("* likelihood after full optimization", tl.getValue());
  if (abs(tl.getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value.");
}

int main() {
  unique_ptr<TreeTemplate<Node> > tree(TreeTemplateTools::parenthesisToTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);"));
  vector<string> seqNames= tree->getLeavesNames();
  vector<int> ids = tree->getNodesId();
  //-------------

  const NucleicAlphabet* alphabet = &AlphabetTools::DNA_ALPHABET;

  VectorSiteContainer sites(alphabet);
  sites.addSequence(BasicSequence("A", "AAATGGCTGTGCACGTC", alphabet));
  sites.addSequence(BasicSequence("B", "GACTGGATCTGCACGTC", alphabet));
  sites.addSequence(BasicSequence("C", "CTCTGGATGTGCACGTG", alphabet));
  sites.addSequence(BasicSequence("D", "AAATGGCGGTGCGCCTA", alphabet));

  unique_ptr<SubstitutionModel> model(new T92(alphabet, 3.));
  unique_ptr<DiscreteDistribution> rdist(new GammaDiscreteRateDistribution(4, 1.0));
  try {
    cout << "Testing Single Tree Traversal likelihood class..." << endl;
    fitModelHSR(model.get(), rdist.get(), *tree, sites, 85.030942031997312824, 65.72293577214308868406);
  } catch (Exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }  

  model.reset(new T92(alphabet, 3.));
  rdist.reset(new GammaDiscreteRateDistribution(4, 1.0));
  try {
    cout << "Testing Double Tree Traversal likelihood class..." << endl;
    fitModelHDR(model.get(), rdist.get(), *tree, sites, 85.030942031997312824, 65.72293577214308868406);
  } catch (Exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }  

  //Let's compare the derivatives:
  RHomogeneousTreeLikelihood tlsr(*tree, sites, model.get(), rdist.get());
  tlsr.initialize();
  DRHomogeneousTreeLikelihood tldr(*tree, sites, model.get(), rdist.get());
  tldr.initialize();
  vector<string> params = tlsr.getBranchLengthsParameters().getParameterNames();
  for (vector<string>::iterator it = params.begin(); it != params.end(); ++it) {
    double d1sr = tlsr.getFirstOrderDerivative(*it);
    double d1dr = tldr.getFirstOrderDerivative(*it);
    cout << *it << "\t" << d1sr << "\t" << d1dr << endl;
    if (abs(d1sr - d1dr) > 0.000001) return 1;
  }

  return 0;
}
