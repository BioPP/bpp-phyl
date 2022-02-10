//
// File: test_likelihood_clock.cpp
// Created by: Julien Dutheil
// Created on: Mon Jul 12 14:57 2011
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
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <iostream>

#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

using namespace bpp;
using namespace std;

void fitModelH(std::shared_ptr<SubstitutionModel> model, std::shared_ptr<DiscreteDistribution> rdist,
               std::shared_ptr<PhyloTree> tree, const VectorSiteContainer& sites,
               double initialValue, double finalValue)
{
  RateAcrossSitesSubstitutionProcess process(model, rdist, tree);

  Context context;
  
  auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, sites, process);
  
  SingleProcessPhyloLikelihood llh(context, lik);

  ApplicationTools::displayResult("Test model", model->getName());
  double initValue=llh.getValue();

  cout << setprecision(20) << initValue << endl;  
  ApplicationTools::displayResult("* initial likelihood", initValue);
  if (abs(initValue - initialValue) > 0.001)
    throw Exception("Incorrect initial value:" + TextTools::toString(initValue) + "<>" + TextTools::toString(initialValue));
  unique_ptr<OutputStream> messenger(new StlOutputStream(new ofstream("messages.txt", ios::out)));
  unique_ptr<OutputStream> profiler(new StlOutputStream(new ofstream("profile.txt", ios::out)));
  profiler->setPrecision(20);

  
  OptimizationTools::optimizeNumericalParameters2(llh, llh.getParameters(), 0, 0.000001, 10000, messenger.get(), profiler.get(), false, true, 2, OptimizationTools::OPTIMIZATION_NEWTON);
  cout << setprecision(20) << llh.getValue() << endl;
  ApplicationTools::displayResult("* likelihood after full optimization", llh.getValue());
  llh.getParameters().printParameters(cout);
  if (abs(llh.getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value:" + TextTools::toString(llh.getValue()) + "<>" + TextTools::toString(finalValue));
}

void fitModelHClock(std::shared_ptr<SubstitutionModel> model, std::shared_ptr<DiscreteDistribution> rdist,
                    std::shared_ptr<PhyloTree> tree, const VectorSiteContainer& sites,

                    double initialValue, double finalValue)
{
  RateAcrossSitesSubstitutionProcess process(model, rdist, tree);

  Context context;
  
  auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, sites, process);
  lik->setClockLike();
  
  SingleProcessPhyloLikelihood llh(context, lik);
    
  ApplicationTools::displayResult("Test model", model->getName());
  double initValue=llh.getValue();

  cout << setprecision(20) << initValue << endl;  
  ApplicationTools::displayResult("* initial likelihood", initValue);
  if (abs(initValue - initialValue) > 0.001)
    throw Exception("Incorrect initial value:" + TextTools::toString(initValue) + "<>" + TextTools::toString(initialValue));
  unique_ptr<OutputStream> messenger(new StlOutputStream(new ofstream("messages.txt", ios::out)));
  unique_ptr<OutputStream> profiler(new StlOutputStream(new ofstream("profile.txt", ios::out)));
  profiler->setPrecision(20);

  OptimizationTools::optimizeNumericalParameters2(llh, llh.getParameters(), 0, 0.000001, 10000, messenger.get(), profiler.get());
  cout << setprecision(20) << llh.getValue() << endl;
  ApplicationTools::displayResult("* likelihood after full optimization", llh.getValue());
  llh.getParameters().printParameters(cout);
  if (abs(llh.getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value:" + TextTools::toString(llh.getValue()) + "<>" + TextTools::toString(finalValue));
}

int main() {
  bpp::Newick reader;
  auto phyloTree = std::shared_ptr<bpp::PhyloTree>(reader.parenthesisToPhyloTree("(((A:0.01, B:0.01):0.02,C:0.03):0.01,D:0.04);", false, "", false, false));
  
  const NucleicAlphabet* alphabet = &AlphabetTools::DNA_ALPHABET;
  shared_ptr<SubstitutionModel> model(new T92(alphabet, 3.));
  DiscreteDistribution* rdist = new ConstantRateDistribution();

  VectorSiteContainer sites(alphabet);
  sites.addSequence(BasicSequence("A", "AAATGGCTGTGCACGTC", alphabet));
  sites.addSequence(BasicSequence("B", "AACTGGATCTGCATGTC", alphabet));
  sites.addSequence(BasicSequence("C", "ATCTGGACGTGCACGTG", alphabet));
  sites.addSequence(BasicSequence("D", "CAACGGGAGTGCGCCTA", alphabet));

  try {
    fitModelH(std::shared_ptr<SubstitutionModel>(model->clone()), std::shared_ptr<DiscreteDistribution>(rdist->clone()), std::shared_ptr<PhyloTree>(phyloTree->clone()), sites, 94.3957, 71.0564);
  } catch (Exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }

  cout << endl << endl;
  
  try {
    fitModelHClock(model, std::shared_ptr<DiscreteDistribution>(rdist->clone()), std::shared_ptr<PhyloTree>(phyloTree->clone()), sites, 94.395699, 72.7196);
  } catch (Exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }

  //-------------
  delete rdist;

  return 0;
}
