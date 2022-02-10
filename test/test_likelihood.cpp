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

#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Legacy/Likelihood/RHomogeneousTreeLikelihood.h>

#include <Bpp/Phyl/Legacy/OptimizationTools.h>

#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/SimpleSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>

#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

#include <iostream>


using namespace bpp;
using namespace std;

void fitModelHSR(std::shared_ptr<SubstitutionModel> model, std::shared_ptr<DiscreteDistribution> rdist,
                 const Tree& tree,
                 const ParametrizablePhyloTree& partree,
                 const SiteContainer& sites,
                 double initialValue, double finalValue)
{
  ApplicationTools::startTimer();

  RHomogeneousTreeLikelihood tl(tree, sites, model->clone(), rdist->clone(), false, false);
  tl.initialize();

  ApplicationTools::displayResult("Test model", model->getName());
  cout << "OldTL: " << setprecision(20) << tl.getValue() << endl;
  cout << "OldTL D1: " << setprecision(20) << tl.getFirstOrderDerivative("BrLen2") << endl;
  cout << "OldTL D2: " << setprecision(20) << tl.getSecondOrderDerivative("BrLen2") << endl;
  ApplicationTools::displayResult("* initial likelihood", tl.getValue());
  if (abs(tl.getValue() - initialValue) > 0.001)
    throw Exception("Incorrect initial value.");
  cout << endl;

  Parameter p1("T92.kappa",0.1);
  Parameter p2("T92.kappa",0.2);

  ParameterList pl1;pl1.addParameter(p1);
  ParameterList pl2;pl2.addParameter(p2);

  Parameter p3("BrLen1",0.1);
  Parameter p4("BrLen1",0.2);

  ParameterList pl3;pl3.addParameter(p3);
  ParameterList pl4;pl4.addParameter(p4);

  unsigned int n = 100000;

  ApplicationTools::startTimer();
  for (size_t i = 0; i < n; ++i) { 
    ApplicationTools::displayGauge(i, n-1);
    tl.matchParametersValues(pl1);
    tl.getValue();
    tl.matchParametersValues(pl2);
    tl.getValue();
  }
  cout << endl;
  ApplicationTools::displayTime("Old Likelihood: model upgrade");

  ApplicationTools::startTimer();
  for (size_t i = 0; i < n; ++i) { 
    ApplicationTools::displayGauge(i, n-1);
    tl.matchParametersValues(pl3);
    tl.getValue();
    tl.matchParametersValues(pl4);
    tl.getValue();
  }
  cout << endl;

  ApplicationTools::displayTime("Old Likelihood: brlen upgrade");

  cout << "=============================" << endl;

  auto process = std::make_shared<RateAcrossSitesSubstitutionProcess>(std::shared_ptr<SubstitutionModel>(model->clone()), std::shared_ptr<DiscreteDistribution>(rdist->clone()), std::shared_ptr<ParametrizablePhyloTree>(partree.clone()));

  Context context;                        
  auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, sites, *process);
  SingleProcessPhyloLikelihood llh(context, lik);

  llh.getFirstOrderDerivative("BrLen0");
  llh.getFirstOrderDerivative("BrLen1");
  llh.getFirstOrderDerivative("BrLen2");
  llh.getFirstOrderDerivative("BrLen3");
  llh.getFirstOrderDerivative("BrLen4");
  
  cout << "NewTL: " << setprecision(20) << llh.getValue() << endl;
  cout << "NewTL D1: " << setprecision(20) << llh.getFirstOrderDerivative("BrLen2") << endl;
  cout << "NewTL D2: " << setprecision(20) << llh.getSecondOrderDerivative("BrLen2") << endl;
  
  ApplicationTools::displayResult("* initial likelihood", llh.getValue());
  if (abs(llh.getValue() - initialValue) > 0.001)
    throw Exception("Incorrect initial value.");
  cout << endl;

  for (size_t i = 0; i < n; ++i) { 
    ApplicationTools::displayGauge(i, n-1);
    llh.matchParametersValues(pl1);
    llh.getValue();
    llh.matchParametersValues(pl2);
    llh.getValue();
  }

  cout << endl;
  ApplicationTools::displayTime("New Likelihood: model upgrade");

  ApplicationTools::startTimer();
  for (size_t i = 0; i < n; ++i) { 
    ApplicationTools::displayGauge(i, n-1);
    llh.matchParametersValues(pl3);
    llh.getValue();
    llh.matchParametersValues(pl4);
    llh.getValue();
  }
  cout << endl;
  
  ApplicationTools::displayTime("New Likelihood: brlen upgrade");

  cout << endl;
  
  cout << "==========================================" << endl;
  cout << "==========================================" << endl;
  cout << endl;
  
  cout << "Optimization : " << endl;
  cout << endl;

  uint nboptim=1000;

  ///////////////
  
  RHomogeneousTreeLikelihood tlop(tree, sites, model->clone(), rdist->clone(), false, false);
  tlop.initialize();

  OptimizationToolsOld::optimizeNumericalParameters2(&tlop, tlop.getParameters(), 0, 0.000001, nboptim, 0, 0);
  cout << setprecision(20) << tlop.getValue() << endl;
  ApplicationTools::displayResult("* lnL after full optimization (old)", tlop.getValue());
  if (abs(tlop.getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value.");
  tlop.getParameters().printParameters(cout);

  
  process = std::make_shared<RateAcrossSitesSubstitutionProcess>(model, std::shared_ptr<DiscreteDistribution>(rdist->clone()), std::shared_ptr<ParametrizablePhyloTree>(partree.clone()));

  Context context2;
  
  lik.reset(new LikelihoodCalculationSingleProcess(context2, sites, *process));
  
  SingleProcessPhyloLikelihood llh2(context2, lik);

  ParameterList opln1=process->getBranchLengthParameters(true);

  OptimizationTools::optimizeNumericalParameters2(llh2, llh2.getParameters(), 0, 0.000001, nboptim, 0, 0);
  cout << setprecision(20) << llh2.getValue() << endl;
  ApplicationTools::displayResult("* lnL after full optimization (new)", llh2.getValue());
  if (abs(llh2.getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value.");
  llh2.getParameters().printParameters(cout);
}


int main() {
  unique_ptr<TreeTemplate<Node> > tree(TreeTemplateTools::parenthesisToTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);"));
  vector<string> seqNames= tree->getLeavesNames();
  vector<int> ids = tree->getNodesId();

  Newick reader;
  auto pTree = unique_ptr<PhyloTree>(reader.parenthesisToPhyloTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);", false, "", false, false));
  ParametrizablePhyloTree paramphyloTree(*pTree);

  //-------------

  const NucleicAlphabet* alphabet = &AlphabetTools::DNA_ALPHABET;

  VectorSiteContainer sites(alphabet);
  sites.addSequence(BasicSequence("A", "ATCCAGACATGCCGGGACTTTGCAGAGAAGGAGTTGTTTCCCATTGCAGCCCAGGTGGATAAGGAACAGC", alphabet));
  sites.addSequence(BasicSequence("B", "CGTCAGACATGCCGTGACTTTGCCGAGAAGGAGTTGGTCCCCATTGCGGCCCAGCTGGACAGGGAGCATC", alphabet));
  sites.addSequence(BasicSequence("C", "GGTCAGACATGCCGGGAATTTGCTGAAAAGGAGCTGGTTCCCATTGCAGCCCAGGTAGACAAGGAGCATC", alphabet));
  sites.addSequence(BasicSequence("D", "TTCCAGACATGCCGGGACTTTACCGAGAAGGAGTTGTTTTCCATTGCAGCCCAGGTGGATAAGGAACATC", alphabet));

  shared_ptr<SubstitutionModel> model(new T92(alphabet, 3.));
  std::shared_ptr<DiscreteDistribution> rdist(new GammaDiscreteRateDistribution(4, 1.0));
  try {
    cout << "Testing Single Tree Traversal likelihood class..." << endl;
    fitModelHSR(model, rdist, *tree, paramphyloTree, sites, 228.6333642493463, 198.47216106233);
  } catch (Exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }  

  return 0;
}


