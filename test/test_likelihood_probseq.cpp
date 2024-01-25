//
// File: test_likelihood_probseq.cpp
// Created by: Laurent Guéguen
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

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Io/Pasta.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Legacy/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Legacy/OptimizationTools.h>

#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/SimpleSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>

#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

#include <iostream>


using namespace bpp;
using namespace std;

void fitModelHSR(
    shared_ptr<SubstitutionModelInterface> model,
    shared_ptr<DiscreteDistributionInterface> rdist,
    const Tree& tree, const ParametrizablePhyloTree&  new_tree,
    shared_ptr<const ProbabilisticSiteContainerInterface> sites,
    double initialValue, double finalValue)
{
  RHomogeneousTreeLikelihood tl(tree, *sites, model, rdist, false, false, false);
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
  Parameter p2("T92.kappa",3);
  
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
  
  ApplicationTools::startTimer();
  
  auto process = make_shared<RateAcrossSitesSubstitutionProcess>(
      shared_ptr<SubstitutionModelInterface>(model->clone()),
      shared_ptr<DiscreteDistributionInterface>(rdist->clone()),
      shared_ptr<ParametrizablePhyloTree>(new_tree.clone()));

  Context context;                        
  auto lik = make_shared<LikelihoodCalculationSingleProcess>(context, sites, process);
  auto newTl = make_shared<SingleProcessPhyloLikelihood>(context, lik);

  
  cout << "NewTL:    " << setprecision(20) << newTl->getValue() << endl;
  cout << "NewTL D1: " << setprecision(20) << newTl->getFirstOrderDerivative("BrLen2") << endl;
  cout << "NewTL D2: " << setprecision(20) << newTl->getSecondOrderDerivative("BrLen2") << endl;
  ApplicationTools::displayResult("* initial likelihood", newTl->getValue());
  if (abs(newTl->getValue() - initialValue) > 0.001)
    throw Exception("Incorrect initial value.");
  cout << endl;

  for (size_t i = 0; i < n; ++i) { 
    ApplicationTools::displayGauge(i, n-1);
    newTl->matchParametersValues(pl1);
    newTl->getValue();
    newTl->matchParametersValues(pl2);
    newTl->getValue();
  }

  cout << endl;
  ApplicationTools::displayTime("New Likelihood: model upgrade");

  ApplicationTools::startTimer();
  for (size_t i = 0; i < n; ++i) { 
    ApplicationTools::displayGauge(i, n-1);
    newTl->matchParametersValues(pl3);
    newTl->getValue();
    newTl->matchParametersValues(pl4);
    newTl->getValue();
  }
  cout << endl;
  
  ApplicationTools::displayTime("New Likelihood: brlen upgrade");

  cout << endl;
  
  cout << "==========================================" << endl;
  cout << "==========================================" << endl;
  cout << endl;
  
  cout << "Optimization : " << endl;
  cout << endl;

  uint nboptim = 1000;
  
  auto tlop = make_shared<RHomogeneousTreeLikelihood>(tree, *sites,
      shared_ptr<SubstitutionModelInterface>(model->clone()),
      shared_ptr<DiscreteDistributionInterface>(rdist->clone()),
      false, false);
  tlop->initialize();

  LegacyOptimizationTools::optimizeNumericalParameters2(
      tlop, tlop->getParameters(),
      0, 0.000001, nboptim, 0, 0);
  cout << setprecision(20) << tlop->getValue() << endl;
  ApplicationTools::displayResult("* lnL after full optimization (old)", tlop->getValue());
  if (abs(tlop->getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value.");
  tlop->getParameters().printParameters(cout);


  process.reset(new RateAcrossSitesSubstitutionProcess(model, std::shared_ptr<DiscreteDistributionInterface>(rdist->clone()), std::shared_ptr<ParametrizablePhyloTree>(new_tree.clone())));  
  lik.reset(new LikelihoodCalculationSingleProcess(context, sites, process));
  newTl.reset(new SingleProcessPhyloLikelihood(context, lik));

  ParameterList opln1 = process->getBranchLengthParameters(true);
  
  OptimizationTools::optimizeNumericalParameters2(
      newTl, newTl->getParameters(),
      0, 0.000001, nboptim, 0, 0);
  cout << setprecision(20) << newTl->getValue() << endl;
  ApplicationTools::displayResult("* lnL after full optimization (new)", newTl->getValue());
  if (abs(newTl->getValue() - finalValue) > 0.001)
    throw Exception("Incorrect final value.");
  newTl->getParameters().printParameters(cout);
}


int main() {
  unique_ptr<TreeTemplate<Node> > tree(TreeTemplateTools::parenthesisToTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);"));
  vector<string> seqNames= tree->getLeavesNames();
  vector<int> ids = tree->getNodesId();

  Newick reader;
  unique_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);", false, "", false, false));
  
  //-------------

  shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET;
  shared_ptr<const NucleicAlphabet> nucAlphabet = AlphabetTools::DNA_ALPHABET;

  Pasta pasta;
  
  auto sites = make_shared<ProbabilisticVectorSiteContainer>(alphabet);
  pasta.readAlignment("exemple1.pa", *sites);
  
  auto model = make_shared<T92>(nucAlphabet, 3.);
  auto rdist = make_shared<ConstantRateDistribution>();
  try {
    cout << "Testing Single Tree Traversal likelihood class..." << endl;
    fitModelHSR(model, rdist, *tree, *pTree, sites, 222.26297478, 215.7976882);
  } catch (Exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }

  return 0;
}


