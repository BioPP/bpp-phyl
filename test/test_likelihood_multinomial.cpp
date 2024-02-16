//
// File: test_likelihood_multinomial.cpp
// Created by: Laurent Guéguen
// Created on: jeudi 3 février 2022, à 18h 43
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
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/MultinomialFromTransitionModel.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/OptimizationTools.h>

#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>

#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

#include <iostream>


using namespace bpp;
using namespace std;

static bool enableDotOutput = true;
static void dotOutput(const std::string& testName, const std::vector<const Node_DF*>& nodes)
{
  if (enableDotOutput)
  {
    using bpp::DotOptions;
    writeGraphToDot(
      "debug_" + testName + ".dot", nodes);//, DotOptions::DetailedNodeInfo | DotOptions::ShowDependencyIndex);
  }
}


int main() {
  

  Newick reader;
  shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((A:0.01, B:0.02):0.03,C:0.01):0.01,D:0.1);", false, "", false, false));

  //-------------

  shared_ptr<const NucleicAlphabet> nucAlphabet = AlphabetTools::DNA_ALPHABET;
  shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET;

  Pasta pasta;
  auto sites = make_shared<ProbabilisticVectorSiteContainer>(alphabet);
  pasta.readAlignment("counts.pa", *sites);
  

  // model
  auto t92 = make_shared<T92>(nucAlphabet, 3.);

  auto multimodel = make_shared<MultinomialFromTransitionModel>(t92); // t92 is copied there


  auto rdist = make_shared<ConstantRateDistribution>();
  auto process = make_shared<NonHomogeneousSubstitutionProcess> (rdist, pTree);

  
  // internal leaves
  process->addModel(t92, {2,4}); // internal branches
  process->addModel(multimodel, {0,1,3,5}); // leaves

  if (!process->isFullySetUp())
    throw Exception("test_likelihood_multinomial: process not fully set up.");

  process->getParameters().printParameters(cerr);

  process->aliasParameters("MultinomialFrom.T92.kappa_2","T92.kappa_1");
  process->aliasParameters("MultinomialFrom.T92.theta_2","T92.theta_1");

  process->getIndependentParameters().printParameters(cerr);

  cerr << endl;
  
  Context context;
  
  auto lik = make_shared<LikelihoodCalculationSingleProcess>(context, sites, process);
  auto newtl = make_shared<SingleProcessPhyloLikelihood>(context, lik);

  cerr << "StartLik: " << setprecision(20) << newtl->getValue() << endl;
  newtl->getParameters().printParameters(cout);
  dotOutput("lik_multinomial_value", {lik->getLikelihoodNode().get()});

  OptimizationTools::optimizeNumericalParameters2(
      newtl,
      newtl->getParameters(),
      0, 0.000001, 1000, 0, 0);
  cout << "NewLik: " << newtl->getValue() << endl;
  newtl->getParameters().printParameters(cout);

  return 0;
}


