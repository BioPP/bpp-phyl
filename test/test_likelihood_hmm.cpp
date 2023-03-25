//
// File: test_likelihood_hmm.cpp
// Created by: laurent Guéguen
// Created on: samedi 5 septembre 2020, à 10h 50
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
#include <Bpp/Numeric/Prob/Simplex.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/OptimizationTools.h>

#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>

#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/SubstitutionProcessCollection.h>
#include <Bpp/Phyl/Likelihood/MixtureSequenceEvolution.h>

#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/HmmProcessPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/HmmOfAlignedPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/AutoCorrelationProcessPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/AutoCorrelationOfAlignedPhyloLikelihood.h>

#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>

#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  Context context;
  
  Newick reader;
  auto tree1 = reader.parenthesisToPhyloTree("(((A:0.1, B:0.2):0.3,C:0.1):0.2,D:0.3);");
  auto tree2 = reader.parenthesisToPhyloTree("((A:0.05, C:0.02):0.1,(D:0.01,B:0.03):0.05);");

  vector<shared_ptr<PhyloNode>> vl = tree1->getAllLeaves();
   
  vector<string> seqNames;
  for (auto& vi : vl)
    seqNames.push_back(vi->getName());
  
  vector<unsigned int> ids = tree1->getNodeIndexes(tree1->getAllNodes());
  //-------------

  shared_ptr<const NucleicAlphabet> nucAlphabet = AlphabetTools::DNA_ALPHABET;
  shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET;
  auto rootFreqs = make_shared<GCFrequencySet>(nucAlphabet);
  
  auto model1 = make_shared<T92>(nucAlphabet, 3.,0.9);
  auto model2 = make_shared<T92>(nucAlphabet, 2., 0.1);
  auto model3 = make_shared<T92>(nucAlphabet, 5., 0.5);

  auto rdist1 = make_shared<ConstantRateDistribution>();//GammaDiscreteRateDistribution>(4, 2.0);
  auto rdist2 = rdist1;//std::make_shared<GammaDiscreteRateDistribution>(3, 1.0);
  
  auto parTree1 = make_shared<ParametrizablePhyloTree>(*tree1);
  auto parTree2 = make_shared<ParametrizablePhyloTree>(*tree2);

  ///////////////////////////////////////////
  // Collection Processes

  auto modelColl = make_shared<SubstitutionProcessCollection>();
  
  modelColl->addModel(model1, 1);
  modelColl->addModel(model2, 2);
  modelColl->addModel(model3, 3);

  modelColl->addFrequencies(rootFreqs, 1);
  modelColl->addDistribution(rdist1, 1);
  modelColl->addDistribution(rdist2, 2);

  modelColl->addTree(parTree1, 1);
  modelColl->addTree(parTree2, 2);

  /////////////////////////////////////////
  // First Process

  Vuint vP1m1{0, 3, 4};
  Vuint vP1m2{1, 2, 5};

  map<size_t, Vuint> mModBr1;
  mModBr1[1]=vP1m1;
  mModBr1[2]=vP1m2;

  modelColl->addSubstitutionProcess(1, mModBr1, 1, 1);
                                   
  ///////////////////////////////////////////
  // Second Process

  Vuint vP2m1{0, 1, 3};
  Vuint vP2m2{2, 4, 5};

  map<size_t, Vuint> mModBr2;
  mModBr2[1]=vP2m1;
  mModBr2[3]=vP2m2;

  modelColl->addSubstitutionProcess(2, mModBr2, 2, 2, 1);

  // Data

  auto sites = make_shared<VectorSiteContainer>(alphabet);
  auto seqA = make_unique<Sequence>("A", "ATCCAGACATGCCGGGACTTTGCAGAGAAGGAGTTGTTTCCCATTGCAGCCCAGGTGGATAAGGAACAGC", alphabet);
  sites->addSequence("A", seqA);
  auto seqB = make_unique<Sequence>("B", "CGTCAGACATGCCGTGACTTTGCCGAGAAGGAGTTGGTCCCCATTGCGGCCCAGCTGGACAGGGAGCATC", alphabet);
  sites->addSequence("B", seqB);
  auto seqC = make_unique<Sequence>("C", "GGTCAGACATGCCGGGAATTTGCTGAAAAGGAGCTGGTTCCCATTGCAGCCCAGGTAGACAAGGAGCATC", alphabet);
  sites->addSequence("C", seqC);
  auto seqD = make_unique<Sequence>("D", "TTCCAGACATGCCGGGACTTTACCGAGAAGGAGTTGTTTTCCATTGCAGCCCAGGTGGATAAGGAACATC", alphabet);
  sites->addSequence("D", seqD);

  //  Hmm of process
  
  vector<size_t> vp(2);
  vp[0]=1; vp[1]=2;

  auto hse = make_shared<AutoCorrelationSequenceEvolution>(modelColl, vp);

  auto collNodes = make_shared<CollectionNodes>(context, modelColl);
  
  auto hppl = make_shared<AutoCorrelationProcessPhyloLikelihood>(sites, hse, collNodes);

  hppl->getParameters().printParameters(std::cerr);

  using bpp::DotOptions;
  bpp::writeGraphToDot("hppl.dot", {hppl->getLikelihoodNode().get()});//, DotOptions::DetailedNodeInfo | DotOp
  cerr << "Hppl: " << hppl->getValue() << endl;


  // Derivative Graph

  // Manual access to dkappa
  auto dkappa = dynamic_cast<ConfiguredParameter*>(hppl->getSharedParameter("T92.kappa_3").get());

  auto dlogLik_dkappa = hppl->getLikelihoodNode()->deriveAsValue(context, *dkappa->dependency(0));

  bpp::writeGraphToDot("dhppl_kappa3.dot", {dlogLik_dkappa.get()});
  
  // check posterior probabilities

  auto postprob = hppl->getPosteriorProbabilitiesPerSitePerProcess();

  for (auto& v : postprob)
  {
    VectorTools::print(v);
    std::cerr << VectorTools::sum(v) << std::endl;
  }

  cout << endl;

  cout << "==========================================" << endl;
  cout << "==========================================" << endl;
  cout << endl;
  
  cout << "Optimization : " << endl;
  cout << endl;

  auto profiler  = make_shared<StlOutputStream>(make_unique<ofstream>("profile.txt", ios::out));
  auto messenger = make_shared<StlOutputStream>(make_unique<ofstream>("messages.txt", ios::out));

  unsigned int cM = OptimizationTools::optimizeNumericalParameters2(
    hppl, hppl->getParameters(), 0,
    0.0001, 10000,
    messenger, profiler,
    false, false,
    1, OptimizationTools::OPTIMIZATION_NEWTON);
  
  cerr << "Opt M rounds: " << cM << endl;

  cerr << "Hppl: " << hppl->getValue() << endl;

  hppl->getParameters().printParameters(std::cout);

  return 0;
}
