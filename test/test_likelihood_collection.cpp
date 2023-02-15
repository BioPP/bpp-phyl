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

#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/MixtureProcessPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/MixtureOfAlignedPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/FormulaOfPhyloLikelihood.h>

#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  Context context;
  
  Newick reader;

  shared_ptr<PhyloTree> tree1(reader.parenthesisToPhyloTree("(((A:0.1, B:0.2):0.3,C:0.1):0.2,D:0.3);"));
  shared_ptr<PhyloTree> tree2(reader.parenthesisToPhyloTree("((A:0.05, C:0.02):0.1,(D:0.01,B:0.03):0.05);"));

  vector<shared_ptr<PhyloNode> > vl= tree1->getAllLeaves();
  
  auto parTree1 = std::make_shared<ParametrizablePhyloTree>(*tree1);
  auto parTree2 = std::make_shared<ParametrizablePhyloTree>(*tree2);

  vector<string> seqNames;
  for (size_t i=0; i<vl.size(); i++)
    seqNames.push_back(vl[i]->getName());
  
  vector<unsigned int> ids = tree1->getNodeIndexes(tree1->getAllNodes());
  //-------------

  shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET;
  shared_ptr<const NucleicAlphabet> nucAlphabet = AlphabetTools::DNA_ALPHABET;
  auto rootFreqs = make_shared<GCFrequencySet>(nucAlphabet);
  
  auto model1 = make_shared<T92>(nucAlphabet, 3.,0.9);
  auto model2 = make_shared<T92>(nucAlphabet, 2., 0.1);
  auto model3 = make_shared<T92>(nucAlphabet, 5., 0.5);

  auto rdist1 = make_shared<ConstantRateDistribution>();//GammaDiscreteRateDistribution>(4, 2.0);
  auto rdist2 = rdist1;//std::make_shared<GammaDiscreteRateDistribution>(3, 1.0);
  
  /////////////////////////////////////////
  // First Process

  auto subPro1 = make_shared<NonHomogeneousSubstitutionProcess>(
      shared_ptr<DiscreteDistribution>(rdist1->clone()),
      shared_ptr<PhyloTree>(tree1->clone()));

  Vuint vP1m1{0, 3, 4};
  Vuint vP1m2{1, 2, 5};

  subPro1->addModel(std::shared_ptr<T92>(model1->clone()), vP1m1);
  subPro1->addModel(std::shared_ptr<T92>(model2->clone()), vP1m2);

  ///////////////////////////////////////////
  // Second Process

  auto subPro2 = make_shared<NonHomogeneousSubstitutionProcess>(
      shared_ptr<DiscreteDistribution>(rdist2->clone()),
      shared_ptr<PhyloTree>(tree2->clone()),
      shared_ptr<FrequencySetInterface>(rootFreqs->clone()));
  
  Vuint vP2m1{0, 1, 3};
  Vuint vP2m2{2, 4, 5};

  subPro2->addModel(std::shared_ptr<T92>(model1->clone()), vP2m1);
  subPro2->addModel(std::shared_ptr<T92>(model3->clone()), vP2m2);

  ///////////////////////////////////////////
  // Similar Collection Processes

  auto modelColl = make_shared<SubstitutionProcessCollection>();
  
  modelColl->addModel(model1, 1);
  modelColl->addModel(model2, 2);
  modelColl->addModel(model3, 3);

  modelColl->addFrequencies(rootFreqs, 1);
  modelColl->addDistribution(rdist1, 1);
  modelColl->addDistribution(rdist2, 2);

  modelColl->addTree(parTree1, 1);
  modelColl->addTree(parTree2, 2);

  map<size_t, Vuint> mModBr1;
  mModBr1[1] = vP1m1;
  mModBr1[2] = vP1m2;

  modelColl->addSubstitutionProcess(1, mModBr1, 1, 1);
                                   
  map<size_t, Vuint> mModBr2;
  mModBr2[1] = vP2m1;
  mModBr2[3] = vP2m2;

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

  // Likelihoods
  auto pc = make_shared<PhyloLikelihoodContainer>(context, modelColl);

  shared_ptr<SubstitutionProcessInterface> sP1c(subPro1->clone());
  shared_ptr<SubstitutionProcessInterface> sP2c(subPro2->clone());

  auto lik1 = make_shared<LikelihoodCalculationSingleProcess>(context, sites, sP1c);
  pc->addPhyloLikelihood(1, make_shared<SingleProcessPhyloLikelihood>(context, lik1));
    
  auto lik2 = make_shared<LikelihoodCalculationSingleProcess>(context, sites, sP2c);
  pc->addPhyloLikelihood(2, make_shared<SingleProcessPhyloLikelihood>(context, lik2));
  
  auto spl1 = dynamic_pointer_cast<AlignedPhyloLikelihoodInterface>((*pc)[1]);
  auto spl2 = dynamic_pointer_cast<AlignedPhyloLikelihoodInterface>((*pc)[2]);

  cerr << setprecision(10) << "TL1:"  << spl1->getValue() << "\tTL2:" << spl2->getValue() << endl;

  auto collNodes = pc->getCollectionNodes();
  
  //  Mixture of process
  
  std::vector<size_t> vp(2);
  vp[0]=1; vp[1]=2;

  auto mse = make_shared<MixtureSequenceEvolution>(modelColl, vp);

  auto mlc = make_shared<MixtureProcessPhyloLikelihood>(sites, mse, collNodes);

  using bpp::DotOptions;
  bpp::writeGraphToDot("mlc.dot", {mlc->getLikelihoodNode().get()});//, DotOptions::DetailedNodeInfo | DotOp
  cerr << "Mlc: " << mlc->getValue() << endl;

  for (size_t pos=0; pos < sites->getNumberOfSites(); pos++){
    DataLik x = spl1->getLikelihoodForASite(pos) * mlc->getSubProcessProb(0) + spl2->getLikelihoodForASite(pos) * mlc->getSubProcessProb(1);
    if (convert(abs(x-mlc->getLikelihoodForASite(pos))) > 0.001)
      cerr << "Mixture Process : Problem on site " << x << endl;
  }

  //  Mixture of phylo

  // New calculation to ensure independence
  
  Context context2;
  
  auto pc2 = make_shared<PhyloLikelihoodContainer>(context2, modelColl);

  auto sP1c2 = shared_ptr<SubstitutionProcessInterface>(subPro1->clone());
  auto sP2c2 = shared_ptr<SubstitutionProcessInterface>(subPro2->clone());

  auto lik12 = make_shared<LikelihoodCalculationSingleProcess>(context2, sites, sP1c2);
  pc2->addPhyloLikelihood(2, make_shared<SingleProcessPhyloLikelihood>(context2, lik12));
  
  auto lik22 = make_shared<LikelihoodCalculationSingleProcess>(context2, sites, sP2c2);
  pc2->addPhyloLikelihood(1, make_shared<SingleProcessPhyloLikelihood>(context2, lik22));

  vector<size_t> nPhylo = {1,2};
  auto moap = make_shared<MixtureOfAlignedPhyloLikelihood>(context2, pc2, nPhylo, false);

  bpp::writeGraphToDot("moap.dot", {moap->getLikelihoodNode().get()});//, DotOptions::DetailedNodeInfo | DotOp
  cerr << "Moap: " << moap->getValue() << endl;

  for (size_t pos=0; pos < sites->getNumberOfSites(); pos++){
    DataLik x = spl1->getLikelihoodForASite(pos) * moap->getPhyloProb(0) + spl2->getLikelihoodForASite(pos) * moap->getPhyloProb(1);
    if (convert(abs(x-moap->getLikelihoodForASite(pos))) > 0.001)
      cerr << "Mixture Alignment: Problem on site " << x << endl;
  }

  cout << endl;
  
  cout << "==========================================" << endl;
  cout << "==========================================" << endl;
  cout << endl;
  
  cout << "Optimization : " << endl;
  cout << endl;

  auto profiler  = make_shared<StlOutputStream>(new ofstream("profile.txt", ios::out));
  auto messenger = make_shared<StlOutputStream>(new ofstream("messages.txt", ios::out));

  unsigned int c1 = OptimizationTools::optimizeNumericalParameters2(
    spl1, spl1->getParameters(), 0,
    0.0001, 10000,
    messenger, profiler,
    false, false,
    1, OptimizationTools::OPTIMIZATION_NEWTON);

  cerr << "Opt 1: rounds " << c1 << endl;

  cerr << "--------------------------------" << endl;
  
  unsigned int c2 = OptimizationTools::optimizeNumericalParameters2(
    spl2, spl2->getParameters(), 0,
    0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

  spl1->getParameters().printParameters(std::cout);
  
  cerr << "Opt 2: rounds " << c2 << endl;

  spl2->getParameters().printParameters(std::cout);

  
  cerr << setprecision(10) << "Ml1:"  << spl1->getValue() << "\tMl2:" << spl2->getValue() << endl;

  cerr << "--------------------------------" << endl;

  unsigned int cM = OptimizationTools::optimizeNumericalParameters2(
    mlc, mlc->getParameters(), 0,
    0.0001, 10000,
    messenger, profiler,
    false, false,
    1, OptimizationTools::OPTIMIZATION_NEWTON);
  
  cerr << "Opt M rounds: " << cM << endl;

  cerr << "Mlc: " << mlc->getValue() << endl;

  mlc->getParameters().printParameters(std::cout);

  cerr << "--------------------------------" << endl;
  
  unsigned int cM2 = OptimizationTools::optimizeNumericalParameters2(
    moap, moap->getParameters(), 0,
    0.0001, 10000,
    messenger, profiler,
    false, false,
    1, OptimizationTools::OPTIMIZATION_NEWTON);
  
  cerr << "Opt MOAP rounds: " << cM2 << endl;

  cerr << "Moap: " << moap->getValue() << endl;

  moap->getParameters().printParameters(std::cout);

  
// Formula
  
  string formula="(phylo2 - phylo1) * (phylo1 - phylo2)";      

  auto tl = make_shared<FormulaOfPhyloLikelihood>(context, pc, formula, false);

  cerr << formula << " : " << tl->getValue() << endl;

  bpp::writeGraphToDot("formula.dot", {tl->getLikelihoodNode().get()});//, DotOptions::DetailedNodeInfo | DotOp

  unsigned int cMtl = OptimizationTools::optimizeNumericalParameters2(
    tl, tl->getParameters(), 0,
    0.0001, 10000,
    messenger, profiler,
    false, false,
    1, OptimizationTools::OPTIMIZATION_NEWTON);

  cerr << "Opt tl rounds: " << cMtl << endl;

  cerr << formula << " " << tl->getValue() << endl;

  return 0;
}
