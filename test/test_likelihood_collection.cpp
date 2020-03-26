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
#include <Bpp/Phyl/Model/FrequenciesSet/NucleotideFrequenciesSet.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/OptimizationTools.h>

#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>

#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/SubstitutionProcessCollection.h>
#include <Bpp/Phyl/NewLikelihood/MixtureSequenceEvolution.h>

#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/MixtureProcessPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/FormulaOfPhyloLikelihood.h>

#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  Context context;
  
  Newick reader;

  unique_ptr<PhyloTree> tree1(reader.parenthesisToPhyloTree("(((A:0.1, B:0.2):0.3,C:0.1):0.2,D:0.3);"));
  unique_ptr<PhyloTree> tree2(reader.parenthesisToPhyloTree("((A:0.05, C:0.02):0.1,(D:0.01,B:0.03):0.05);"));


  vector<shared_ptr<PhyloNode> > vl= tree1->getAllLeaves();
  
  vector<string> seqNames;
  for (size_t i=0; i<vl.size(); i++)
    seqNames.push_back(vl[i]->getName());
  
  vector<unsigned int> ids = tree1->getNodeIndexes(tree1->getAllNodes());
  //-------------

  const NucleicAlphabet* alphabet = &AlphabetTools::DNA_ALPHABET;
  auto rootFreqs = std::make_shared<GCFrequenciesSet>(alphabet);
  
  auto model1 = std::make_shared<T92>(alphabet, 3.,0.9);
  auto model2 = std::make_shared<T92>(alphabet, 2., 0.1);
  auto model3 = std::make_shared<T92>(alphabet, 5., 0.5);

  auto rdist1 = std::make_shared<GammaDiscreteRateDistribution>(4, 2.0);
  auto rdist2 = std::make_shared<GammaDiscreteRateDistribution>(3, 1.0);
  
  auto parTree1 = std::make_shared<ParametrizablePhyloTree>(*tree1);
  auto parTree2 = std::make_shared<ParametrizablePhyloTree>(*tree2);

  // First Process

  NonHomogeneousSubstitutionProcess* subPro1=new NonHomogeneousSubstitutionProcess(rdist1->clone(), parTree1->clone());

  unsigned int P1[] = {0,3,4,5,1,2};
  Vuint vP1m1(&P1[0], &P1[4]);
  Vuint vP1m2(&P1[4], &P1[sizeof(P1)/sizeof(P1[0])]);

  subPro1->addModel(std::shared_ptr<T92>(model1->clone()),vP1m1);
  subPro1->addModel(std::shared_ptr<T92>(model2->clone()),vP1m2);
  
  // Second Process

  NonHomogeneousSubstitutionProcess* subPro2= new NonHomogeneousSubstitutionProcess(rdist2->clone(), parTree2->clone(), rootFreqs->clone());

  unsigned int P2[] = {0,1,2,3,4,5};
  Vuint vP2m1(&P2[0], &P2[5]);
  Vuint vP2m2(&P2[5], &P2[sizeof(P2)/sizeof(P2[0])]);

  subPro2->addModel(std::shared_ptr<T92>(model1->clone()),vP2m1);
  subPro2->addModel(std::shared_ptr<T92>(model3->clone()),vP2m2);

  // Similar Collection Process

  SubstitutionProcessCollection* modelColl=new SubstitutionProcessCollection();
  
  modelColl->addModel(std::shared_ptr<T92>(model1), 1);
  modelColl->addModel(std::shared_ptr<T92>(model2), 2);
  modelColl->addModel(std::shared_ptr<T92>(model3), 3);
  modelColl->addFrequencies(rootFreqs, 1);
  modelColl->addDistribution(rdist1, 1);
  modelColl->addDistribution(rdist2, 2);

  modelColl->addTree(parTree1, 1);
  modelColl->addTree(parTree2, 2);

  map<size_t, Vuint> mModBr1;
  mModBr1[1]=vP1m1;
  mModBr1[2]=vP1m2;

  modelColl->addSubstitutionProcess(1, mModBr1, 1, 1);
                                   
  map<size_t, Vuint> mModBr2;
  mModBr2[1]=vP2m1;
  mModBr2[3]=vP2m2;

  modelColl->addSubstitutionProcess(2, mModBr2, 2, 2, 1);

  // Data

  VectorSiteContainer sites(alphabet);
  sites.addSequence(BasicSequence("A", "GAACACGAAAGCATGAATGTTCAGTGAGTAGATCAAATATGTCATTTCTGAATTATTATA", alphabet));
  sites.addSequence(BasicSequence("B", "TTTGAACTGTTTGAATATAAGAAAGTTAAATATCTTATAACCAAGTAATATGTTTTAAGA", alphabet));
  sites.addSequence(BasicSequence("C", "GTAATACTTTATAAATACTGATCAATTCAGATAATTTTCAGAACTAACATATATATTATG", alphabet));
  sites.addSequence(BasicSequence("D", "TCGATCGAAAGCCAGGATCAACAATCTTTAACTTATATCGAAAATCATTTATGTGAAGGC", alphabet));

  // Likelihoods
  PhyloLikelihoodContainer pc;

  SubstitutionProcess* sP1c=subPro1->clone();
  SubstitutionProcess* sP2c=subPro2->clone();

  auto lik1 = std::make_shared<LikelihoodCalculationSingleProcess>(context, sites, *sP1c);

  pc.addPhyloLikelihood(1, new SingleProcessPhyloLikelihood(context, lik1));
    
  auto lik2 = std::make_shared<LikelihoodCalculationSingleProcess>(context, sites, *sP2c);

  pc.addPhyloLikelihood(2, new SingleProcessPhyloLikelihood(context, lik2));
  
  AlignedPhyloLikelihood* spl1=dynamic_cast<AlignedPhyloLikelihood*>(pc[1]);
  AlignedPhyloLikelihood* spl2=dynamic_cast<AlignedPhyloLikelihood*>(pc[2]);

  cerr << setprecision(10) << "TL1:"  << spl1->getValue() << "\tTL2:" << spl2->getValue() << endl;


  //  Mixture
  
  std::vector<size_t> vp(2);
  vp[0]=1; vp[1]=2;

  MixtureSequenceEvolution mse(modelColl, vp);

  MixtureProcessPhyloLikelihood mlc(context, *sites.clone(), mse);
  
  cerr << "Mlc: " << mlc.getValue() << endl;

  for (size_t pos=0; pos < sites.getNumberOfSites(); pos++){
    double x=spl1->getLikelihoodForASite(pos) * mlc.getSubProcessProb(0) + spl2->getLikelihoodForASite(pos) * mlc.getSubProcessProb(1);
    if (abs(x-mlc.getLikelihoodForASite(pos))>0.001)
      cerr << "Problem on site " << x << endl;
  }

  cout << endl;
  
  cout << "==========================================" << endl;
  cout << "==========================================" << endl;
  cout << endl;
  
  cout << "Optimization : " << endl;
  cout << endl;

  OutputStream* profiler  = new StlOutputStream(new ofstream("profile.txt", ios::out));
  OutputStream* messenger = new StlOutputStream(new ofstream("messages.txt", ios::out));

  unsigned int c1 = OptimizationTools::optimizeNumericalParameters2(
    spl1, spl1->getParameters(), 0,
    0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

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
    &mlc, mlc.getParameters(), 0,
    0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);
  
  cerr << "Opt M rounds: " << cM << endl;

  cerr << "Mlc: " << mlc.getValue() << endl;

  mlc.getParameters().printParameters(std::cout);


  
// Formula
  
  string formula="phylo1";  
    
  unique_ptr<FormulaOfPhyloLikelihood> tl(new FormulaOfPhyloLikelihood(context, &pc,formula));

  cerr << "tl " << tl->getValue() << endl;

  unsigned int cMtl = OptimizationTools::optimizeNumericalParameters2(
    tl.get(), tl->getParameters(), 0,
    0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

  cerr << "Opt tl rounds: " << cMtl << endl;

  cerr << "tl " << tl->getValue() << endl;

  tl.release();
  

  return 0;
}
