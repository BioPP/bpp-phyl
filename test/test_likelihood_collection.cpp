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
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/FrequenciesSet/NucleotideFrequenciesSet.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/OptimizationTools.h>

#include <Bpp/Phyl/NewLikelihood/ParametrizableTree.h>

#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/SubstitutionProcessCollection.h>
#include <Bpp/Phyl/NewLikelihood/MixturePhyloLikelihood.h>

#include <iostream>

using namespace bpp;
using namespace newlik;
using namespace std;

int main() {
  TreeTemplate<Node>* tree1 = TreeTemplateTools::parenthesisToTree("(((A:0.1, B:0.2):0.3,C:0.1):0.2,D:0.3);");
  TreeTemplate<Node>* tree2 = TreeTemplateTools::parenthesisToTree("((A:0.1, C:0.2):0.2,(D:0.1,B:0.2):0.1);");

  vector<string> seqNames= tree1->getLeavesNames();
  vector<int> ids = tree1->getNodesId();
  //-------------

  const NucleicAlphabet* alphabet = &AlphabetTools::DNA_ALPHABET;
  FrequenciesSet* rootFreqs = new GCFrequenciesSet(alphabet);
  
  SubstitutionModel* model1 = new T92(alphabet, 3.,0.9);
  SubstitutionModel* model2 = new T92(alphabet, 2., 0.1);
  SubstitutionModel* model3 = new T92(alphabet, 5., 0.5);

  DiscreteDistribution* rdist1 = new GammaDiscreteRateDistribution(4, 2.0);
  DiscreteDistribution* rdist2 = new GammaDiscreteRateDistribution(3, 1.0);
  
  ParametrizableTree* parTree1 = new ParametrizableTree(*tree1);
  ParametrizableTree* parTree2 = new ParametrizableTree(*tree2);

  // First Process

  NonHomogeneousSubstitutionProcess* subPro1=new NonHomogeneousSubstitutionProcess(rdist1->clone(), parTree1->clone());

  int P1[] = {0,3,4,5,1,2};
  Vint vP1m1(&P1[0], &P1[4]);
  Vint vP1m2(&P1[4], &P1[sizeof(P1)/sizeof(P1[0])]);

  subPro1->addModel(model1->clone(),vP1m1);
  subPro1->addModel(model2->clone(),vP1m2);
  
  // Second Process

  NonHomogeneousSubstitutionProcess* subPro2= new NonHomogeneousSubstitutionProcess(rdist2->clone(), parTree2->clone(), rootFreqs->clone());

  int P2[] = {0,1,2,3,4,5};
  Vint vP2m1(&P2[0], &P2[5]);
  Vint vP2m2(&P2[5], &P2[sizeof(P2)/sizeof(P2[0])]);

  subPro2->addModel(model1->clone(),vP2m1);
  subPro2->addModel(model3->clone(),vP2m2);

  // Similar Collection Process

  SubstitutionProcessCollection* modelColl=new SubstitutionProcessCollection();
  
  modelColl->addModel(model1, 1);
  modelColl->addModel(model2, 2);
  modelColl->addModel(model3, 3);
  modelColl->addFrequencies(rootFreqs, 1);
  modelColl->addDistribution(rdist1, 1);
  modelColl->addDistribution(rdist2, 2);

  modelColl->addTree(parTree1, 1);
  modelColl->addTree(parTree2, 2);

  map<size_t, Vint> mModBr1;
  mModBr1[1]=vP1m1;
  mModBr1[2]=vP1m2;

  modelColl->addSubstitutionProcess(mModBr1, 1, 1);
                                   
  map<size_t, Vint> mModBr2;
  mModBr2[1]=vP2m1;
  mModBr2[3]=vP2m2;

  modelColl->addSubstitutionProcess(mModBr2, 2, 2, 1);

  // Data

  VectorSiteContainer sites(alphabet);
  sites.addSequence(BasicSequence("A", "GAACACGAAAGCATGAATGTTCAGTGAGTAGATCAAATATGTCATTTCTGAATTATTATA", alphabet));
  sites.addSequence(BasicSequence("B", "TTTGAACTGTTTGAATATAAGAAAGTTAAATATCTTATAACCAAGTAATATGTTTTAAGA", alphabet));
  sites.addSequence(BasicSequence("C", "GTAATACTTTATAAATACTGATCAATTCAGATAATTTTCAGAACTAACATATATATTATG", alphabet));
  sites.addSequence(BasicSequence("D", "TCGATCGAAAGCCAGGATCAACAATCTTTAACTTATATCGAAAATCATTTATGTGAAGGC", alphabet));

  // Likelihoods

  SubstitutionProcess* sP1c=subPro1->clone();
  SubstitutionProcess* sP2c=subPro2->clone();

  SingleRecursiveTreeLikelihoodCalculation* rtl1=new SingleRecursiveTreeLikelihoodCalculation(*sites.clone(), sP1c, true, true);
  SinglePhyloLikelihood spl1(sP1c, rtl1, true);
  spl1.computeTreeLikelihood();
  
  SingleRecursiveTreeLikelihoodCalculation* rtl2=new SingleRecursiveTreeLikelihoodCalculation(*sites.clone(), sP2c, true, true);
  SinglePhyloLikelihood spl2(sP2c, rtl2, true);
  spl2.computeTreeLikelihood();
  
  cerr << setprecision(10) << "TL1:"  << spl1.getValue() << "\tTL2:" << spl2.getValue() << endl;

  MixturePhyloLikelihood mlc(*sites.clone(), modelColl);

  cerr << "Mlc: " << mlc.getValue() << endl;

  for (size_t pos=0; pos < sites.getNumberOfSites(); pos++){
    double x=spl1.getLikelihoodForASite(pos) * mlc.getSubProcessProb(0) + spl2.getLikelihoodForASite(pos) * mlc.getSubProcessProb(1);
    if (abs(x-mlc.getLikelihoodForASite(pos))>0.001)
      cerr << "Problem on site " << x << endl;
  }

  OutputStream* profiler  = new StlOutputStream(new ofstream("profile.txt", ios::out));
  OutputStream* messenger = new StlOutputStream(new ofstream("messages.txt", ios::out));

  unsigned int c1 = OptimizationTools::optimizeNumericalParameters2(
   &spl1, spl1.getParameters(), 0,
   0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

  cerr << "Opt 1: " << c1 << endl;
  unsigned int c2 = OptimizationTools::optimizeNumericalParameters2(
         &spl2, spl2.getParameters(), 0,
         0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

  cerr << "Opt 2: " << c2 << endl;
  unsigned int cM = OptimizationTools::optimizeNumericalParameters2(
                                                                    &mlc, mlc.getParameters(), 0,
                                                                    0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

  
  cerr << setprecision(10) << "TL1:"  << spl1.getValue() << "\tTL2:" << spl2.getValue() << endl;
  cerr << "Opt M: " << cM << endl;

  cerr << "Mlc: " << mlc.getValue() << endl;

  mlc.getParameters().printParameters(std::cout);
  
  //   for (size_t i = 0; i < nmodels; ++i) {
  //     cout << modelSet2->getModel(i)->getParameter("theta").getValue() << "\t" << modelSet3->getModel(i)->getParameter("theta").getValue() << endl;
  //     //if (abs(modelSet2->getModel(i)->getParameter("theta").getValue() - modelSet3->getModel(i)->getParameter("theta").getValue()) > 0.1)
  //     //  return 1;
  //     thetasEst1[i] +=  modelSet2->getModel(i)->getParameter("theta").getValue();
  //     thetasEst2[i] +=  modelSet3->getModel(i)->getParameter("theta").getValue();
  //   }
  // }
  // thetasEst1 /= static_cast<double>(nrep);
  // thetasEst2 /= static_cast<double>(nrep);

  // //Now compare estimated values to real ones:
  // for (size_t i = 0; i < thetas.size(); ++i) {
  //    cout << thetas[i] << "\t" << thetasEst1[i] << "\t" << thetasEst2[i] << endl;
  //    double diff1 = abs(thetas[i] - thetasEst1[i]);
  //    double diff2 = abs(thetas[i] - thetasEst2[i]);
  //    if (diff1 > 0.2 || diff2 > 0.2)
  //       return 1;
      //}

  // //-------------
  // delete tree;
  // delete modelSet;
  // delete rdist;

  return 0;
}
