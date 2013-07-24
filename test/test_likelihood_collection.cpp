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
#include <Bpp/Phyl/Model/RateDistribution.all>
#include <Bpp/Phyl/OptimizationTools.h>

#include <Bpp/Phyl/NewLikelihood/ParametrizableTree.h>

#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/SubstitutionProcessCollection.h>
#include <Bpp/Phyl/NewLikelihood/SingleRecursiveTreeLikelihoodCalculation.h>
#include <Bpp/Phyl/NewLikelihood/MixtureLikelihoodCollection.h>

#include <iostream>

using namespace bpp;
using namespace newlik;
using namespace std;

int main() {
  TreeTemplate<Node>* tree1 = TreeTemplateTools::parenthesisToTree("(((A:0.1, B:0.2):0.3,C:0.1):0.2,(D:0.3,(E:0.2,F:0.05):0.1):0.1);");
  TreeTemplate<Node>* tree2 = TreeTemplateTools::parenthesisToTree("(((A:0.1, C:0.2):0.4,E:0.5):0.2,(D:0.1,(B:0.2,F:0.15):0.1):0.1);");

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

  int P1m1[] = {0,1,2,3,4,5};
  Vint vP1m1(&P1m1[0], &P1m1[6]);
  int P1m2[] = {6,7,8,9};
  Vint vP1m2(&P1m2[0], &P1m2[4]);

  subPro1->addModel(model1->clone(),vP1m1);
  subPro1->addModel(model2->clone(),vP1m2);
  
  // Second Process

  NonHomogeneousSubstitutionProcess* subPro2= new NonHomogeneousSubstitutionProcess(rdist2->clone(), parTree2->clone());

  int P2m1[] = {0,1,2,3,4};
  Vint vP2m1(&P2m1[0], &P2m1[5]);
  int P2m2[] = {5,6,7,8,9};
  Vint vP2m2(&P2m2[0], &P2m2[5]);

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

  ParameterList pl=modelColl->getParameters();
  for (size_t i=0; i<pl.size(); i++)
    cout << pl[i].getName() << " " << pl[i].getValue() << endl;

  map<size_t, Vint> mModBr1;
  mModBr1[1]=vP1m1;
  mModBr1[2]=vP1m2;

  modelColl->addSubstitutionProcess(mModBr1, 1, 1);
                                   
  map<size_t, Vint> mModBr2;
  mModBr2[1]=vP2m1;
  mModBr2[3]=vP2m2;

  modelColl->addSubstitutionProcess(mModBr2, 2, 2);

  // Data

  VectorSiteContainer sites(alphabet);
  sites.addSequence(BasicSequence("A", "AAATG", alphabet));
  sites.addSequence(BasicSequence("B", "GACTG", alphabet));
  sites.addSequence(BasicSequence("C", "CTCTG", alphabet));
  sites.addSequence(BasicSequence("D", "AGATG", alphabet));
  sites.addSequence(BasicSequence("E", "ACGTG", alphabet));
  sites.addSequence(BasicSequence("F", "CAGTT", alphabet));


  // Likelihoods
  
  SingleRecursiveTreeLikelihoodCalculation rtl1(*sites.clone(), subPro1, true, false);

  
  SingleRecursiveTreeLikelihoodCalculation rtl2(*sites.clone(), subPro2, true, true);

  
  //   for (size_t i = 0; i < nmodels; ++i) {
  //     ntl.setParameterValue("T92.theta_" + TextTools::toString(i + 1), thetas[i]);
  //     ntl2.setParameterValue("T92.theta_" + TextTools::toString(i + 1), thetas[i]);
  //   }

   cout << setprecision(10) << "TL1:"  << rtl1.getValue() << "\tTL2:" << rtl2.getValue() << endl;

  //  cerr << modelColl.numberOfSubstitutionProcess() << endl;
  MixtureLikelihoodCollection mlc(*sites.clone(),modelColl);

  //  cout << "Mlc: " << mlc.getValue() << endl;
  
  //   break;
    // unsigned int c1 = OptimizationTools::optimizeNumericalParameters2(
    //     &tl, tl.getParameters(), 0,
    //     0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON)
      ;

    // unsigned int c2 = OptimizationTools::optimizeNumericalParameters2(
    //     &tl2, tl2.getParameters(), 0,
    //     0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);

    // cout << "OldTL: " << c1 << ": " << tl.getValue() << "\t" << c2 << ": " << tl2.getValue() << endl;
                        
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
