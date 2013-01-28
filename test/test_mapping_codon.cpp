//
// File: test_mapping_codon.cpp
// Created by: Julien Dutheil
// Created on: Sat Mar 19 9:36 2011
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

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/StandardCodonAlphabet.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Codon/CodonRateSubstitutionModel.h>
#include <Bpp/Phyl/Model/FrequenciesSet/CodonFrequenciesSet.h>
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Mapping/SubstitutionRegister.h>
#include <Bpp/Phyl/Mapping/NaiveSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/LaplaceSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/UniformizationSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/DecompositionSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/ProbabilisticSubstitutionMapping.h>
#include <Bpp/Phyl/Mapping/SubstitutionMappingTools.h>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree("((A:0.001, B:0.002):0.008,C:0.01,D:0.1);");
  vector<int> ids = tree->getNodesId();
  ids.pop_back(); //Ignore root

  //-------------

  CodonAlphabet* alphabet = new StandardCodonAlphabet(&AlphabetTools::DNA_ALPHABET);
  GeneticCode* gc = new StandardGeneticCode(&AlphabetTools::DNA_ALPHABET);
  //SubstitutionModel* model = new YN98(gc, CodonFrequenciesSet::getFrequenciesSetForCodons(CodonFrequenciesSet::F0, *alphabet));
  SubstitutionModel* model = new CodonRateSubstitutionModel(
        dynamic_cast<const CodonAlphabet*>(gc->getSourceAlphabet()),
        new JCnuc(dynamic_cast<CodonAlphabet*>(alphabet)->getNucleicAlphabet()));
  MatrixTools::printForR(model->getGenerator(), "g");
  DiscreteDistribution* rdist = new ConstantDistribution(1.0);
  HomogeneousSequenceSimulator simulator(model, rdist, tree);
  TotalSubstitutionRegister* totReg = new TotalSubstitutionRegister(alphabet);
  DnDsSubstitutionRegister* dndsReg = new DnDsSubstitutionRegister(gc);

  unsigned int n = 20000;
  vector< vector<double> > realMap(n);
  vector< vector< vector<double> > > realMapTotal(n);
  vector< vector< vector<double> > > realMapDnDs(n);
  VectorSiteContainer sites(tree->getLeavesNames(), alphabet);
  for (unsigned int i = 0; i < n; ++i) {
    ApplicationTools::displayGauge(i, n-1, '=');
    RASiteSimulationResult* result = simulator.dSimulate();
    realMap[i].resize(ids.size());
    realMapTotal[i].resize(ids.size());
    realMapDnDs[i].resize(ids.size());
    for (size_t j = 0; j < ids.size(); ++j) {
      realMap[i][j] = static_cast<double>(result->getSubstitutionCount(ids[j]));
      realMapTotal[i][j].resize(totReg->getNumberOfSubstitutionTypes());
      realMapDnDs[i][j].resize(dndsReg->getNumberOfSubstitutionTypes());
      result->getSubstitutionCount(ids[j], *totReg, realMapTotal[i][j]);
      result->getSubstitutionCount(ids[j], *dndsReg, realMapDnDs[i][j]);
      if (realMapTotal[i][j][0] != realMap[i][j]) {
        cerr << "Error, total substitution register provides wrong result." << endl;
        return 1;
      }
      //if (abs(VectorTools::sum(realMapDetailed[i][j]) - realMap[i][j]) > 0.000001) {
      //  cerr << "Error, detailed substitution register provides wrong result." << endl;
      //  return 1;
      //}
    }
    auto_ptr<Site> site(result->getSite());
    site->setPosition(i);
    sites.addSite(*site, false);
    delete result;
  }
  ApplicationTools::displayTaskDone();
  
  //-------------
  //Now build the substitution vectors with the true model:
  //Fasta fasta;
  //fasta.write("Simulations.fasta", sites);
  DRHomogeneousTreeLikelihood drhtl(*tree, sites, model, rdist);
  drhtl.initialize();
  cout << drhtl.getValue() << endl;
 
  SubstitutionCount* sCountAna = new LaplaceSubstitutionCount(model, 10);
  Matrix<double>* m = sCountAna->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Analytical total count:" << endl;
  MatrixTools::print(*m);
  delete m;
  ProbabilisticSubstitutionMapping* probMapAna = 
    SubstitutionMappingTools::computeSubstitutionVectors(drhtl, *sCountAna);

  SubstitutionCount* sCountTot = new NaiveSubstitutionCount(totReg);
  m = sCountTot->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Simple total count:" << endl;
  MatrixTools::print(*m);
  delete m;
  ProbabilisticSubstitutionMapping* probMapTot = 
    SubstitutionMappingTools::computeSubstitutionVectors(drhtl, *sCountTot);

  SubstitutionCount* sCountDnDs = new NaiveSubstitutionCount(dndsReg);
  m = sCountDnDs->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Detailed count, type 1:" << endl;
  MatrixTools::print(*m);
  delete m;
  ProbabilisticSubstitutionMapping* probMapDnDs = 
    SubstitutionMappingTools::computeSubstitutionVectors(drhtl, *sCountDnDs);

  SubstitutionCount* sCountUniTot = new UniformizationSubstitutionCount(model, totReg);
  m = sCountUniTot->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Total count, uniformization method:" << endl;
  MatrixTools::print(*m);
  delete m;
  ProbabilisticSubstitutionMapping* probMapUniTot = 
    SubstitutionMappingTools::computeSubstitutionVectors(drhtl, *sCountUniTot);

  SubstitutionCount* sCountUniDnDs = new UniformizationSubstitutionCount(model, dndsReg);
  m = sCountUniDnDs->getAllNumbersOfSubstitutions(0.001,2);
  cout << "Detailed count, uniformization method, type 2:" << endl;
  MatrixTools::print(*m);
  delete m;
  ProbabilisticSubstitutionMapping* probMapUniDnDs = 
    SubstitutionMappingTools::computeSubstitutionVectors(drhtl, *sCountUniDnDs);

  //Check per branch:
  
  /*
  //1. Total:
  for (unsigned int j = 0; j < ids.size(); ++j) {
    double totalReal = 0;
    double totalObs1 = 0;
    double totalObs2 = 0;
    double totalObs3 = 0;
    double totalObs4 = 0;
    double totalObs5 = 0;
    for (unsigned int i = 0; i < n; ++i) {
      totalReal += realMap[i][j];
      totalObs1 += probMapAna->getNumberOfSubstitutions(ids[j], i, 0);
      totalObs2 += probMapTot->getNumberOfSubstitutions(ids[j], i, 0);
      //totalObs3 += VectorTools::sum(probMapDet->getNumberOfSubstitutions(ids[j], i));
      totalObs4 += probMapDecTot->getNumberOfSubstitutions(ids[j], i, 0);
      //totalObs5 += VectorTools::sum(probMapDecDet->getNumberOfSubstitutions(ids[j], i));
    }
    if (tree->isLeaf(ids[j])) cout << tree->getNodeName(ids[j]) << "\t";
    cout << tree->getDistanceToFather(ids[j]) << "\t" << totalReal << "\t" << totalObs1 << "\t" << totalObs2 << "\t" << totalObs3 << "\t" << totalObs4 << "\t" << totalObs5 << endl;
    if (abs(totalReal - totalObs1) / totalReal > 0.1) return 1;
    if (abs(totalReal - totalObs2) / totalReal > 0.1) return 1;
    if (abs(totalReal - totalObs3) / totalReal > 0.1) return 1;
    if (abs(totalReal - totalObs4) / totalReal > 0.1) return 1;
  }
  //2. Detail:
  for (unsigned int j = 0; j < ids.size(); ++j) {
    vector<double> real(4, 0);
    vector<double> obs1(4, 0);
    vector<double> obs2(4, 0);
    for (unsigned int i = 0; i < n; ++i) {
      real += realMapDetailed[i][j];
      //VectorTools::print(real);
      //vector<double> c = probMapDet->getNumberOfSubstitutions(ids[j], i);
      //VectorTools::print(c);
      //obs1 += probMapDet->getNumberOfSubstitutions(ids[j], i);
      //obs2 += probMapDecDet->getNumberOfSubstitutions(ids[j], i);
    }
    if (tree->isLeaf(ids[j])) cout << tree->getNodeName(ids[j]) << "\t";
    cout << tree->getDistanceToFather(ids[j]) << "\t";
    for (unsigned int t = 0; t < 4; ++t) {
      cout << obs1[t] << "/" << real[t] << "\t";
      cout << obs2[t] << "/" << real[t] << "\t";
    }
    cout << endl;
    //if (abs(totalReal - totalObs) / totalReal > 0.1) return 1;
  }
  */

  //-------------
  delete tree;
  delete alphabet;
  delete model;
  delete rdist;
  delete sCountTot;
  delete sCountDnDs;
  delete probMapAna;
  delete probMapTot;
  delete probMapDnDs;
  delete probMapUniTot;
  delete probMapUniDnDs;
  //return (abs(obs - 0.001) < 0.001 ? 0 : 1);
  return 0;
}
