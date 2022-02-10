//
// File: test_mapping.cpp
// Created by: Julien Dutheil
// Created on: Thu Nov 25 16:02 2010
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
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/GTR.h>
#include <Bpp/Phyl/Model/Protein/JTT92.h>
#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSiteSimulator.h>
#include <Bpp/Phyl/Mapping/SubstitutionRegister.h>
#include <Bpp/Phyl/Mapping/SubstitutionCount.h>
#include <Bpp/Phyl/Mapping/LaplaceSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/DecompositionSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/UniformizationSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/NaiveSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/SubstitutionMappingTools.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/SimpleSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Seq/AlphabetIndex/GranthamAAVolumeIndex.h>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  try {
  Newick reader;
  Context context;
  
  shared_ptr<PhyloTree> new_tree(reader.parenthesisToPhyloTree("((A:0.001, B:0.002):0.008,C:0.01,D:0.02);", false, "", false, false));

  vector<uint> ids = {0, 1, 2, 3, 4};

  //-------------

  NucleicAlphabet* alphabet = new DNA();
  auto model = std::make_shared<GTR>(alphabet, 1, 0.2, 0.3, 0.4, 0.4, 0.1, 0.35, 0.35, 0.2);
//  DiscreteDistribution* rdist = new ConstantDistribution(1);
  auto rdist = std::make_shared<GammaDiscreteDistribution>(4, 0.4, 0.4);

  shared_ptr<RateAcrossSitesSubstitutionProcess> process(new RateAcrossSitesSubstitutionProcess(model, std::shared_ptr<DiscreteDistribution>(rdist->clone()), std::shared_ptr<PhyloTree>(new_tree->clone())));

  SimpleSubstitutionProcessSiteSimulator simulator(*process);
  
  TotalSubstitutionRegister* totReg = new TotalSubstitutionRegister(model->getStateMap());
  ComprehensiveSubstitutionRegister* detReg = new ComprehensiveSubstitutionRegister(model->getStateMap());

  size_t n = 50000;
  vector< vector<double> > realMap(n);
  vector< vector< vector<double> > > realMapTotal(n);
  vector< vector< vector<double> > > realMapDetailed(n);
  VectorSiteContainer sites(new_tree->getAllLeavesNames(), alphabet);
  for (size_t i = 0; i < n; ++i) {
    ApplicationTools::displayGauge(i, n - 1, '=');
    unique_ptr<SiteSimulationResult> result(simulator.dSimulateSite());
    realMap[i].resize(ids.size());
    realMapTotal[i].resize(ids.size());
    realMapDetailed[i].resize(ids.size());
    for (size_t j = 0; j < ids.size(); ++j) {
      realMap[i][j] = static_cast<double>(result->getSubstitutionCount(ids[j]));
      realMapTotal[i][j].resize(totReg->getNumberOfSubstitutionTypes());
      realMapDetailed[i][j].resize(detReg->getNumberOfSubstitutionTypes());
      result->getSubstitutionCount(ids[j], *totReg, realMapTotal[i][j]);
      result->getSubstitutionCount(ids[j], *detReg, realMapDetailed[i][j]);
      if (realMapTotal[i][j][0] != realMap[i][j]) {
        throw Exception("Error, total substitution register provides wrong result.");
        return 1;
      }
      if (abs(VectorTools::sum(realMapDetailed[i][j]) - realMap[i][j]) > 0.000001) {
        throw Exception("Error, detailed substitution register provides wrong result.");
        return 1;
      }
    }
    unique_ptr<Site> site(result->getSite(*model));
    site->setPosition(static_cast<int>(i));    
    sites.addSite(*site, false);
  }
  ApplicationTools::displayTaskDone();

  //-------------
  //Now build the substitution vectors with the true model:

  // Newlik

  auto tmComp=make_shared<LikelihoodCalculationSingleProcess>(context, sites, *process);
  SingleProcessPhyloLikelihood newTl(context, tmComp);
  cout << "LogLik: " << newTl.getValue() << endl;
    
  
  
  SubstitutionCount* sCountAna = new LaplaceSubstitutionCount(model.get(), 10);
  Matrix<double>* m = sCountAna->getAllNumbersOfSubstitutions(0.001, 1);
  cout << "Analytical (Laplace) total count:" << endl;
  MatrixTools::print(*m);
  delete m;
  ProbabilisticSubstitutionMapping* probNEWMapAna = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountAna);

  cout << endl;

  //Simple:
  SubstitutionCount* sCountTot = new NaiveSubstitutionCount(model.get(), totReg);
  m = sCountTot->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Simple total count:" << endl;
  MatrixTools::print(*m);
  delete m;
  ProbabilisticSubstitutionMapping* probNEWMapTot = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountTot);
  cout << endl;

  SubstitutionCount* sCountDet = new NaiveSubstitutionCount(model.get(), detReg);
  m = sCountDet->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Detailed count, type 1:" << endl;
  MatrixTools::print(*m);
  delete m;
  ProbabilisticSubstitutionMapping* probNEWMapDet = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountDet);
  cout << endl;

  //Decomposition:
  SubstitutionCount* sCountDecTot = new DecompositionSubstitutionCount(model.get(), totReg);
  m = sCountDecTot->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Total count, decomposition method:" << endl;
  MatrixTools::print(*m);
  delete m;
  ProbabilisticSubstitutionMapping* probNEWMapDecTot = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountDecTot);

  SubstitutionCount* sCountDecDet = new DecompositionSubstitutionCount(model.get(), detReg);
  m = sCountDecDet->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Detailed count, decomposition method, type 1:" << endl;
  MatrixTools::print(*m);
  delete m;
  ProbabilisticSubstitutionMapping* probNEWMapDecDet = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountDecDet);
  cout << endl;

  //Uniformization
  SubstitutionCount* sCountUniTot = new UniformizationSubstitutionCount(model.get(), totReg);
  m = sCountUniTot->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Total count, uniformization method:" << endl;
  MatrixTools::print(*m);
  delete m;
  ProbabilisticSubstitutionMapping* probNEWMapUniTot = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountUniTot);  

  SubstitutionCount* sCountUniDet = new UniformizationSubstitutionCount(model.get(), detReg);
  m = sCountUniDet->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Detailed count, uniformization method, type 1:" << endl;
  MatrixTools::print(*m);
  ProbabilisticSubstitutionMapping* probNEWMapUniDet = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountUniDet);
  cout << endl;
  
  //Check saturation:
  cout << "checking saturation..." << endl;
  double td[] = {0.001, 0.01, 0.1, 1, 2, 3, 4, 10};
  Vdouble vd(td, td+sizeof(td)/sizeof(double));

  for (auto d : vd)
  {
    m = sCountUniDet->getAllNumbersOfSubstitutions(d,1);
    cout << "Total count, uniformization method for " << d << endl;
    MatrixTools::print(*m);
    cout << MatrixTools::sumElements(*m) << endl;
    delete m;
  }
  cout << endl;

  //Check per branch:
  //1. Total:

  for (size_t j = 0; j < ids.size(); ++j) {
    double totalReal = 0;
    double totalObs1 = 0;
    double totalObs2 = 0;
    double totalObs3 = 0;
    double totalObs4 = 0;
    double totalObs5 = 0;
    double totalObs6 = 0;
    double totalObs7 = 0;
    for (size_t i = 0; i < n; ++i) {
      totalReal += realMap[i][j];
      totalObs1 += probNEWMapAna->getCount(ids[j], i, 0);
      totalObs2 += probNEWMapTot->getCount(ids[j], i, 0);
      totalObs3 += VectorTools::sum(probNEWMapDet->getCounts(ids[j], i));
      totalObs4 += probNEWMapDecTot->getCount(ids[j], i, 0);
      totalObs5 += VectorTools::sum(probNEWMapDecDet->getCounts(ids[j], i));
      totalObs6 += probNEWMapUniTot->getCount(ids[j], i, 0);
      totalObs7 += VectorTools::sum(probNEWMapUniDet->getCounts(ids[j], i));
    }
    shared_ptr<PhyloNode> node=new_tree->getNode(ids[j]);
    
    cout << (new_tree->isLeaf(node)?node->getName():" ") << "\t";
    
    cout << new_tree->getEdgeToFather(node)->getLength() << "\t" << totalReal << "\t" << totalObs1 << "\t" << totalObs2 << "\t" << totalObs3 << "\t" << totalObs4 << "\t" << totalObs5 << "\t" << totalObs6 << "\t" << totalObs7 << endl;
    if (abs(totalReal - totalObs1) / totalReal > 0.1) throw Exception("Laplace substitution mapping failed, observed: " + TextTools::toString(totalObs1) + ", expected " + TextTools::toString(totalReal));
    //if (abs(totalReal - totalObs2) / totalReal > 0.1) return 1; //We do not expect the naive mapping to actually give an accurate result!
    //if (abs(totalReal - totalObs3) / totalReal > 0.1) return 1;
    if (abs(totalReal - totalObs4) / totalReal > 0.1) throw Exception("Uniformization (total) substitution mapping failed, observed: " + TextTools::toString(totalObs4) + ", expected " + TextTools::toString(totalReal));
    if (abs(totalReal - totalObs5) / totalReal > 0.1) throw Exception("Uniformization (detailed) substitution mapping failed, observed: " + TextTools::toString(totalObs5) + ", expected " + TextTools::toString(totalReal));
    if (abs(totalReal - totalObs6) / totalReal > 0.1) throw Exception("Decomposition (total) substitution mapping failed, observed: " + TextTools::toString(totalObs6) + ", expected " + TextTools::toString(totalReal));
    if (abs(totalReal - totalObs7) / totalReal > 0.1) throw Exception("Decomposition (detailed) substitution mapping failed, observed: " + TextTools::toString(totalObs7) + ", expected " + TextTools::toString(totalReal));
  }

  cout << endl;
  cout << "Details:" << endl;
  cout << "-------" << endl << endl;
  
//2. Detail:
  for (size_t j = 0; j < ids.size(); ++j) {
    vector<double> real(4, 0);
    vector<double> obs1(4, 0);
    vector<double> obs2(4, 0);
    vector<double> obs3(4, 0);
    for (size_t i = 0; i < n; ++i) {
      real += realMapDetailed[i][j];
      //VectorTools::print(real);
      vector<double> c = probNEWMapDet->getCounts(ids[j], i);
      //VectorTools::print(c);
      obs1 += probNEWMapDet->getCounts(ids[j], i);
      obs2 += probNEWMapDecDet->getCounts(ids[j], i);
      obs3 += probNEWMapUniDet->getCounts(ids[j], i);
    }
    shared_ptr<PhyloNode> node=new_tree->getNode(ids[j]);
    
    cout << (new_tree->isLeaf(node)?node->getName():" ") << "\t";
    
    cout << new_tree->getEdgeToFather(node)->getLength() << "\t";
    for (unsigned int t = 0; t < 4; ++t) {
      cout << obs1[t] << "/" << real[t] << "\t";
      cout << obs2[t] << "/" << real[t] << "\t";
      cout << obs3[t] << "/" << real[t] << "\t";
    }
    cout << endl;
    //if (abs(totalReal - totalObs) / totalReal > 0.1) return 1;
  }

  
  //-------------
  delete alphabet;
  delete sCountTot;
  delete sCountDet;
  delete probNEWMapTot;
  delete probNEWMapDet;
  delete probNEWMapDecTot;
  delete probNEWMapDecDet;
  delete probNEWMapUniTot;
  delete probNEWMapUniDet;
  //return (abs(obs - 0.001) < 0.001 ? 0 : 1);
  } catch (exception& e) {
    cout << "Test failed. Reason:" << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
