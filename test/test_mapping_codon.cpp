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
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Codon/YN98.h>
#include <Bpp/Phyl/Model/FrequenciesSet/CodonFrequenciesSet.h>
#include <Bpp/Phyl/Simulation/SubstitutionProcessSequenceSimulator.h>
#include <Bpp/Phyl/Mapping/SubstitutionRegister.h>
#include <Bpp/Phyl/Mapping/NaiveSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/LaplaceSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/UniformizationSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/DecompositionSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/SubstitutionMappingTools.h>
#include <Bpp/Phyl/Mapping/ProbabilisticSubstitutionMapping.h>
#include <Bpp/Phyl/Mapping/SubstitutionMappingTools.h>
#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/NewLikelihood/SimpleSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/RateAcrossSitesSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  Newick reader;
  Context context;
  
  unique_ptr<PhyloTree> new_tree(reader.parenthesisToPhyloTree("((A:0.001, B:0.002):0.008,C:0.01,D:0.1);", false, "", false, false));

  vector<int> ids = {0, 1, 2, 3, 4};

  //-------------

  auto gc = std::make_shared<StandardGeneticCode>(AlphabetTools::DNA_ALPHABET);
  
  auto model = std::make_shared<YN98>(gc.get(), CodonFrequenciesSet::getFrequenciesSetForCodons(CodonFrequenciesSet::F0, gc.get()));
  DiscreteDistribution* rdist = new ConstantDistribution(1.0);
  std::shared_ptr<ParametrizablePhyloTree> pTree(new ParametrizablePhyloTree(*new_tree));
  unique_ptr<RateAcrossSitesSubstitutionProcess> process(new RateAcrossSitesSubstitutionProcess(model, rdist->clone(), pTree->clone()));

  SimpleSubstitutionProcessSequenceSimulator simulator(*process);

  TotalSubstitutionRegister* totReg = new TotalSubstitutionRegister(model->getStateMap());
  DnDsSubstitutionRegister* dndsReg = new DnDsSubstitutionRegister(model->getStateMap(), *gc);

  unsigned int n = 10000;
  vector< vector<double> > realMap(n);
  vector< vector< vector<double> > > realMapTotal(n);
  vector< vector< vector<double> > > realMapDnDs(n);
  VectorSiteContainer sites(new_tree->getAllLeavesNames(), gc->getSourceAlphabet());
  for (unsigned int i = 0; i < n; ++i) {
    ApplicationTools::displayGauge(i, n-1, '=');
    unique_ptr<New_SiteSimulationResult> result(simulator.dSimulateSite());
    realMap[i].resize(ids.size());
    realMapTotal[i].resize(ids.size());
    realMapDnDs[i].resize(ids.size());
    for (size_t j = 0; j < ids.size(); ++j) {
      realMap[i][j] = static_cast<double>(result->getSubstitutionCount((uint)ids[j]));
      realMapTotal[i][j].resize(totReg->getNumberOfSubstitutionTypes());
      realMapDnDs[i][j].resize(dndsReg->getNumberOfSubstitutionTypes());
      result->getSubstitutionCount((uint)ids[j], *totReg, realMapTotal[i][j]);
      result->getSubstitutionCount((uint)ids[j], *dndsReg, realMapDnDs[i][j]);
      if (realMapTotal[i][j][0] != realMap[i][j]) {
        cerr << "Error, total substitution register provides wrong result." << endl;
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
  // Matrix<double>* m = sCountAna->getAllNumbersOfSubstitutions(0.001,1);
  // cout << "Analytical total count:" << endl;
  // MatrixTools::print(*m);
  // delete m;
  ProbabilisticSubstitutionMapping* probMapAna = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountAna);

  SubstitutionCount* sCountTot = new NaiveSubstitutionCount(model.get(), totReg);
  // m = sCountTot->getAllNumbersOfSubstitutions(0.001,1);
  // cout << "Simple total count:" << endl;
  // MatrixTools::print(*m);
  // delete m;
  ProbabilisticSubstitutionMapping* probMapTot = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountTot);

  SubstitutionCount* sCountDnDs = new NaiveSubstitutionCount(model.get(), dndsReg);
  // m = sCountDnDs->getAllNumbersOfSubstitutions(0.001,1);
  // cout << "Detailed count, type 1:" << endl;
  // MatrixTools::print(*m);
  // delete m;
  ProbabilisticSubstitutionMapping* probMapDnDs = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountDnDs);
  
  SubstitutionCount* sCountUniTot = new UniformizationSubstitutionCount(model.get(), totReg);
  // m = sCountUniTot->getAllNumbersOfSubstitutions(0.001,1);
  // cout << "Total count, uniformization method:" << endl;
  // MatrixTools::print(*m);
  // delete m;
  ProbabilisticSubstitutionMapping* probMapUniTot = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountUniTot);

  SubstitutionCount* sCountUniDnDs = new UniformizationSubstitutionCount(model.get(), dndsReg);
  // m = sCountUniDnDs->getAllNumbersOfSubstitutions(0.001,2);
  // cout << "Detailed count, uniformization method, type 2:" << endl;
  // MatrixTools::print(*m);
  // delete m;
  ProbabilisticSubstitutionMapping* probMapUniDnDs = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountUniDnDs);

  //Check per branch:
  
  //1. Total:
  for (unsigned int j = 0; j < ids.size(); ++j) {
    double totalReal = 0;
    double totalDnDs = 0;
    double totalObs1 = 0;
    double totalObs2 = 0;
    double totalObs3 = 0;
    double totalObs4 = 0;
    double totalObs5 = 0;
    for (unsigned int i = 0; i < n; ++i) {
      totalReal += realMap[i][j];
      totalDnDs += VectorTools::sum(realMapDnDs[i][j]);
      totalObs1 += probMapAna->getCount(ids[j], i, 0);
      totalObs2 += probMapTot->getCount(ids[j], i, 0);
      totalObs3 += VectorTools::sum(probMapDnDs->getCounts(ids[j], i));
      totalObs4 += probMapUniTot->getCount(ids[j], i, 0);
      totalObs5 += VectorTools::sum(probMapUniDnDs->getCounts(ids[j], i));
    }
    shared_ptr<PhyloNode> node=new_tree->getNode(ids[j]);
    
    cout << (new_tree->isLeaf(node)?node->getName():" ") << "\t";
    
    cout << new_tree->getEdgeToFather(node)->getLength() << "\t";
    cout << totalReal << "\t" << totalObs1 << "\t" << totalObs2 
         << "\t" << totalObs4 << "\t\t" << totalDnDs<< "\t" << totalObs3 << "\t" << totalObs5
         << endl;
    if (abs(totalReal - totalObs1) / totalReal > 0.1) return 1;
    if (abs(totalReal - totalObs2) / totalReal > 0.1) return 1;
    if (abs(totalDnDs - totalObs3) / totalReal > 0.1) return 1;
    if (abs(totalReal - totalObs4) / totalReal > 0.1) return 1;
    if (abs(totalDnDs - totalObs5) / totalReal > 0.1) return 1;
  }

  //-------------
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
