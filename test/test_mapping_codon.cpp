// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
#include <Bpp/Phyl/Model/FrequencySet/CodonFrequencySet.h>
#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSiteSimulator.h>
#include <Bpp/Phyl/Mapping/SubstitutionRegister.h>
#include <Bpp/Phyl/Mapping/NaiveSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/LaplaceSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/UniformizationSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/DecompositionSubstitutionCount.h>
#include <Bpp/Phyl/Mapping/SubstitutionMappingTools.h>
#include <Bpp/Phyl/Mapping/ProbabilisticSubstitutionMapping.h>
#include <Bpp/Phyl/Mapping/SubstitutionMappingTools.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/SimpleSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  Newick reader;
  Context context;
  
  shared_ptr<PhyloTree> new_tree(reader.parenthesisToPhyloTree("((A:0.001, B:0.002):0.008,C:0.01,D:0.1);", false, "", false, false));

  vector<uint> ids = {0, 1, 2, 3, 4};

  //-------------

  auto gc = std::make_shared<StandardGeneticCode>(AlphabetTools::DNA_ALPHABET);
  
  auto model = std::make_shared<YN98>(gc, CodonFrequencySetInterface::getFrequencySetForCodons(CodonFrequencySetInterface::F0, gc));
  auto rdist = std::make_shared<ConstantDistribution>(1.0);
  shared_ptr<RateAcrossSitesSubstitutionProcess> process(new RateAcrossSitesSubstitutionProcess(model, rdist, std::shared_ptr<PhyloTree>(new_tree->clone())));

  SimpleSubstitutionProcessSiteSimulator simulator(process);

  auto totReg = make_shared<TotalSubstitutionRegister>(model->getStateMap());
  auto dndsReg = make_shared<DnDsSubstitutionRegister>(model->getStateMap(), gc);

  unsigned int n = 10000;
  vector< vector<double> > realMap(n);
  vector< vector< vector<double> > > realMapTotal(n);
  vector< vector< vector<double> > > realMapDnDs(n);
  auto sites = make_shared<VectorSiteContainer>(new_tree->getAllLeavesNames(), gc->getSourceAlphabet());
  for (unsigned int i = 0; i < n; ++i) {
    ApplicationTools::displayGauge(i, n-1, '=');
    unique_ptr<SiteSimulationResult> result(simulator.dSimulateSite());
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
    auto simSite = result->getSite(*model);
    unique_ptr<Site> site(dynamic_cast<Site*>(simSite.release()));
    site->setCoordinate(static_cast<int>(i));
    sites->addSite(site, false);
  }
  ApplicationTools::displayTaskDone();
  
  //-------------
  //Now build the substitution vectors with the true model:

  // Newlik
  auto tmComp = make_shared<LikelihoodCalculationSingleProcess>(context, sites, process);
  SingleProcessPhyloLikelihood newTl(context, tmComp);
  cout << "LogLik: " << newTl.getValue() << endl;
    

  auto sCountAna = make_shared<LaplaceSubstitutionCount>(model, 10);
  // Matrix<double>* m = sCountAna->getAllNumbersOfSubstitutions(0.001,1);
  // cout << "Analytical total count:" << endl;
  // MatrixTools::print(*m);
  // delete m;
  auto probMapAna = SubstitutionMappingTools::computeCounts(*tmComp, *sCountAna);

  auto sCountTot = make_shared<NaiveSubstitutionCount>(model, totReg);
  // m = sCountTot->getAllNumbersOfSubstitutions(0.001,1);
  // cout << "Simple total count:" << endl;
  // MatrixTools::print(*m);
  // delete m;
  auto probMapTot = SubstitutionMappingTools::computeCounts(*tmComp, *sCountTot);

  auto sCountDnDs = make_shared<NaiveSubstitutionCount>(model, dndsReg);
  // m = sCountDnDs->getAllNumbersOfSubstitutions(0.001,1);
  // cout << "Detailed count, type 1:" << endl;
  // MatrixTools::print(*m);
  // delete m;
  auto probMapDnDs = SubstitutionMappingTools::computeCounts(*tmComp, *sCountDnDs);
  
  auto sCountUniTot = make_shared<UniformizationSubstitutionCount>(model, totReg);
  // m = sCountUniTot->getAllNumbersOfSubstitutions(0.001,1);
  // cout << "Total count, uniformization method:" << endl;
  // MatrixTools::print(*m);
  // delete m;
  auto probMapUniTot = SubstitutionMappingTools::computeCounts(*tmComp, *sCountUniTot);

  auto sCountUniDnDs = make_shared<UniformizationSubstitutionCount>(model, dndsReg);
  // m = sCountUniDnDs->getAllNumbersOfSubstitutions(0.001,2);
  // cout << "Detailed count, uniformization method, type 2:" << endl;
  // MatrixTools::print(*m);
  // delete m;
  auto probMapUniDnDs = SubstitutionMappingTools::computeCounts(*tmComp, *sCountUniDnDs);

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
    
    cout << "\t\t\t";

    cout << 
      (abs(totalReal - totalObs1) / totalReal)  << "\t"
         << (abs(totalReal - totalObs2) / totalReal ) << "\t"
         << (abs(totalReal - totalObs4) / totalReal ) << "\t"
         << (abs(totalDnDs - totalObs3) / totalDnDs ) << "\t"
         << (abs(totalDnDs - totalObs5) / totalDnDs ) << endl;

    if (abs(totalReal - totalObs1) / totalReal > 0.1) return 1;
    if (abs(totalReal - totalObs2) / totalReal > 0.1) return 1;
    if (abs(totalReal - totalObs4) / totalReal > 0.1) return 1;
    if (abs(totalDnDs - totalObs3) / totalDnDs > 0.1) return 1;
    if (abs(totalDnDs - totalObs5) / totalDnDs > 0.1) return 1;
  }

  //-------------
  //return (abs(obs - 0.001) < 0.001 ? 0 : 1);
  return 0;
}
