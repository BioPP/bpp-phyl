// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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

  shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET;
  shared_ptr<const NucleicAlphabet> nucAlphabet = AlphabetTools::DNA_ALPHABET;

  auto model = make_shared<GTR>(nucAlphabet, 1, 0.2, 0.3, 0.4, 0.4, 0.1, 0.35, 0.35, 0.2);
  auto rdist = make_shared<GammaDiscreteDistribution>(4, 0.4, 0.4);

  auto process = make_shared<RateAcrossSitesSubstitutionProcess>(model, shared_ptr<DiscreteDistributionInterface>(rdist->clone()), std::shared_ptr<PhyloTree>(new_tree->clone()));

  SimpleSubstitutionProcessSiteSimulator simulator(process);
  
  auto totReg = make_shared<TotalSubstitutionRegister>(model->getStateMap());
  auto detReg = make_shared<ComprehensiveSubstitutionRegister>(model->getStateMap());

  size_t n = 50000;
  vector< vector<double>> realMap(n);
  vector<vector<vector<double>>> realMapTotal(n);
  vector<vector<vector<double>>> realMapDetailed(n);
  auto sites = make_shared<VectorSiteContainer>(new_tree->getAllLeavesNames(), alphabet);
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
    auto simSite = result->getSite(*model);
    unique_ptr<Site> site(dynamic_cast<Site*>(simSite.release()));
    site->setCoordinate(static_cast<int>(i));    
    sites->addSite(site, false);
  }
  ApplicationTools::displayTaskDone();
cout << "ok0" << endl;
  //-------------
  //Now build the substitution vectors with the true model:

  // Newlik

  auto tmComp = make_shared<LikelihoodCalculationSingleProcess>(context, sites, process);
cout << "ok1" << endl;
  SingleProcessPhyloLikelihood newTl(context, tmComp);
cout << "ok2" << endl;
  cout << "LogLik: " << newTl.getValue() << endl;
    
  
  
  auto sCountAna = make_shared<LaplaceSubstitutionCount>(model, 10);
  auto m = sCountAna->getAllNumbersOfSubstitutions(0.001, 1);
  cout << "Analytical (Laplace) total count:" << endl;
  MatrixTools::print(*m);
  auto probNEWMapAna = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountAna);

  cout << endl;

  //Simple:
  auto sCountTot = make_shared<NaiveSubstitutionCount>(model, totReg);
  m = sCountTot->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Simple total count:" << endl;
  MatrixTools::print(*m);
  auto probNEWMapTot = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountTot);
  cout << endl;

  auto sCountDet = make_shared<NaiveSubstitutionCount>(model, detReg);
  m = sCountDet->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Detailed count, type 1:" << endl;
  MatrixTools::print(*m);
  auto probNEWMapDet = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountDet);
  cout << endl;

  //Decomposition:
  auto sCountDecTot = make_shared<DecompositionSubstitutionCount>(model, totReg);
  m = sCountDecTot->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Total count, decomposition method:" << endl;
  MatrixTools::print(*m);
  auto probNEWMapDecTot = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountDecTot);

  auto sCountDecDet = make_shared<DecompositionSubstitutionCount>(model, detReg);
  m = sCountDecDet->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Detailed count, decomposition method, type 1:" << endl;
  MatrixTools::print(*m);
  auto probNEWMapDecDet = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountDecDet);
  cout << endl;

  //Uniformization
  auto sCountUniTot = make_shared<UniformizationSubstitutionCount>(model, totReg);
  m = sCountUniTot->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Total count, uniformization method:" << endl;
  MatrixTools::print(*m);
  auto probNEWMapUniTot = 
    SubstitutionMappingTools::computeCounts(*tmComp, *sCountUniTot);  

  auto sCountUniDet = make_shared<UniformizationSubstitutionCount>(model, detReg);
  m = sCountUniDet->getAllNumbersOfSubstitutions(0.001,1);
  cout << "Detailed count, uniformization method, type 1:" << endl;
  MatrixTools::print(*m);
  auto probNEWMapUniDet = 
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
  } catch (exception& e) {
    cout << "Test failed. Reason:" << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
