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
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSequenceSimulator.h>

#include <Bpp/Phyl/Likelihood/RNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>

#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/RateAcrossSitesSubstitutionProcess.h>

#include <Bpp/Phyl/NewLikelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

#include <iostream>

using namespace bpp;
using namespace std;

int main() {

  TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree("(((A:0.1, B:0.2):0.3,C:0.15):0.25,(D:0.35,(E:0.26,F:0.05):0.12):0.16);");

  Newick reader;
  unique_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((A:0.1, B:0.2):0.3,C:0.15):0.25,(D:0.35,(E:0.26,F:0.05):0.12):0.16);", false, "", false, false));
  ParametrizablePhyloTree parTree(*pTree);

  vector<string> seqNames= tree->getLeavesNames();
  vector<int> ids = tree->getNodesId();
  //-------------

  const NucleicAlphabet* alphabet = &AlphabetTools::DNA_ALPHABET;
  
  auto rootFreqs = std::make_shared<GCFrequencySet>(alphabet);
  auto model = std::make_shared<T92>(alphabet, 3.);
  std::map<std::string, std::vector<Vint>> globalParameterVectors;
  globalParameterVectors["T92.kappa"]=std::vector<Vint>();
  
  //Very difficult to optimize on small datasets:
  DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(4, 1.0);
  
  auto rootFreqs2 = std::shared_ptr<FrequencySet>(dynamic_cast<FrequencySet*>(rootFreqs->clone()));
  DiscreteDistribution* rdist2 = rdist->clone();
  std::shared_ptr<SubstitutionModel> model2(model->clone());

  map<string, string> alias;

  SubstitutionModelSet* modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model->clone(), rootFreqs, tree, alias, globalParameterVectors);

  std::vector<std::string> globalParameterNames;
  globalParameterNames.push_back("T92.kappa");

  NonHomogeneousSubstitutionProcess* subProSim= NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(model2, rdist2, parTree.clone(), rootFreqs2, globalParameterNames);

  SubstitutionProcess* nsubPro=subProSim->clone();

  // Simulation
  size_t nsites = 1000;
  unsigned int nrep = 3;
  size_t nmodels = modelSet->getNumberOfModels();
  vector<double> thetas(nmodels);
  vector<double> thetasEst1(nmodels);
  vector<double> thetasEst1n(nmodels);

  for (size_t i = 0; i < nmodels; ++i) {
    double theta = RandomTools::giveRandomNumberBetweenZeroAndEntry(0.9) + 0.05;
    cout << "Theta" << i << " set to " << theta << endl; 
    subProSim->setParameterValue("T92.theta_" + TextTools::toString(i + 1), theta);
    thetas[i] = theta;
  }

  SimpleSubstitutionProcessSequenceSimulator simulator(*subProSim);

  nrep=20;
  
  for (unsigned int j = 0; j < nrep; j++) {

    OutputStream* profiler  = new StlOutputStream(new ofstream("profile.txt", ios::out));
    OutputStream* messenger = new StlOutputStream(new ofstream("messages.txt", ios::out));

    //Simulate data:
    auto sites(simulator.simulate(nsites));

    //Now fit model:

    RNonHomogeneousTreeLikelihood tl(*tree, *sites.get(), modelSet, rdist, true, true, false);
    tl.initialize();

    Context context;
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, *sites->clone(), *nsubPro);
    
    SingleProcessPhyloLikelihood ntl(context, lik, lik->getParameters());

    cout << setprecision(10) << "OldTL init: "  << tl.getValue()  << endl;
    cout << setprecision(10) << "NewTL init: "  << ntl.getValue()  << endl;

    unsigned int c1 = OptimizationTools::optimizeNumericalParameters2(
      &tl, tl.getParameters(), 0,
      0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);
    
    unsigned int nc1 = OptimizationTools::optimizeNumericalParameters2(
      ntl, ntl.getParameters(), 0,
      0.0001, 10000, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);


    cout << "OldTL optim: " << c1 << ": " << tl.getValue()  << endl;
    cout << "NewTL optim: " << nc1 << ": " << ntl.getValue() << endl;

    cout << "Thetas : " << endl;
    
    for (size_t i = 0; i < nmodels; ++i) {
      cout << tl.getSubstitutionModelSet()->getModel((int)i)->getParameter("theta").getValue() << "\t" << ntl.getLikelihoodCalculation()->getParameter("T92.theta_"+to_string(i+1)).getValue() << endl;
      //if (abs(modelSet2->getModel(i)->getParameter("theta").getValue() - modelSet3->getModel(i)->getParameter("theta").getValue()) > 0.1)
      //  return 1;
      thetasEst1[i] += tl.getSubstitutionModelSet()->getModel((int)i)->getParameter("theta").getValue();
      thetasEst1n[i] += ntl.getLikelihoodCalculation()->getParameter("T92.theta_"+to_string(i+1)).getValue();
    }
  }
  thetasEst1 /= static_cast<double>(nrep);
  thetasEst1n /= static_cast<double>(nrep);

  //Now compare estimated values to real ones:
  cout << "Real" << "\t" << "Est_Old1" << "\t";
  cout << "Est_New1" <<  endl;
  for (size_t i = 0; i < thetas.size(); ++i) {
    cout << thetas[i] << "\t" << thetasEst1[i] << "\t";
    cout << thetasEst1n[i] << endl;
     double diff1 = abs(thetas[i] - thetasEst1[i]);
     double diffn1 = abs(thetas[i] - thetasEst1n[i]);
     if (diff1 > 0.2  || diffn1 > 0.2 )
       return 1;
  }

  return 0;
}
