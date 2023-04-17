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
#include <Bpp/Phyl/Legacy/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSequenceSimulator.h>

#include <Bpp/Phyl/Legacy/Likelihood/RNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Legacy/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Legacy/OptimizationTools.h>

#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>

#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

#include <iostream>

using namespace bpp;
using namespace std;

int main() {

  auto tree = TreeTemplateTools::parenthesisToTree("(((A:0.1, B:0.2):0.3,C:0.15):0.25,(D:0.35,(E:0.26,F:0.05):0.12):0.16);");

  Newick reader;
  shared_ptr<PhyloTree> pTree = reader.parenthesisToPhyloTree("(((A:0.1, B:0.2):0.3,C:0.15):0.25,(D:0.35,(E:0.26,F:0.05):0.12):0.16);", false, "", false, false);

  vector<string> seqNames = tree->getLeavesNames();
  vector<int> ids = tree->getNodesId();
  //-------------

  shared_ptr<const NucleicAlphabet> nucAlphabet = AlphabetTools::DNA_ALPHABET;
  shared_ptr<const Alphabet> alphabet = AlphabetTools::DNA_ALPHABET;
  
  auto rootFreqs = make_shared<GCFrequencySet>(nucAlphabet);
  auto model = make_shared<T92>(nucAlphabet, 3., .1);
  map<string, vector<Vint>> globalParameterVectors;
  globalParameterVectors["T92.kappa"] = vector<Vint>();
  
  //Very difficult to optimize on small datasets:
  auto rdist = make_shared<GammaDiscreteRateDistribution>(4, 1.0);
  
  auto rootFreqs2 = shared_ptr<FrequencySetInterface>(rootFreqs->clone());

  auto rdist2 = shared_ptr<DiscreteDistribution>(rdist->clone());
  std::shared_ptr<SubstitutionModelInterface> model2(model->clone());

  map<string, string> alias;

  shared_ptr<SubstitutionModelSet> modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(shared_ptr<SubstitutionModelInterface>(model->clone()), rootFreqs, *tree, alias, globalParameterVectors);

  std::vector<std::string> globalParameterNames;
  globalParameterNames.push_back("T92.kappa");

  shared_ptr<SubstitutionProcessInterface> subProSim = NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(model2, rdist2, pTree, rootFreqs2, globalParameterNames);


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

  SimpleSubstitutionProcessSequenceSimulator simulator(subProSim);

  nrep = 20;
  
  for (unsigned int j = 0; j < nrep; ++j) {

    auto profiler  = make_shared<StlOutputStream>(make_unique<ofstream>("profile.txt", ios::out));
    auto messenger = make_shared<StlOutputStream>(make_unique<ofstream>("messages.txt", ios::out));

    //Simulate data:
    shared_ptr<SiteContainerInterface> sites = simulator.simulate(nsites);

    //Now fit model:

    auto tl = make_shared<RNonHomogeneousTreeLikelihood>(*tree, *sites, modelSet, rdist, true, true, false);
    tl->initialize();

    Context context;
    auto lik = make_shared<LikelihoodCalculationSingleProcess>(context, sites, subProSim);
    auto ntl = make_shared<SingleProcessPhyloLikelihood>(context, lik);

    cout << setprecision(10) << "OldTL init: "  << tl->getValue()  << endl;
    cout << setprecision(10) << "NewTL init: "  << ntl->getValue()  << endl;

    unsigned int c1 = OptimizationToolsOld::optimizeNumericalParameters2(
      tl, tl->getParameters(), 0,
      0.0001, 10000,
      messenger, profiler,
      false, false,
      1, OptimizationTools::OPTIMIZATION_NEWTON);
    
    unsigned int nc1 = OptimizationTools::optimizeNumericalParameters2(
      ntl, ntl->getParameters(), 0,
      0.0001, 10000,
      messenger, profiler,
      false, false,
      1, OptimizationTools::OPTIMIZATION_NEWTON);


    cout << "OldTL optim: " << c1 << ": " << tl->getValue()  << endl;
    cout << "NewTL optim: " << nc1 << ": " << ntl->getValue() << endl;

    cout << "Thetas : " << endl;
    
    for (size_t i = 0; i < nmodels; ++i) {
      cout << tl->substitutionModelSet().model(i).getParameter("theta").getValue() << "\t" << ntl->getLikelihoodCalculation()->getParameter("T92.theta_"+to_string(i+1)).getValue() << endl;
      //if (abs(modelSet2->getModel(i)->getParameter("theta").getValue() - modelSet3->getModel(i)->getParameter("theta").getValue()) > 0.1)
      //  return 1;
      thetasEst1[i] += tl->substitutionModelSet().model(i).getParameter("theta").getValue();
      thetasEst1n[i] += ntl->likelihoodCalculation().getParameter("T92.theta_"+to_string(i+1)).getValue();
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
