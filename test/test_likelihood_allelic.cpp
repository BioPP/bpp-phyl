//
// File: test_likelihood.cpp
// Created by: Julien Dutheil
// Created on: Mon Apr 04 10:18 2011
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

#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/POMO.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequencySet/FrequencySet.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorProbabilisticSiteContainer.h>
#include <Bpp/Phyl/Legacy/OptimizationTools.h>

#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/SimpleSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>

#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

#include <iostream>


using namespace bpp;
using namespace std;


int main() {
  
  Newick reader;
  auto phyloTree = std::shared_ptr<PhyloTree>(reader.readPhyloTree("lysozymeLarge.dnd"));

  //-------------

  string nameSeq = "lysozymeLarge.fasta";
  Fasta fasta;
 
  auto alphaall = std::make_shared<AllelicAlphabet>(AlphabetTools::DNA_ALPHABET, 3);
  
  auto sites = std::shared_ptr<AlignedSequenceContainer>(fasta.readAlignment(nameSeq , &alphaall->getStateAlphabet()));
  
  SiteContainerTools::changeGapsToUnknownCharacters(*sites);

  vector<std::shared_ptr<ProbabilisticSequence>> vseq;
  for (size_t ns=0;ns < sites->getNumberOfSequences(); ns++)
    vseq.push_back(std::shared_ptr<ProbabilisticSequence>(alphaall->convertFromStateAlphabet(sites->getSequence(ns))));

  VectorProbabilisticSiteContainer sites2(alphaall.get());

  for (const auto& seq:vseq)
    sites2.addSequence(*seq);

  auto t92 = std::make_shared<T92>(&AlphabetTools::DNA_ALPHABET);
  
  auto fitness = std::make_shared<FullNucleotideFrequencySet>(&AlphabetTools::DNA_ALPHABET);

  auto model = std::make_shared<POMO>(alphaall.get(), t92, fitness);

  auto rootFreqs = std::make_shared<FullFrequencySet>(model->shareStateMap());

  auto distribution = std::make_shared<ConstantRateDistribution>();
  
  auto process  = NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model, distribution, phyloTree, rootFreqs);

//  auto process= std::make_shared<SimpleSubstitutionProcess>(model, phyloTree);

  Context context;                        

  auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, sites2, *process);

  SingleProcessPhyloLikelihood llh(context, lik);

  
  cout << "NewTL: " << setprecision(20) << llh.getValue() << endl;

  unique_ptr<OutputStream> messenger(new StlOutputStream(new ofstream("messages.txt", ios::out)));
  unique_ptr<OutputStream> profiler(new StlOutputStream(new ofstream("profile.txt", ios::out)));
  profiler->setPrecision(20);

  OptimizationTools::optimizeNumericalParameters2(llh, llh.getParameters(), 0, 0.000001, 10000, messenger.get(), profiler.get(), false, true, 2, OptimizationTools::OPTIMIZATION_NEWTON);

  llh.getParameters().printParameters(cerr);

  cout << setprecision(20) << llh.getValue() << endl;

  return 0;
}


