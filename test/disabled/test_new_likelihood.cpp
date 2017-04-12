// File: test_likelihood.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 23/02/2017

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

#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>

#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/RateAcrossSitesSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/SimpleSubstitutionProcess.h>

#include <Bpp/Utils/Cpp14.h>
#include <Bpp/Utils/ForRange.h>
#include <chrono>
#include <iostream>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

using bpp::Cpp14::make_unique;
using bpp::makeRange;

namespace
{
  using TimePoint = typename std::chrono::high_resolution_clock::time_point;
  TimePoint timingStart(void) { return std::chrono::high_resolution_clock::now(); }
  void timingEnd(TimePoint start, const std::string& prefix)
  {
    auto end = timingStart(); // ill named, just to get now()
    std::cout << "[time-ns] " << prefix << " "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << "\n";
  }

  template <typename Lik>
  void do_param_changes_multiple_times(Lik& llh, const std::string& timePrefix, const bpp::ParameterList& p1,
                                       const bpp::ParameterList& p2)
  {
    constexpr std::size_t updatesNbIterations = 1000;
    auto ts = timingStart();
    for (auto i : makeRange(updatesNbIterations))
    {
      (void)i;
      llh.matchParametersValues(p1);
      llh.getValue();
      llh.matchParametersValues(p2);
      llh.getValue();
    }
    timingEnd(ts, timePrefix);
  }

  template <typename Lik>
  void do_optimization(Lik& llh, const std::string& timePrefix)
  {
    constexpr int nbOptim = 100000;
    auto ts = timingStart();
    bpp::OptimizationTools::optimizeNumericalParameters2(&llh, llh.getParameters(), 0, 0.000001, nbOptim, 0, 0);
    timingEnd(ts, timePrefix);
  }
}

struct ValuesToCompare
{
  double initialLikelihood{};
  double initial1DerivativeBr2{};
  double initial2DerivativeBr2{};
  double finalLikelihood{};
  ValuesToCompare() = default;
};

TEST_CASE("comparing results between old and new likelihood (single traversal)")
{
  using namespace bpp;

  // Input sequences
  const auto& alphabet = AlphabetTools::DNA_ALPHABET;
  VectorSiteContainer sites(&alphabet);
  sites.addSequence(
    BasicSequence("A", "ATCCAGACATGCCGGGACTTTGCAGAGAAGGAGTTGTTTCCCATTGCAGCCCAGGTGGATAAGGAACAGC", &alphabet));
  sites.addSequence(
    BasicSequence("B", "CGTCAGACATGCCGTGACTTTGCCGAGAAGGAGTTGGTCCCCATTGCGGCCCAGCTGGACAGGGAGCATC", &alphabet));
  sites.addSequence(
    BasicSequence("C", "GGTCAGACATGCCGGGAATTTGCTGAAAAGGAGCTGGTTCCCATTGCAGCCCAGGTAGACAAGGAGCATC", &alphabet));
  sites.addSequence(
    BasicSequence("D", "TTCCAGACATGCCGGGACTTTACCGAGAAGGAGTTGTTTTCCATTGCAGCCCAGGTGGATAAGGAACATC", &alphabet));

  // Evolution model
  T92 model(&alphabet, 3.);
  GammaDiscreteRateDistribution rateDistribution(4, 1.0);

  // Set of parameters to apply to tree + model
  ParameterList paramModel1;
  paramModel1.addParameter(Parameter("T92.kappa", 0.1));
  ParameterList paramModel2;
  paramModel2.addParameter(Parameter("T92.kappa", 0.2));

  ParameterList paramBrLen1;
  paramBrLen1.addParameter(Parameter("BrLen1", 0.1));
  ParameterList paramBrLen2;
  paramBrLen2.addParameter(Parameter("BrLen1", 0.2));

  // FIXME check values ???
  // double initialValue = 228.6333642493463;
  // double finalValue = 198.47216106233;

  // Old likelihood
  ValuesToCompare oldL;
  if (0){
    auto ts = timingStart();
    auto tree = std::unique_ptr<TreeTemplate<Node>>(
      TreeTemplateTools::parenthesisToTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);"));
    RHomogeneousTreeLikelihood llh(*tree, sites, model.clone(), rateDistribution.clone(), false, false);
    llh.initialize();
    oldL.initialLikelihood = llh.getValue();
    timingEnd(ts, "old_init_value");

    oldL.initial1DerivativeBr2 = llh.getFirstOrderDerivative("BrLen2");
    oldL.initial2DerivativeBr2 = llh.getSecondOrderDerivative("BrLen2");

    do_param_changes_multiple_times(llh, "old_param_model_change", paramModel1, paramModel2);
    do_param_changes_multiple_times(llh, "old_param_brlen_change", paramBrLen1, paramBrLen2);
    do_optimization(llh, "old_optimization");
    oldL.finalLikelihood = llh.getValue();
  }

  // New likelihood
  ValuesToCompare newL;
  {
    auto ts = timingStart();
    Newick reader;
    auto phyloTree = std::unique_ptr<PhyloTree>(
      reader.parenthesisToPhyloTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);", false, "", false, false));
    auto paramPhyloTree = std::make_shared<ParametrizablePhyloTree>(*phyloTree);
    auto process =
      make_unique<RateAcrossSitesSubstitutionProcess>(model.clone(), rateDistribution.clone(), paramPhyloTree->clone());
    auto likelihoodCompStruct = make_unique<RecursiveLikelihoodTreeCalculation>(sites, process.get(), false, true);
    SingleProcessPhyloLikelihood llh(process.get(), likelihoodCompStruct.release());
    llh.computeLikelihood();
    newL.initialLikelihood = llh.getValue();
    timingEnd(ts, "new_init_value");

    newL.initial1DerivativeBr2 = llh.getFirstOrderDerivative("BrLen2");
    newL.initial2DerivativeBr2 = llh.getSecondOrderDerivative("BrLen2");

    do_param_changes_multiple_times(llh, "new_param_model_change", paramModel1, paramModel2);
    do_param_changes_multiple_times(llh, "new_param_brlen_change", paramBrLen1, paramBrLen2);
    do_optimization(llh, "new_optimization");
    newL.finalLikelihood = llh.getValue();
  }

  // TODO newTlop.getParameters().printParameters(cout);

  CHECK(doctest::Approx(oldL.initialLikelihood) == newL.initialLikelihood);
  CHECK(doctest::Approx(oldL.initial1DerivativeBr2) == newL.initial1DerivativeBr2);
  CHECK(doctest::Approx(oldL.initial2DerivativeBr2) == newL.initial2DerivativeBr2);
  CHECK(doctest::Approx(oldL.finalLikelihood).epsilon(0.0001) == newL.finalLikelihood);
}
