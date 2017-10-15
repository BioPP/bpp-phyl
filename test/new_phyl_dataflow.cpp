//
// File: new_phyl_dataflow.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-19 00:00:00
// Last modified: 2017-06-08
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

  This software is a computer program whose purpose is to provide classes
  for numerical calculus. This file is part of the Bio++ project.

  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#define ENABLE_OLD
#define ENABLE_NEW
//#define ENABLE_DF

// Common stuff
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Eigen/Dense>
#include <chrono>
#include <fstream>

// Old likelihood
#ifdef ENABLE_OLD
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#endif
// Newlik
#ifdef ENABLE_NEW
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/SimpleSubstitutionProcess.h>
#endif
// DF
#ifdef ENABLE_DF
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/ImportMaster.h>
#include <Bpp/NewPhyl/ImportNewlik.h>
#include <Bpp/NewPhyl/Model.h>
#include <Bpp/NewPhyl/Phylogeny.h>
#endif

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
  void printLik(double logLik, const std::string& prefix)
  {
    std::cout << "[log-lik] " << prefix << " " << logLik << "\n";
  }

  template<typename Func>
  void do_func_multiple_times(const std::string& timePrefix, Func f)
  {
    constexpr std::size_t updatesNbIterations = 100000;
    auto ts = timingStart();
    for (auto i : bpp::range(updatesNbIterations))
    {
      (void)i;
      f();
    }
    timingEnd(ts, timePrefix);
  }
  template<typename Lik>
  void do_param_changes_multiple_times_legacy(Lik& llh,
                                              const std::string& timePrefix,
                                              const bpp::ParameterList& p1,
                                              const bpp::ParameterList& p2)
  {
    llh.matchParametersValues(p1);
    printLik(llh.getValue(), timePrefix);
    llh.matchParametersValues(p2);
    printLik(llh.getValue(), timePrefix);

    do_func_multiple_times(timePrefix, [&]() {
      llh.matchParametersValues(p1);
      llh.getValue();
      llh.matchParametersValues(p2);
      llh.getValue();
    });
  }
#ifdef ENABLE_DF
  void do_param_changes_multiple_times_df(const bpp::DF::ValueRef<double>& lik,
                                          const std::string& timePrefix,
                                          bpp::DF::ParameterRef<double> param,
                                          double v1,
                                          double v2)
  {
    param->setValue(v1);
    printLik(lik->getValue(), timePrefix);
    param->setValue(v2);
    printLik(lik->getValue(), timePrefix);

    do_func_multiple_times(timePrefix, [&]() {
      param->setValue(v1);
      lik->getValue();
      param->setValue(v2);
      lik->getValue();
    });
  }
#endif

  struct CommonStuff
  {
    const bpp::NucleicAlphabet& alphabet;
    bpp::VectorSiteContainer sites;
    const char* treeStr;
    bpp::ParameterList paramModel1;
    bpp::ParameterList paramModel2;
    bpp::ParameterList paramBrLen1;
    bpp::ParameterList paramBrLen2;

    CommonStuff()
      : alphabet(bpp::AlphabetTools::DNA_ALPHABET)
      , sites(&alphabet)
      , treeStr("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);")
    {
      // Init sequences
      sites.addSequence(
        bpp::BasicSequence("A", "ATCCAGACATGCCGGGACTTTGCAGAGAAGGAGTTGTTTCCCATTGCAGCCCAGGTGGATAAGGAACAGC", &alphabet));
      sites.addSequence(
        bpp::BasicSequence("B", "CGTCAGACATGCCGTGACTTTGCCGAGAAGGAGTTGGTCCCCATTGCGGCCCAGCTGGACAGGGAGCATC", &alphabet));
      sites.addSequence(
        bpp::BasicSequence("C", "GGTCAGACATGCCGGGAATTTGCTGAAAAGGAGCTGGTTCCCATTGCAGCCCAGGTAGACAAGGAGCATC", &alphabet));
      sites.addSequence(
        bpp::BasicSequence("D", "TTCCAGACATGCCGGGACTTTACCGAGAAGGAGTTGTTTTCCATTGCAGCCCAGGTGGATAAGGAACATC", &alphabet));

      // Set of parameters to apply to tree + model
      paramModel1.addParameter(bpp::Parameter("T92.kappa", 0.1));
      paramModel2.addParameter(bpp::Parameter("T92.kappa", 0.2));
      paramBrLen1.addParameter(bpp::Parameter("BrLen1", 0.1));
      paramBrLen2.addParameter(bpp::Parameter("BrLen1", 0.2));

#ifndef NO_WARMUP
      // Warm up with eigen dummy computation
      double d = 0;
      for (int i = 0; i < 500; ++i)
        d += (Eigen::MatrixXd::Random(200, 200) * Eigen::MatrixXd::Random(200, 200)).determinant();
      static_cast<void>(d);
#endif
    }
  };
}

#ifdef ENABLE_OLD
TEST_CASE("old")
{
  const CommonStuff c;
  auto ts = timingStart();
  auto model = new bpp::T92(&c.alphabet, 3.);
  auto distribution = new bpp::ConstantDistribution(1.0);
  auto tree = std::unique_ptr<bpp::TreeTemplate<bpp::Node>>(bpp::TreeTemplateTools::parenthesisToTree(c.treeStr));
  bpp::RHomogeneousTreeLikelihood llh(*tree, c.sites, model, distribution, false, false);
  timingEnd(ts, "old_setup");
  ts = timingStart();
  llh.initialize();
  auto logLik = llh.getValue();
  timingEnd(ts, "old_init_value");
  printLik(logLik, "old_init_value");

  do_param_changes_multiple_times_legacy(llh, "old_param_model_change", c.paramModel1, c.paramModel2);
  do_param_changes_multiple_times_legacy(llh, "old_param_brlen_change", c.paramBrLen1, c.paramBrLen2);
}
#endif

#ifdef ENABLE_NEW
TEST_CASE("new")
{
  const CommonStuff c;
  auto ts = timingStart();
  auto model = new bpp::T92(&c.alphabet, 3.);
  bpp::Newick reader;
  auto phyloTree = std::unique_ptr<bpp::PhyloTree>(reader.parenthesisToPhyloTree(c.treeStr, false, "", false, false));
  auto paramPhyloTree = new bpp::ParametrizablePhyloTree(*phyloTree);
  auto process =
    std::unique_ptr<bpp::SimpleSubstitutionProcess>(new bpp::SimpleSubstitutionProcess(model, paramPhyloTree));
  auto likelihoodCompStruct = std::unique_ptr<bpp::RecursiveLikelihoodTreeCalculation>(
    new bpp::RecursiveLikelihoodTreeCalculation(c.sites, process.get(), false, true));
  bpp::SingleProcessPhyloLikelihood llh(process.get(), likelihoodCompStruct.release());
  timingEnd(ts, "new_setup");
  ts = timingStart();
  llh.computeLikelihood();
  auto logLik = llh.getValue();
  timingEnd(ts, "new_init_value");
  printLik(logLik, "new_init_value");

  do_param_changes_multiple_times_legacy(llh, "new_param_model_change", c.paramModel1, c.paramModel2);
  do_param_changes_multiple_times_legacy(llh, "new_param_brlen_change", c.paramBrLen1, c.paramBrLen2);
}
#endif

#ifdef ENABLE_DF
TEST_CASE("df")
{
  const CommonStuff c;
  auto ts = timingStart();
  // Read tree structure
  auto tree = std::unique_ptr<bpp::TreeTemplate<bpp::Node>>(bpp::TreeTemplateTools::parenthesisToTree(c.treeStr));
  auto treeData = bpp::Phyl::convertTreeTemplate(*tree);

  // Model
  auto model = bpp::Phyl::ModelNode::create(std::unique_ptr<bpp::SubstitutionModel>(new bpp::T92(&c.alphabet, 3.)));

  // Create a specification
  auto branchLengthMap =
    bpp::make_frozen(bpp::Topology::make_branch_parameter_map_from_value_map(*treeData.branchLengths));
  auto modelMap = bpp::make_frozen(bpp::Topology::make_uniform_branch_value_map(
    treeData.topology, bpp::DF::convertRef<bpp::DF::Value<const bpp::SubstitutionModel*>>(model)));
  auto process = bpp::Phyl::Process{treeData.topology, branchLengthMap, modelMap, c.alphabet.getSize()};
  auto sequenceMap = bpp::Phyl::makeSequenceMap(*treeData.nodeNames, c.sites);
  auto likParams = bpp::Phyl::LikelihoodParameters{process, sequenceMap};

  auto logLikNode =
    bpp::DF::convertRef<bpp::DF::Value<double>>(bpp::DF::instantiateNodeSpec(bpp::Phyl::LogLikelihoodSpec{likParams}));
  timingEnd(ts, "df_setup");
  ts = timingStart();
  auto logLik = logLikNode->getValue();
  timingEnd(ts, "df_init_value");
  printLik(logLik, "df_init_value");

  std::ofstream fd("df_debug");
  bpp::DF::debugDag(fd, logLikNode);

  // Change parameters
  do_param_changes_multiple_times_df(logLikNode, "df_param_model_change", model->getParameter("kappa"), 0.1, 0.2);

  auto topologyIdOfPhyloNode1 = treeData.treeTemplateNodeIndexes->index(1).value();
  auto& brlen1Param = branchLengthMap->access(treeData.topology->node(topologyIdOfPhyloNode1).fatherBranch()).value();
  do_param_changes_multiple_times_df(logLikNode, "df_param_brlen_change", brlen1Param, 0.1, 0.2);
}
#endif
