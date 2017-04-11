// File: test_extendable_tree.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 20/03/2017

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

#undef NDEBUG

// Shared stuff
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

// Old system
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
// Newlik
#include <Bpp/Phyl/Io/Newick.h>
// DF
#include <Bpp/Phyl/DF/PhylogenyTree.h>

#include <Bpp/Utils/Cpp14.h>
#include <Bpp/Utils/ForRange.h>
#include <chrono>
#include <iostream>
#include <map>

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
    constexpr std::size_t updatesNbIterations = 10;
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
}

struct ValuesToCompare
{
  double initialLikelihood{};
  ValuesToCompare() = default;
};

TEST_CASE("Compare likelihood computations with 3 methods")
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
  ConstantDistribution distribution(1.0);

  // Set of parameters to apply to tree + model
  ParameterList paramModel1;
  paramModel1.addParameter(Parameter("T92.kappa", 0.1));
  ParameterList paramModel2;
  paramModel2.addParameter(Parameter("T92.kappa", 0.2));

  ParameterList paramBrLen1;
  paramBrLen1.addParameter(Parameter("BrLen1", 0.1));
  ParameterList paramBrLen2;
  paramBrLen2.addParameter(Parameter("BrLen1", 0.2));

  // Old likelihood
  ValuesToCompare oldL;
  {
    auto ts = timingStart();
    auto tree = std::unique_ptr<TreeTemplate<Node>>(
      TreeTemplateTools::parenthesisToTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);"));
    RHomogeneousTreeLikelihood llh(*tree, sites, model.clone(), distribution.clone(), false, false);
    llh.initialize();
    oldL.initialLikelihood = llh.getValue();
    timingEnd(ts, "old_init_value");

    do_param_changes_multiple_times(llh, "old_param_model_change", paramModel1, paramModel2);
    do_param_changes_multiple_times(llh, "old_param_brlen_change", paramBrLen1, paramBrLen2);
  }
  // DF likelihood
  ValuesToCompare newL;
  {
    auto ts = timingStart();
    Newick reader;
    auto phyloTree = std::unique_ptr<PhyloTree>(
      reader.parenthesisToPhyloTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);", false, "", false, false));

    //
    auto& seq = sites.getSequence("A");
    auto nbSites = seq.size();
    auto nbStates = seq.getAlphabet()->getSize();
    auto* localModel = model.clone();

    // Create problem tree
    New::PhylogenyConditionalLikelihood lik;
    auto manipulator = lik.tree();

    // Fill tree
    std::map<PhyloTree::NodeIndex, New::PhylogenyProcess::IndexType> nodeConv;
    for (auto phylNodeId : phyloTree->getAllNodesIndexes())
    {
      auto& node = manipulator.addNode();
      std::cout << "node_create(id=" << node.getIndex() << ",phylid=" << phylNodeId << ")\n";
      nodeConv[phylNodeId] = node.getIndex();
      auto phylNode = phyloTree->getNode(phylNodeId);
      if (phylNode->hasName())
      {
        // Add sequence to corresponding node
        std::cout << "add_seq(id=" << node.getIndex() << ",seqName=" << phylNode->getName() << ")\n";
        node.setSequence(&sites.getSequence(phylNode->getName()));
      }
    }
    for (auto phylNodeId : phyloTree->getAllNodesIndexes())
    {
      // Add links (to father, doesn't matter as we iterate over all nodes)
      if (phyloTree->hasFather(phylNodeId))
      {
        auto phylNode = phyloTree->getNode(phylNodeId);
        auto phylFatherId = phyloTree->getNodeIndex(phyloTree->getFather(phylNode));
        auto phylBranch = phyloTree->getEdgeToFather(phylNodeId);
        auto& branch = manipulator.addBranch(nodeConv.at(phylFatherId), nodeConv.at(phylNodeId));
        branch.setLength(phylBranch->getLength());
        branch.setModel(localModel);
        std::cout << "branch_create(id=" << branch.getIndex() << ",from=" << branch.getFather()
                  << ",to=" << branch.getChild() << ",len=" << branch.getLength() << ")\n";
      }
    }
    std::cout << "root_id=" << nodeConv[phyloTree->getRootIndex()] << "\n";
    // Final init... (crappy but temporary)
    for (auto i : makeRange(manipulator.nbNodes()))
      manipulator.node(i).initStuff(nbSites, nbStates);
    for (auto i : makeRange(manipulator.nbBranches()))
      manipulator.branch(i).initStuff(nbSites, nbStates);

    newL.initialLikelihood = manipulator.node(nodeConv[phyloTree->getRootIndex()]).getLogLik();
    timingEnd(ts, "df_init");

    auto& v = manipulator.node(nodeConv[phyloTree->getRootIndex()]).conditionalLikelihood_.getValue();
    for (auto siteId : makeRange(v.size()))
    {
      std::cout << "Lik[site " << siteId << "] = { ";
      for (auto l : v[siteId])
        std::cout << l << " ";
      std::cout << "}\n";
    }
  }

  CHECK(doctest::Approx(oldL.initialLikelihood) == newL.initialLikelihood);
}
