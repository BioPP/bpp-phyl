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

#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/Model.h>
#include <Bpp/NewPhyl/Phylogeny.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <fstream>

TEST_CASE("test")
{
  // Read tree structure
  bpp::Newick reader;
  auto phyloTree = std::unique_ptr<bpp::PhyloTree>(
    reader.parenthesisToPhyloTree("((A:0.01, B:0.02):0.03,C:0.01,D:0.1);", false, "", false, false));
  auto phyloTreeData = bpp::Phyl::convertPhyloTree(*phyloTree);

  std::ofstream ft("topology_debug");
  bpp::Topology::debugTree(ft, phyloTreeData.topology);

  // Init sequences
  const auto& alphabet = bpp::AlphabetTools::DNA_ALPHABET;
  bpp::VectorSiteContainer sites(&alphabet);
  sites.addSequence(
    bpp::BasicSequence("A", "ATCCAGACATGCCGGGACTTTGCAGAGAAGGAGTTGTTTCCCATTGCAGCCCAGGTGGATAAGGAACAGC", &alphabet));
  sites.addSequence(
    bpp::BasicSequence("B", "CGTCAGACATGCCGTGACTTTGCCGAGAAGGAGTTGGTCCCCATTGCGGCCCAGCTGGACAGGGAGCATC", &alphabet));
  sites.addSequence(
    bpp::BasicSequence("C", "GGTCAGACATGCCGGGAATTTGCTGAAAAGGAGCTGGTTCCCATTGCAGCCCAGGTAGACAAGGAGCATC", &alphabet));
  sites.addSequence(
    bpp::BasicSequence("D", "TTCCAGACATGCCGGGACTTTACCGAGAAGGAGTTGTTTTCCATTGCAGCCCAGGTGGATAAGGAACATC", &alphabet));

  // Model
  auto model = bpp::DF::Value<const bpp::SubstitutionModel*>::create<bpp::Phyl::ModelNode>(
    std::unique_ptr<bpp::SubstitutionModel>(new bpp::T92(&alphabet, 3.)));

  // Init
  auto branchLengthMap =
    bpp::make_frozen(bpp::Topology::make_branch_parameter_map_from_value_map(*phyloTreeData.branchLengths));
  auto modelMap = bpp::make_frozen(bpp::Topology::make_uniform_branch_value_map(phyloTreeData.topology, model));
  auto process = bpp::Phyl::Process{phyloTreeData.topology, branchLengthMap, modelMap, alphabet.getSize()};

  // Make leaf data
  auto leafDataMapTmp =
    bpp::make_freezable<bpp::Topology::NodeValueMap<bpp::DF::Parameter<const bpp::Sequence*>>>(phyloTreeData.topology);
  for (auto i : bpp::index_range(*leafDataMapTmp))
    leafDataMapTmp->access(i) = phyloTreeData.nodeNames->access(i).map([&sites](const std::string& name) {
      return bpp::DF::Parameter<const bpp::Sequence*>::create(&sites.getSequence(name));
    });
  auto leafDataMap = std::move(leafDataMapTmp).freeze();

  // Finally, likelihood parameters
  auto likParams = bpp::Phyl::LikelihoodParameters{process, leafDataMap, sites.getNumberOfSites()};

  bpp::DF::Value<double> logLikNode{bpp::DF::instantiateNodeSpec(bpp::Phyl::LogLikelihoodSpec{likParams})};

  std::ofstream fd("df_debug");
  bpp::DF::debugDag(fd, logLikNode);

  std::cout << "Log lik = " << logLikNode.getValue() << "\n";
}
