//
// File: Phylogeny.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-12
// Last modified: 2018-03-10
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

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

#ifndef BPP_NEWPHYL_LIKELIHOODEXAMPLE_H
#define BPP_NEWPHYL_LIKELIHOODEXAMPLE_H

#if 0
namespace bpp {
  /* Phylogeny DF graph construction functions.
   * TODO doc
   */
  DF::ValueRef<MatrixDouble> makeConditionalLikelihoodNode (const TreeTopologyInterface & tree,
                                                            const SequenceNodeAccess & sequenceNodes,
                                                            const BranchLengthNodeAccess & brlenNodes,
                                                            const ModelNodeAccess & modelNodes,
                                                            TopologyNodeIndex node);
  DF::ValueRef<MatrixDouble> makeForwardLikelihoodNode (const TreeTopologyInterface & tree,
                                                        const SequenceNodeAccess & sequenceNodes,
                                                        const BranchLengthNodeAccess & brlenNodes,
                                                        const ModelNodeAccess & modelNodes,
                                                        TopologyBranchIndex branch);

} // namespace bpp
#endif

#include <Bpp/Exceptions.h>

#include "Bpp/NewPhyl/Likelihood.h"
#include "Bpp/Phyl/Tree/PhyloTree.h"
#include "Bpp/Phyl/Model/SubstitutionModel.h"

#include <map>

namespace bpp {
  // Store interesting nodes of the likelihood example
  struct SimpleLikelihoodNodes {
    dataflow::ValueRef<double> totalLogLikelihood;

    std::map<PhyloTree::EdgeIndex, std::shared_ptr<dataflow::NumericMutable<double>>> branchLengthValues;
  };

  /* Build a likelihood computation dataflow graph for a simple example.
   *
   * The same model is used everywhere for simplicity.
   * In a real case, something like a map<EdgeIndex, shared_ptr<Model>> would give the model for each branch.
   *
   * In this example, a new leaf NumericMutable is generated for each branch length.
   * The set of parameters (branch lengths) is returned in the branchLengthValues map.
   * In a real case, something like a map<EdgeIndex, ValueRef<double>> would provide branch lengths.
   * The branch length values can be provided by any computation, or as a leaf NumericMutable node.
   */
  inline SimpleLikelihoodNodes makeSimpleLikelihoodNodes (const PhyloTree & tree,
                                                          std::shared_ptr<dataflow::ConfiguredModel> model,
                                                          std::size_t nbSite) {
    SimpleLikelihoodNodes r;

    dataflow::Context c;

    const auto nbState = model->getValue()->getNumberOfStates();
    const auto likelihoodMatrixDim = conditionalLikelihoodDimension (nbState, nbSite);

    // Recursive intermediate functions TODO
    auto makeConditionalLikelihoodNode = [&](PhyloTree::NodeIndex index) {};
    auto makeForwardLikelihoodNode = [&](PhyloTree::EdgeIndex index) {};

    // Build conditional likelihoods up to root recursively.
    if (!tree.isRooted ()) {
      throw Exception ("PhyloTree must be rooted");
    }
    auto rootConditionalLikelihoods = makeConditionalLikelihoodNode (tree.getRootIndex ());

    // Combine them to equilibrium frequencies to get the log likelihood
    auto equFreqs = dataflow::EquilibriumFrequenciesFromModel::create (
      c, {model}, equilibriumFrequenciesDimension (nbState));
    auto siteLikelihoods = dataflow::LikelihoodFromRootConditional::create (
      c, {equFreqs, rootConditionalLikelihood}, rowVectorDimension (nbSite));
    auto totalLogLikelihood =
      dataflow::TotalLogLikelihood::create (c, {siteLikelihoods}, rowVectorDimension (nbSite));

    // We want -log(likelihood)
    r.totalLogLikelihood = dataflow::CWiseNegate<double> (c, {totalLogLikelihood}, Dimension<double> ());

    return r;
  }

} // namespace bpp

#endif // BPP_NEWPHYL_LIKELIHOODEXAMPLE_H
