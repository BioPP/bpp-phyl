//
// File: HomogeneousLikelihoodExample.h
// Authors:
//   Laurent Guéguen (2018)
// Created: mercredi 24 octobre 2018, à 18h 37
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

#ifndef BPP_NEWPHYL_HOMOGENEOUS_LIKELIHOOD_EXAMPLE_H
#define BPP_NEWPHYL_HOMOGENEOUS_LIKELIHOOD_EXAMPLE_H

#include "Bpp/NewPhyl/DataFlowWrappers.h"
#include "Bpp/NewPhyl/Model.h"
#include "Bpp/NewPhyl/DiscreteDistribution.h"
#include "Bpp/NewPhyl/FrequenciesSet.h"
#include <Bpp/Seq/Container/AlignedValuesContainer.h>

#include "Bpp/Phyl/Model/SubstitutionModel.h"
#include "Bpp/Phyl/Model/FrequenciesSet/FrequenciesSet.h"
#include "Bpp/Phyl/Tree/PhyloTree.h"

#include "Bpp/NewPhyl/PhyloTree_BrRef.h"
#include "Bpp/NewPhyl/ForwardLikelihoodTree.h"

#include <unordered_map>
#include <list>

/* This file contains temporary helpers and wrappers.
 * They are used to bridge the gap between bpp::dataflow stuff and the rest of bpp.
 * They have only been used (and thus tested) for a single likelihood example.
 * They do not deal with all of bpp features, which is why they are only temporary.
 *
 * Ultimately, stuff in this file should be changed to a new system to describe phylogenic computations, which
 * would generate dataflow graphs to do the actual computations.
 */

namespace bpp {

  // Recursion helper class.
  // This stores state used by the two mutually recursive functions used to generate cond lik nodes.
  // The struct is similar to how a lambda is done internally, and allow the function definitions to be short.
  // The pure function equivalent has seven arguments, which is horrible.
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
  
  inline dataflow::ValueRef<double> makeHomogeneousLikelihoodNodes (dataflow::Context & c,
                                                                    const AlignedValuesContainer & sites,
                                                                    std::shared_ptr<dataflow::PhyloTree_BrRef> tree,
                                                                    std::shared_ptr<dataflow::ConfiguredModel> model,
                                                                    std::shared_ptr<dataflow::ConfiguredFrequenciesSet> rootFreqs = 0,
                                                                    std::shared_ptr<dataflow::ConfiguredDistribution> rates = 0) {

    std::size_t nbSite = sites.getNumberOfSites();
    
    const auto nbState = model->getTargetValue()->getNumberOfStates (); // Number of stored state values !
    // Build conditional likelihoods up to root recursively.
    if (!tree->isRooted ()) {
      throw Exception ("PhyloTree must be rooted");
    }

    // Set root frequencies 
    auto rFreqs = rootFreqs?dataflow::ConfiguredParametrizable::createVector<dataflow::ConfiguredFrequenciesSet, dataflow::FrequenciesFromFrequenciesSet> (
      c, {rootFreqs}, rowVectorDimension (Eigen::Index (nbState))):
      dataflow::ConfiguredParametrizable::createVector<dataflow::ConfiguredModel, dataflow::EquilibriumFrequenciesFromModel> (
        c, {model}, rowVectorDimension (Eigen::Index (nbState)));

    dataflow::ValueRef<Eigen::RowVectorXd> siteLikelihoods;
    
    if (rates)
    {
      auto brlenmap = bpp::dataflow::createBrLenMap(c, *tree);
      
      uint nbCat=(uint)rates->getTargetValue()->getNumberOfCategories();

      std::vector<std::shared_ptr<bpp::dataflow::Node>> vLogRoot;
      
      for (uint nCat=0; nCat<nbCat; nCat++)
      {
        auto catRef = bpp::dataflow::CategoryFromDiscreteDistribution::create(c, {rates}, nCat);
  
        auto rateBrLen = bpp::dataflow::multiplyBrLenMap(c, *tree, std::move(catRef));

        auto treeCat = std::shared_ptr<bpp::dataflow::PhyloTree_BrRef>(new bpp::dataflow::PhyloTree_BrRef(*tree, std::move(rateBrLen)));

        bpp::dataflow::ForwardLikelihoodTree flt(c, treeCat, model->getTargetValue()->getStateMap());
        flt.initialize(sites);
        auto rootConditionalLikelihoods=flt.getForwardLikelihoodArray(treeCat->getRootIndex ());

        auto siteLikelihoodsCat = dataflow::LikelihoodFromRootConditional::create (
          c, {rFreqs, rootConditionalLikelihoods}, rowVectorDimension (Eigen::Index (nbSite)));

        vLogRoot.push_back(std::move(siteLikelihoodsCat));
      }

      auto catProb = bpp::dataflow::ProbabilitiesFromDiscreteDistribution::create(c, {rates});

      vLogRoot.push_back(catProb);
      
      siteLikelihoods = bpp::dataflow::CWiseMean<Eigen::RowVectorXd, bpp::dataflow::ReductionOf<Eigen::RowVectorXd>, Eigen::RowVectorXd>::create(c, std::move(vLogRoot), rowVectorDimension (Eigen::Index(nbSite)));
    }
    else
    {
      bpp::dataflow::ForwardLikelihoodTree flt(c, tree, model->getTargetValue()->getStateMap());
      flt.initialize(sites);

      auto rootConditionalLikelihoods=flt.getForwardLikelihoodArray(tree->getRootIndex ());

      siteLikelihoods = dataflow::LikelihoodFromRootConditional::create (
        c, {rFreqs, rootConditionalLikelihoods}, rowVectorDimension (Eigen::Index (nbSite)));
    }
    
    auto totalLogLikelihood =
      dataflow::SumOfLogarithms<Eigen::RowVectorXd>::create (c, {siteLikelihoods}, rowVectorDimension (Eigen::Index (nbSite)));

// We want -log(likelihood)
    auto totalNegLogLikelihood =
      dataflow::CWiseNegate<double>::create (c, {totalLogLikelihood}, Dimension<double> ());

    return totalNegLogLikelihood;
  }

} // namespace bpp

#endif // BPP_NEWPHYL_HOMOGENEOUS_LIKELIHOOD_EXAMPLE_H

