//
// File: LikelihoodExample.h
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

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

#include "Bpp/NewPhyl/Likelihood.h"
#include "Bpp/Phyl/Model/SubstitutionModel.h"
#include "Bpp/Phyl/Tree/PhyloTree.h"

#include <unordered_map>

/* This file contains temporary helpers and wrappers.
 * They are used to bridge the gap between bpp::dataflow stuff and the rest of bpp.
 * They have only been used (and thus tested) for a single likelihood example.
 * They do not deal with all of bpp features, which is why they are only temporary.
 *
 * Ultimately, stuff in this file should be changed to a new system to describe phylogenic computations, which
 * would generate dataflow graphs to do the actual computations.
 */
namespace bpp {
  // Store interesting nodes of the likelihood example
  struct SimpleLikelihoodNodes {
    dataflow::ValueRef<double> totalLogLikelihood;

    std::unordered_map<PhyloTree::EdgeIndex, std::shared_ptr<dataflow::NumericMutable<double>>>
      branchLengthValues;
  };

  // Recursion helper class.
  // This stores state used by the two mutually recursive functions used to generate cond lik nodes.
  // The struct is similar to how a lambda is done internally, and allow the function definitions to be short.
  // The pure function equivalent has seven arguments, which is horrible.
  struct SimpleLikelihoodNodesHelper {
    dataflow::Context & c;
    SimpleLikelihoodNodes & r;
    std::shared_ptr<dataflow::ConfiguredModel> model;
    const PhyloTree & tree;
    const VectorSiteContainer & sites;
    MatrixDimension likelihoodMatrixDim;
    std::size_t nbState;
    std::size_t nbSite;

    dataflow::NodeRef makeInitialConditionalLikelihood (const std::string & sequenceName) {
      /* FIXME Generate the matrix of {0,1} for each (state, site).
       * I am not sure what is the current interface class to use in the bpp_seq stuff.
       * Thus I selected VectorSiteContainer which is an implementation class.
       * This is certainly not the most permissive interface.
       * In addition, I iterate on raw states, but newlik code seems to use StateMap conversion.
       * I do not handle any of that in this example !
       * It should also be checked that edge models have the same state space, etc...
       */
      const auto sequenceIndex = sites.getSequencePosition (sequenceName);
      Eigen::MatrixXd initCondLik (nbState, nbSite);
      for (std::size_t site = 0; site < nbSite; ++site) {
        for (std::size_t state = 0; state < nbState; ++state) {
          initCondLik (Eigen::Index (state), Eigen::Index (site)) =
            sites.getStateValueAt (site, sequenceIndex, int(state));
        }
      }
      return dataflow::NumericConstant<Eigen::MatrixXd>::create (c, std::move (initCondLik));
    }

    dataflow::NodeRef makeForwardLikelihoodNode (PhyloTree::EdgeIndex index) {
      const auto initBrlen = tree.getEdge (index)->getLength ();
      auto brlen = dataflow::NumericMutable<double>::create (c, initBrlen);
      r.branchLengthValues.emplace (index, brlen);

      auto childConditionalLikelihood = makeConditionalLikelihoodNode (tree.getSon (index));
      auto transitionMatrix =
        dataflow::TransitionMatrixFromModel::create (c, {model, brlen}, transitionMatrixDimension (nbState));
      return dataflow::ForwardLikelihoodFromConditional::create (
        c, {transitionMatrix, childConditionalLikelihood}, likelihoodMatrixDim);
    }

    dataflow::NodeRef makeConditionalLikelihoodNode (PhyloTree::NodeIndex index) {
      const auto childBranchIndexes = tree.getBranches (index);
      if (childBranchIndexes.empty ()) {
        return makeInitialConditionalLikelihood (tree.getNode (index)->getName ());
      } else {
        dataflow::NodeRefVec deps (childBranchIndexes.size ());
        for (std::size_t i = 0; i < childBranchIndexes.size (); ++i) {
          deps[i] = makeForwardLikelihoodNode (childBranchIndexes[i]);
        }
        return dataflow::ConditionalLikelihoodFromChildrenForward::create (c, std::move (deps),
                                                                           likelihoodMatrixDim);
      }
    }
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
  
  inline SimpleLikelihoodNodes makeSimpleLikelihoodNodes (dataflow::Context & c, const PhyloTree & tree,
                                                          const VectorSiteContainer & sites,
                                                          std::shared_ptr<dataflow::ConfiguredModel> model) {
    const auto nbState = model->getValue()->getNumberOfStates (); // Number of stored state values !
    const auto nbSite = sites.getNumberOfSites ();
    const auto likelihoodMatrixDim = conditionalLikelihoodDimension (nbState, nbSite);
    SimpleLikelihoodNodes r;

    // Build conditional likelihoods up to root recursively.
    if (!tree.isRooted ()) {
      throw Exception ("PhyloTree must be rooted");
    }

    // Recursively generate dataflow graph for conditional likelihood using helper struct.
    SimpleLikelihoodNodesHelper helper{c, r, model, tree, sites, likelihoodMatrixDim, nbState, nbSite};
    auto rootConditionalLikelihoods = helper.makeConditionalLikelihoodNode (tree.getRootIndex ());

    // Combine them to equilibrium frequencies to get the log likelihood
    auto equFreqs = dataflow::EquilibriumFrequenciesFromModel::create (
      c, {model}, rowVectorDimension (Eigen::Index (nbState)));
    auto siteLikelihoods = dataflow::LikelihoodFromRootConditional::create (
      c, {equFreqs, rootConditionalLikelihoods}, rowVectorDimension (Eigen::Index (nbSite)));
    auto totalLogLikelihood =
      dataflow::TotalLogLikelihood::create (c, {siteLikelihoods}, rowVectorDimension (Eigen::Index (nbSite)));

    // We want -log(likelihood)
    r.totalLogLikelihood =
      dataflow::CWiseNegate<double>::create (c, {totalLogLikelihood}, Dimension<double> ());
    return r;
  }

  /* Wraps a dataflow::NumericMutable<double> as a bpp::Parameter.
   * 2 values exist: the one in the node, and the one in bpp::Parameter.
   * The dataflow one is considered to be the reference.
   *
   * FIXME This is a temporary system:
   * - It ignores all the Parameter listener system.
   * - Synchronization between the 2 values is a best effort.
   * - bpp::Parameter should be improved with respect to semantics, because it currently is a mess of multiple
   * systems with sharing, some shared_ptr stuff, etc...
   */
  class DataFlowParameter : public Parameter {
  public:
    DataFlowParameter (const std::string & name,
                       std::shared_ptr<dataflow::NumericMutable<double>> existingNode)
      : Parameter (name, existingNode->getValue ()), node_ (std::move (existingNode)) {}

    // Parameter boilerplate
    DataFlowParameter * clone () const override { return new DataFlowParameter (*this); }

    // Override value access
    double getValue () const override { return node_->getValue (); }
    void setValue (double v) override {
      Parameter::setValue (v);                  // Will apply possible constraints
      node_->setValue (Parameter::getValue ()); // Get constrained value
    }

    dataflow::NumericMutable<double> & node () const { return *node_; }

  private:
    std::shared_ptr<dataflow::NumericMutable<double>> node_;
  };

  /* Wraps a dataflow graph as a function: resultNode = f(variableNodes).
   *
   * FIXME This temporary interface to bpp::DerivableSecondOrder stuff should be improved.
   *
   * Any bpp::Parameter can be given in the bpp::ParameterList, but only DataFlowParameter are supported
   * because we need to have dataflow nodes.
   * No specific check is done, in case of error it will be a std::bad_cast from dynamic_cast.
   *
   * In addition, as we need a context for derivation but the bpp legacy API does not support it, a reference
   * is stored in the class which is dangerous with respect to lifetime.
   */
  class DataFlowFunction : public DerivableSecondOrder {
  private:
    dataflow::Context & context_;

    // Store nodes
    dataflow::ValueRef<double> resultNode_;
    ParameterList variableNodes_;

    // Cache generated nodes representing derivatives, to avoid recreating them every time.
    // Using the mutable keyword because the table must be changed even in const methods.
    struct StringPairHash {
      std::size_t operator() (const std::pair<std::string, std::string> & p) const {
        std::hash<std::string> strHash{};
        return strHash (p.first) ^ (strHash (p.second) << 1);
      }
    };
    mutable std::unordered_map<std::string, dataflow::ValueRef<double>> firstOrderDerivativeNodes_;
    mutable std::unordered_map<std::pair<std::string, std::string>, dataflow::ValueRef<double>,
                               StringPairHash>
      secondOrderDerivativeNodes_;

  public:
    DataFlowFunction (dataflow::Context & context, dataflow::ValueRef<double> resultNode,
                      const ParameterList & variableNodes)
      : context_ (context), resultNode_ (std::move (resultNode)), variableNodes_ (variableNodes) {}

    // Legacy boilerplate
    DataFlowFunction * clone () const override { return new DataFlowFunction (*this); }

    // bpp::Parametrizable (prefix unused FIXME?)
    bool hasParameter (const std::string & name) const override { return variableNodes_.hasParameter (name); }
    const ParameterList & getParameters () const override { return variableNodes_; }
    const Parameter & getParameter (const std::string & name) const override {
      return variableNodes_.getParameter (name);
    }
    double getParameterValue (const std::string & name) const override {
      return variableNodes_.getParameterValue (name);
    }
    void setAllParametersValues (const ParameterList & params) override {
      return variableNodes_.setAllParametersValues (params);
    }
    void setParameterValue (const std::string & name, double value) override {
      variableNodes_.setParameterValue (name, value);
    }
    void setParametersValues (const ParameterList & params) override {
      variableNodes_.setParametersValues (params);
    }
    bool matchParametersValues (const ParameterList & params) override {
      return variableNodes_.matchParametersValues (params);
    }
    std::size_t getNumberOfParameters () const override { return variableNodes_.size (); }
    void setNamespace (const std::string &) override {}
    std::string getNamespace () const override { return {}; }
    std::string getParameterNameWithoutNamespace (const std::string & name) const override { return name; }

    // bpp::Function
    void setParameters (const ParameterList & params) override {
      variableNodes_.setParametersValues (params);
    }
    double getValue () const override { return resultNode_->getValue (); }

    // bpp::DerivableFirstOrder
    void enableFirstOrderDerivatives (bool) override {}
    bool enableFirstOrderDerivatives () const override { return true; }
    double getFirstOrderDerivative (const std::string & variable) const override {
      return firstOrderDerivativeNode (variable)->getValue ();
    }

    // bpp::DerivableSecondOrder
    void enableSecondOrderDerivatives (bool) override {}
    bool enableSecondOrderDerivatives () const override { return true; }
    double getSecondOrderDerivative (const std::string & variable) const override {
      return getSecondOrderDerivative (variable, variable);
    }
    double getSecondOrderDerivative (const std::string & variable1,
                                     const std::string & variable2) const override {
      return secondOrderDerivativeNode (variable1, variable2)->getValue ();
    }

    // Get nodes of derivatives directly
    dataflow::ValueRef<double> firstOrderDerivativeNode (const std::string & variable) const {
      const auto it = firstOrderDerivativeNodes_.find (variable);
      if (it != firstOrderDerivativeNodes_.end ()) {
        return it->second;
      } else {
        auto node = resultNode_->deriveAsValue (context_, accessVariableNode (variable));
        firstOrderDerivativeNodes_.emplace (variable, node);
        return node;
      }
    }
    dataflow::ValueRef<double> secondOrderDerivativeNode (const std::string & variable1,
                                                          const std::string & variable2) const {
      const auto key = std::make_pair (variable1, variable2);
      const auto it = secondOrderDerivativeNodes_.find (key);
      if (it != secondOrderDerivativeNodes_.end ()) {
        return it->second;
      } else {
        // Reuse firstOrderDerivative() to generate the first derivative with caching
        auto node =
          firstOrderDerivativeNode (variable1)->deriveAsValue (context_, accessVariableNode (variable2));
        secondOrderDerivativeNodes_.emplace (key, node);
        return node;
      }
    }

  private:
    static dataflow::NumericMutable<double> & accessVariableNode (const Parameter & param) {
      return dynamic_cast<const DataFlowParameter &> (param).node ();
    }
    dataflow::NumericMutable<double> & accessVariableNode (const std::string & name) const {
      return accessVariableNode (getParameter (name));
    }
  };
} // namespace bpp

#endif // BPP_NEWPHYL_LIKELIHOODEXAMPLE_H
