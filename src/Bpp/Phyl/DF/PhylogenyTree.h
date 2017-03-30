// File: PhylogenyTree.h
// Authors:
//   Francois Gindraud (2017)
// Created: 15/03/2017

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

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

#ifndef _DF_PHYLOGENY_TREE_H_
#define _DF_PHYLOGENY_TREE_H_

//#include <Bpp/Phyl/NewLikelihood/SubstitutionProcess.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Seq/Sequence.h>

#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Phyl/DF/DataFlowComputationClasses.h>
#include <Bpp/Utils/Cpp14.h>
#include <Bpp/Utils/ForRange.h>

#include <cstdint>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

namespace bpp
{
  namespace New
  {
    // Forward declare manipulator
    template <typename NodeType, typename BranchType>
    class TreeManipulator;

    // Base class of extended trees.
    // Contains a tree topology only.
    class ExtendableTree
    {
    public:
      using IndexType = std::uint32_t;
      static constexpr IndexType invalidIndex = std::numeric_limits<IndexType>::max();

      // Base node type
      class Node
      {
      private:
        IndexType id_;
        IndexType fatherBranch_{invalidIndex};
        std::vector<IndexType> childBranches_;

      public:
        void setIndex(IndexType id) { id_ = id; }
        IndexType getIndex() const { return id_; }
        virtual ~Node() = default;
      };
      // Base branch type
      class Branch
      {
      private:
        IndexType id_;
        IndexType fatherNode_;
        IndexType childNode_;

      public:
        virtual ~Branch() = default;
      };
      // Store tree.
      // Ptr to node and branches to allow virtual dispatching.
      // Also because tree elements will ultimately store DF nodes that are non movable.
      struct Storage
      {
        std::vector<std::unique_ptr<Node>> nodes_;
        std::vector<std::unique_ptr<Branch>> branches_;
      };

      // Method that allow to access the tree.
      // Should be overriden in each extended tree subclass to return a manipulator with extended Node and Branch.
      TreeManipulator<Node, Branch> tree(void);

    protected:
      Storage tree_;
    };

    // Class that manipulates the tree (access, modification, etc).
    // Access the extensions of basic Node and Branch given in the templates.
    template <typename NodeType, typename BranchType>
    class TreeManipulator
    {
    public:
      using IndexType = ExtendableTree::IndexType;
      static_assert(std::is_base_of<ExtendableTree::Node, NodeType>::value,
                    "NodeType must inherit from ExtendableTree::Node");
      static_assert(std::is_base_of<ExtendableTree::Branch, BranchType>::value,
                    "BranchType must inherit from ExtendableTree::Branch");

      TreeManipulator(ExtendableTree::Storage& storage)
        : tree_(storage)
      {
      }

      // Downcast access (cannot be static due to language rules)
      NodeType& node(IndexType index) const { return dynamic_cast<NodeType&>(*tree_.nodes_[index]); }
      BranchType& branch(IndexType index) const { return dynamic_cast<BranchType&>(*tree_.branches_[index]); }

      std::size_t nbNodes() const { return tree_.nodes_.size(); }
      std::size_t nbBranches() const { return tree_.branches_.size(); }

      NodeType& addNode(void) const
      {
        auto id = IndexType(tree_.nodes_.size());
        tree_.nodes_.emplace_back(Cpp14::make_unique<NodeType>());
        auto& n = node(id);
        n.setIndex(id);
        return n;
      }

    private:
      ExtendableTree::Storage& tree_;
    };

    // Need to put it after TreeManipulator definition
    auto ExtendableTree::tree(void) -> TreeManipulator<Node, Branch> { return {tree_}; }

    // Adds branch lengths
    class PhylogenyTree : public virtual ExtendableTree
    {
    public:
      using ExtendableTree::Node;
      class Branch : public virtual ExtendableTree::Branch
      {
      public:
        void setLength(double length) { length_.setValue(length); }
        double getLength(void) { return length_.getValue(); }
      protected:
        DF::ParameterNode<double> length_;
      };

      TreeManipulator<Node, Branch> tree(void) { return {tree_}; }
    };

    // Add a model per branch
    class PhylogenyProcess : public virtual PhylogenyTree
    {
    public:
      using PhylogenyTree::Node;
      class Branch : public virtual PhylogenyTree::Branch
      {
      public:
        Branch()
          : PhylogenyTree::Branch()
        {
          // Connect to brlen and model (value, diff, diff2)
          modelMatrixForBranchLength_.setDependency<ModelMatrixComputation::BranchLen>(length_);
          modelMatrixForBranchLength_.setDependency<ModelMatrixComputation::Model>(model_);
          modelDMatrixForBranchLength_.setDependency<ModelDMatrixComputation::BranchLen>(length_);
          modelDMatrixForBranchLength_.setDependency<ModelDMatrixComputation::Model>(model_);
          modelDDMatrixForBranchLength_.setDependency<ModelDDMatrixComputation::BranchLen>(length_);
          modelDDMatrixForBranchLength_.setDependency<ModelDDMatrixComputation::Model>(model_);
        }

        void setModel(SubstitutionModel* model) { model_.setValue(model); }
        SubstitutionModel* getModel(void) { return model_.getValue(); }

      protected:
        DF::ParameterNode<SubstitutionModel*> model_;

        // Compute matrix from model
        struct ModelMatrixComputation
        {
          using ResultType = DefaultMatrix<double>;
          enum
          {
            BranchLen,
            Model
          };
          using ArgumentTypes = std::tuple<double, SubstitutionModel*>;
          static void compute(ResultType& result, double brLen, SubstitutionModel* model)
          {
            result = model->getPij_t(brLen);
          }
        };
        DF::HeterogeneousComputationNode<ModelMatrixComputation> modelMatrixForBranchLength_;

        // Compute diff matrix from model
        struct ModelDMatrixComputation
        {
          using ResultType = DefaultMatrix<double>;
          enum
          {
            BranchLen,
            Model
          };
          using ArgumentTypes = std::tuple<double, SubstitutionModel*>;
          static void compute(ResultType& result, double brLen, SubstitutionModel* model)
          {
            result = model->getdPij_dt(brLen);
          }
        };
        DF::HeterogeneousComputationNode<ModelDMatrixComputation> modelDMatrixForBranchLength_;

        // Compute diff2 matrix from model
        struct ModelDDMatrixComputation
        {
          using ResultType = DefaultMatrix<double>;
          enum
          {
            BranchLen,
            Model
          };
          using ArgumentTypes = std::tuple<double, SubstitutionModel*>;
          static void compute(ResultType& result, double brLen, SubstitutionModel* model)
          {
            result = model->getd2Pij_dt2(brLen);
          }
        };
        DF::HeterogeneousComputationNode<ModelDDMatrixComputation> modelDDMatrixForBranchLength_;
      };

      TreeManipulator<Node, Branch> tree(void) { return {tree_}; }
    };

    using LikelihoodValues = std::vector<double>;
    using LikelihoodValuesBySite = std::vector<LikelihoodValues>;

    class PhylogenySiteObservations : public virtual ExtendableTree
    {
    public:
      class Node : public virtual ExtendableTree::Node
      {
      public:
        void setSequence(const Sequence* seq) { sequence_ = seq; }
        bool hasSequence() const { return sequence_ != nullptr; }
        const std::string& sequenceName() const { return sequence_->getName(); }

      private:
        const Sequence* sequence_{};
      };
      using Branch = ExtendableTree::Branch;

      TreeManipulator<Node, Branch> tree(void) { return {tree_}; }
    };

    // Add conditional likelihoods for all sites
    class PhylogenyConditionalLikelihood : public virtual PhylogenyProcess, public virtual PhylogenySiteObservations
    {
    public:
      class Node : public virtual PhylogenyProcess::Node, public virtual PhylogenySiteObservations::Node
      {
      protected:
        // Conditionnal lik computed from childrens
        struct ConditionalLikelihoodComputation
        {
          // FIXME matrix ?
          using ResultType = LikelihoodValuesBySite;
          using ArgumentType = LikelihoodValuesBySite;
          static void reset(ResultType& result)
          {
            for (auto& site : result)
              for (auto& likOfChar : site)
                likOfChar = 1.;
          }
          static void reduce(ResultType& result, const ArgumentType& forwardLik)
          {
            for (auto siteIndex : makeRange(forwardLik.size()))
            {
              auto& siteCondLik = result[siteIndex];
              auto& siteFwdLik = forwardLik[siteIndex];
              for (auto c : makeRange(siteCondLik.size()))
                siteCondLik[c] *= siteFwdLik[c];
            }
          }
        };
        DF::ReductionComputationNode<ConditionalLikelihoodComputation> conditionalLikelihood_;

        //
      };
      class Branch : public virtual PhylogenyProcess::Branch, public virtual PhylogenySiteObservations::Branch
      {
      public:
      };

      TreeManipulator<Node, Branch> tree(void) { return {tree_}; }
    };
  }
} // end of namespace bpp.

#endif //_DF_PHYLOGENY_TREE_H_
