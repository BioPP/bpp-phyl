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
        IndexType id_{invalidIndex};
        IndexType fatherBranch_{invalidIndex};
        std::vector<IndexType> childBranches_;

      public:
        void graphInit(IndexType id) { id_ = id; }
        IndexType getIndex() const { return id_; }
        void addChildBranch(IndexType branchId) { childBranches_.emplace_back(branchId); }
        void setFatherBranch(IndexType branchId) { fatherBranch_ = branchId; }
        virtual ~Node() = default;
      };
      // Base branch type
      class Branch
      {
      private:
        IndexType id_{invalidIndex};
        IndexType fatherNode_{invalidIndex};
        IndexType childNode_{invalidIndex};

      public:
        void graphInit(IndexType id, Node& father, Node& child)
        {
          id_ = id;
          fatherNode_ = father.getIndex();
          father.addChildBranch(id);
          childNode_ = child.getIndex();
          child.setFatherBranch(id);
        }
        IndexType getIndex() const { return id_; }
        IndexType getFather() const { return fatherNode_; }
        IndexType getChild() const { return childNode_; }
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

      // Downcast (cannot be static due to language rules)
      static NodeType& node_cast(ExtendableTree::Node& node) { return dynamic_cast<NodeType&>(node); }
      static const NodeType& node_cast(const ExtendableTree::Node& node) { return dynamic_cast<const NodeType&>(node); }
      static BranchType& branch_cast(ExtendableTree::Branch& branch) { return dynamic_cast<BranchType&>(branch); }
      static const BranchType& branch_cast(const ExtendableTree::Branch& branch)
      {
        return dynamic_cast<const BranchType&>(branch);
      }

      // Tree access
      ExtendableTree::Node& basic_node(IndexType index) const { return *tree_.nodes_[index]; }
      ExtendableTree::Branch& basic_branch(IndexType index) const { return *tree_.branches_[index]; }

      NodeType& node(IndexType index) const { return node_cast(basic_node(index)); }
      BranchType& branch(IndexType index) const { return branch_cast(basic_branch(index)); }

      IndexType nbNodes() const { return IndexType(tree_.nodes_.size()); }
      IndexType nbBranches() const { return IndexType(tree_.branches_.size()); }

      NodeType& addNode(void) const
      {
        auto id = IndexType(tree_.nodes_.size());
        tree_.nodes_.emplace_back(Cpp14::make_unique<NodeType>());
        auto& n = node(id);
        n.graphInit(id);
        return n;
      }

      BranchType& addBranch(NodeType& from, NodeType& to)
      {
        auto id = IndexType(tree_.branches_.size());
        tree_.branches_.emplace_back(Cpp14::make_unique<BranchType>());
        auto& br = branch(id);
        br.graphInit(id, from, to);
        return br;
      }
      BranchType& addBranch(IndexType from, IndexType to) { return addBranch(node(from), node(to)); }

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
        const Sequence& getSequence() const { return *sequence_; }
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
      public:
        void initStuff(std::size_t nbSites, std::size_t nbStates)
        {
          conditionalLikelihood_.initValue() = LikelihoodValuesBySite(nbSites, LikelihoodValues(nbStates));
          if (hasSequence())
          {
            // Setup the init likelihoods
            LikelihoodValuesBySite initLikBySite(nbSites);
            for (auto site : makeRange(initLikBySite.size()))
            {
              LikelihoodValues liks(nbStates, 0.);
              liks[getSequence().getValue(site)] = 1.;
              initLikBySite[site] = std::move(liks);
            }
            initLikelihoodFromData_.setValue(std::move(initLikBySite));
            // Link them to cond
            conditionalLikelihood_.addDependencyTo(initLikelihoodFromData_);
          }
        }

        double getLogLik(SubstitutionModel* model)
        {
          const auto& freqs = model->getFrequencies();
          auto& condLiksBySite = conditionalLikelihood_.getValue();
          double logLik = 0.;
          for (auto& likOfSite : condLiksBySite)
          {
            double lik = 0.;
            for (auto c : makeRange(likOfSite.size()))
              lik += likOfSite[c] * freqs[c];
            logLik += std::log(lik);
          }
          return -logLik;
        }

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

        DF::ParameterNode<LikelihoodValuesBySite> initLikelihoodFromData_;
      };

      class Branch : public virtual PhylogenyProcess::Branch, public virtual PhylogenySiteObservations::Branch
      {
      public:
        Branch()
          : PhylogenyProcess::Branch()
        {
          forwardLikelihood_.setDependency<ForwardLikelihoodComputation::BranchModel>(modelMatrixForBranchLength_);
        }

        void graphInit(IndexType id, Node& from, Node& to)
        {
          PhylogenyProcess::Branch::graphInit(id, from, to);
          from.conditionalLikelihood_.addDependencyTo(forwardLikelihood_);
          forwardLikelihood_.setDependency<ForwardLikelihoodComputation::ChildCondLik>(to.conditionalLikelihood_);
        }

        void initStuff(std::size_t nbSites, std::size_t nbStates)
        {
          forwardLikelihood_.initValue() = LikelihoodValuesBySite(nbSites, LikelihoodValues(nbStates));
        }

      protected:
        struct ForwardLikelihoodComputation
        {
          using ResultType = LikelihoodValuesBySite;
          enum
          {
            ChildCondLik,
            BranchModel
          };
          using ArgumentTypes = std::tuple<LikelihoodValuesBySite, DefaultMatrix<double>>;
          static void compute(ResultType& resBySite, const LikelihoodValuesBySite& childCondLikBySite,
                              const DefaultMatrix<double>& transitionMatrix)
          {
            for (auto site : makeRange(resBySite.size()))
            {
              auto& fwdLiks = resBySite[site];
              auto& childCondLik = childCondLikBySite[site];
              for (auto target_c : makeRange(fwdLiks.size()))
              {
                fwdLiks[target_c] = 0;
                for (auto c : makeRange(fwdLiks.size()))
                  fwdLiks[target_c] += transitionMatrix(target_c, c) * childCondLik[c];
              }
            }
          }
        };
        DF::HeterogeneousComputationNode<ForwardLikelihoodComputation> forwardLikelihood_;
      };

      TreeManipulator<Node, Branch> tree(void) { return {tree_}; }
    };
  }
} // end of namespace bpp.

#endif //_DF_PHYLOGENY_TREE_H_
