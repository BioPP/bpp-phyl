// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_ABSTRACTTREELIKELIHOOD_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_ABSTRACTTREELIKELIHOOD_H

#include <Bpp/Numeric/AbstractParametrizable.h>

#include "TreeLikelihood.h"
#include "../../Model/SubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Container/AlignmentData.h>

// From the STL:
#include <memory>

namespace bpp
{
/**
 * @brief Partial implementation of the TreeLikelihood interface.
 *
 * This class implements a few methods useful for most of likelihood
 * computation methods.
 *
 * It includes a tree_ and a data_ pointers.
 * This objects are owned by the class, and hence hard copied when cloning, and destroyed by the destructor.
 *
 * - The Parametrizable interface;
 * - The getTree() method;
 *
 * It also adds an abstract method for recursive computations.
 */
class AbstractTreeLikelihood :
  public virtual TreeLikelihoodInterface,
  public AbstractParametrizable
{
public:
  /**
   * @brief A very simple branch iterator.
   *
   * The constructor takes a vector of nodes id to iterate over.
   */
  class SimpleBranchIterator :
    public BranchIterator
  {
private:
    std::vector<int> nodesId_;
    size_t index_;

public:
    SimpleBranchIterator(const std::vector<int>& nodesId) :
      nodesId_(nodesId), index_(0) {}

public:
    int next()
    {
      if (!hasNext())
        throw Exception("AbstractTreeLikelihood::SimpleBranchIterator::next(). No more branch in the set.");
      return nodesId_[index_++];
    }

    bool hasNext() const { return index_ < nodesId_.size(); }
  };

  /**
   * @brief A very simple site iterator.
   *
   * This iterator loops over a continuous range of sites.
   * The ocnstructor takes as input the number of sites to iterate over,
   * and optionally an offset argument, specifying the index of the first site.
   */
  class SimpleSiteIterator :
    public SiteIterator
  {
private:
    size_t maxIndex_;
    size_t index_;
    size_t offset_;

public:
    SimpleSiteIterator(size_t nbSites, size_t offset = 0) :
      maxIndex_(nbSites), index_(0), offset_(offset) {}

public:
    size_t next()
    {
      if (!hasNext())
        throw Exception("AbstractTreeLikelihood::SimpleSiteIterator::next(). No more site in the set.");
      return offset_ + index_++;
    }

    bool hasNext() const { return index_ < maxIndex_; }
  };

  /**
   * @name Branch iterator for models without site partition.
   *
   * @{
   */
  class ConstNoPartitionBranchModelDescription :
    public ConstBranchModelDescription
  {
private:
    std::shared_ptr<const TransitionModelInterface> model_;
    size_t nbSites_;

public:
    ConstNoPartitionBranchModelDescription(std::shared_ptr<const TransitionModelInterface> model, size_t nbSites) :
      model_(model), nbSites_(nbSites) {}

    ConstNoPartitionBranchModelDescription(const ConstNoPartitionBranchModelDescription& bmd) :
      model_(bmd.model_),
      nbSites_(bmd.nbSites_)
    {}

    ConstNoPartitionBranchModelDescription& operator=(const ConstNoPartitionBranchModelDescription& bmd)
    {
      model_ = bmd.model_;
      nbSites_ = bmd.nbSites_;
      return *this;
    }

public:
    std::shared_ptr<const TransitionModelInterface> getModel() const override  { return model_; }
    const TransitionModelInterface& model() const override  { return *model_; }

    std::shared_ptr<const SubstitutionModelInterface> getSubstitutionModel() const override
    {
      return std::dynamic_pointer_cast<const SubstitutionModelInterface>(model_);
    }

    const SubstitutionModelInterface& substitutionModel() const override
    {
      return dynamic_cast<const SubstitutionModelInterface&>(*model_);
    }

    SiteIterator* getNewSiteIterator() const override { return new SimpleSiteIterator(nbSites_); }
  };

  class ConstNoPartitionBranchModelIterator :
    public ConstBranchModelIterator
  {
private:
    ConstNoPartitionBranchModelDescription branchModelDescription_;
    size_t index_;

public:
    ConstNoPartitionBranchModelIterator(std::shared_ptr<const TransitionModelInterface> model, size_t nbSites) :
      branchModelDescription_(model, nbSites), index_(0) {}

public:
    ConstNoPartitionBranchModelDescription* next()
    {
      if (!hasNext())
        throw Exception("AbstractHomogeneousTreeLikelihood::ConstHomogeneousBranchModelIterator::next(). No more branch in the set.");
      index_++;
      return &branchModelDescription_;
    }

    bool hasNext() const { return index_ == 0; }
  };

  class ConstNoPartitionSiteModelDescription :
    public ConstSiteModelDescription
  {
private:
    std::shared_ptr<const TransitionModelInterface> model_;
    std::vector<int> nodesId_;

public:
    ConstNoPartitionSiteModelDescription(std::shared_ptr<const TransitionModelInterface> model, const std::vector<int> nodesId) :
      model_(model), nodesId_(nodesId) {}

    ConstNoPartitionSiteModelDescription(const ConstNoPartitionSiteModelDescription& smd) :
      model_(smd.model_),
      nodesId_(smd.nodesId_)
    {}

    ConstNoPartitionSiteModelDescription& operator=(const ConstNoPartitionSiteModelDescription& smd)
    {
      model_ = smd.model_;
      nodesId_ = smd.nodesId_;
      return *this;
    }

public:
    std::shared_ptr<const TransitionModelInterface> getModel() const override { return model_; }
    const TransitionModelInterface& model() const override { return *model_; }

    std::shared_ptr<const SubstitutionModelInterface> getSubstitutionModel() const override
    {
      return std::dynamic_pointer_cast<const SubstitutionModelInterface>(model_);
    }

    const SubstitutionModelInterface& substitutionModel() const override
    {
      return dynamic_cast<const SubstitutionModelInterface&>(*model_);
    }

    BranchIterator* getNewBranchIterator() const override { return new SimpleBranchIterator(nodesId_); }
  };

  /** @} */

protected:
  std::unique_ptr<const AlignmentDataInterface> data_;
  mutable std::shared_ptr< TreeTemplate<Node> > tree_;
  bool computeFirstOrderDerivatives_;
  bool computeSecondOrderDerivatives_;
  bool initialized_;

public:
  AbstractTreeLikelihood() :
    AbstractParametrizable(""),
    data_(),
    tree_(),
    computeFirstOrderDerivatives_(true),
    computeSecondOrderDerivatives_(true),
    initialized_(false) {}

  AbstractTreeLikelihood(const AbstractTreeLikelihood& lik) :
    AbstractParametrizable(lik),
    data_(),
    tree_(),
    computeFirstOrderDerivatives_(lik.computeFirstOrderDerivatives_),
    computeSecondOrderDerivatives_(lik.computeSecondOrderDerivatives_),
    initialized_(lik.initialized_)
  {
    if (lik.data_) data_ = std::unique_ptr<AlignmentDataInterface>(lik.data_->clone());
    if (lik.tree_) tree_ = std::unique_ptr< TreeTemplate<Node> >(lik.tree_->clone());
  }

  AbstractTreeLikelihood& operator=(const AbstractTreeLikelihood& lik)
  {
    AbstractParametrizable::operator=(lik);
    if (lik.data_) data_ = std::unique_ptr<AlignmentDataInterface>(lik.data_->clone());
    else data_ = 0;
    if (lik.tree_) tree_ = std::unique_ptr< TreeTemplate<Node> >(lik.tree_->clone());
    else tree_ = 0;
    computeFirstOrderDerivatives_ = lik.computeFirstOrderDerivatives_;
    computeSecondOrderDerivatives_ = lik.computeSecondOrderDerivatives_;
    initialized_ = lik.initialized_;
    return *this;
  }

  /**
   * @brief Abstract class destructor
   *
   * This destructor is empty.
   */
  virtual ~AbstractTreeLikelihood() {}

public:
  /**
   * @name The TreeLikelihood interface.
   *
   * @{
   */
  bool hasLikelihoodData() const { return data_ != nullptr; }
  const AlignmentDataInterface& data() const { return *data_; }
  std::shared_ptr<const Alphabet> getAlphabet() const { return data_->getAlphabet(); }
  Vdouble getLikelihoodPerSite()                 const;
  Vdouble getLogLikelihoodPerSite()              const;
  VVdouble getLikelihoodPerSitePerState()    const;
  VVdouble getLogLikelihoodPerSitePerState() const;
  size_t getNumberOfSites() const { return data_->getNumberOfSites(); }
  const Tree& tree() const { return *tree_; }
  void enableDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
  void enableFirstOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = yn; }
  void enableSecondOrderDerivatives(bool yn) { computeFirstOrderDerivatives_ = computeSecondOrderDerivatives_ = yn; }
  bool enableFirstOrderDerivatives() const { return computeFirstOrderDerivatives_; }
  bool enableSecondOrderDerivatives() const { return computeSecondOrderDerivatives_; }
  bool isInitialized() const { return initialized_; }
  void initialize() { initialized_ = true; }
  /** @} */
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_ABSTRACTTREELIKELIHOOD_H
