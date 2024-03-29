// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H


#include "AbstractDiscreteRatesAcrossSitesTreeLikelihood.h"
#include "HomogeneousTreeLikelihood.h"

// From STL:
#include <memory>

namespace bpp
{
/**
 * @brief Partial implementation for homogeneous model of the TreeLikelihood interface.
 *
 * This class provides a pointer toward a single substitution model + several utilitary variables.
 */
class AbstractHomogeneousTreeLikelihood :
  public virtual HomogeneousTreeLikelihood,
  public AbstractDiscreteRatesAcrossSitesTreeLikelihood
{
public:
  class ConstHomogeneousSiteModelIterator :
    public ConstSiteModelIterator
  {
private:
    ConstNoPartitionSiteModelDescription siteModelDescription_;
    size_t index_;

public:
    ConstHomogeneousSiteModelIterator(const Tree& tree, std::shared_ptr<const TransitionModelInterface> model) :
      siteModelDescription_(model, tree.getBranchesId()),
      index_(0) {}

public:
    ConstSiteModelDescription* next()
    {
      if (!hasNext())
        throw Exception("AbstractHomogeneousTreeLikelihood::ConstHomogeneousSiteModelIterator::next(). No more site in the set.");
      index_++;
      return &siteModelDescription_;
    }

    bool hasNext() const { return index_ == 0; }
  };

protected:
  std::shared_ptr<TransitionModelInterface> model_;
  ParameterList brLenParameters_;

  mutable std::map<int, VVVdouble> pxy_;

  mutable std::map<int, VVVdouble> dpxy_;

  mutable std::map<int, VVVdouble> d2pxy_;

  std::vector<double> rootFreqs_;

  /**
   * @brief Pointer toward all nodes in the tree.
   *
   * The position in the array is the number used in the parameter name.
   * This may be different from the node id, unless you used the resetNodeId method on the input tree.
   */
  std::vector<Node*> nodes_;

  // some values we'll need:
  size_t nbSites_,         // the number of sites in the container
      nbDistinctSites_, // the number of distinct sites
      nbClasses_,       // the number of rate classes
      nbStates_,        // the number of states in the alphabet
      nbNodes_;         // the number of nodes in the tree

  bool verbose_;

  double minimumBrLen_;
  double maximumBrLen_;
  std::shared_ptr<ConstraintInterface> brLenConstraint_;

public:
  AbstractHomogeneousTreeLikelihood(
      const Tree& tree,
      std::shared_ptr<TransitionModelInterface> model,
      std::shared_ptr<DiscreteDistributionInterface> rDist,
      bool checkRooted = true,
      bool verbose = true);

  /**
   * @brief Copy constructor
   *
   * This constructor is to be called by the derived class copy constructor.
   */
  AbstractHomogeneousTreeLikelihood(const AbstractHomogeneousTreeLikelihood& lik);

  /**
   * @brief Assignation operator
   *
   * This operator is to be called by the derived class operator.
   */
  AbstractHomogeneousTreeLikelihood& operator=(const AbstractHomogeneousTreeLikelihood& lik);

  virtual ~AbstractHomogeneousTreeLikelihood() {}

private:
  /**
   * @brief Method called by constructor.
   */
  void init_(const Tree& tree,
      std::shared_ptr<TransitionModelInterface> model,
      std::shared_ptr<DiscreteDistributionInterface> rDist,
      bool checkRooted,
      bool verbose);

public:
  /**
   * @name The TreeLikelihood interface.
   *
   * Other methods are implemented in the AbstractTreeLikelihood class.
   *
   * @{
   */
  size_t getNumberOfStates() const { return model_->getNumberOfStates(); }

  const std::vector<int>& getAlphabetStates() const { return model_->getAlphabetStates(); }

  int getAlphabetStateAsInt(size_t i) const { return model_->getAlphabetStateAsInt(i); }

  std::string getAlphabetStateAsChar(size_t i) const { return model_->getAlphabetStateAsChar(i); }

  void initialize();

  ParameterList getBranchLengthsParameters() const;

  ParameterList getSubstitutionModelParameters() const;

  ParameterList getRateDistributionParameters() const
  {
    return AbstractDiscreteRatesAcrossSitesTreeLikelihood::getRateDistributionParameters();
  }

  const std::vector<double>& getRootFrequencies(size_t siteIndex) const { return model_->getFrequencies(); }

  VVVdouble getTransitionProbabilitiesPerRateClass(int nodeId, size_t siteIndex) const { return pxy_[nodeId]; }

  ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const
  {
    return new ConstNoPartitionBranchModelIterator(model_, nbDistinctSites_);
  }

  ConstSiteModelIterator* getNewSiteModelIterator(size_t siteIndex) const
  {
    return new ConstHomogeneousSiteModelIterator(*tree_, model_);
  }

  /** @} */

  /**
   * @name The HomogeneousTreeLikelihood interface.
   *
   * Other methods are implemented in the AbstractTreeLikelihood class.
   *
   * @{
   */
  std::shared_ptr<const TransitionModelInterface> getModel() const { return model_; }
  std::shared_ptr<const TransitionModelInterface> getModel(int nodeId, size_t siteIndex) const { return model_; }

  std::shared_ptr<TransitionModelInterface> getModel() { return model_; }
  std::shared_ptr<TransitionModelInterface> getModel(int nodeId, size_t siteIndex) { return model_; }

  void setModel(std::shared_ptr<TransitionModelInterface> model);
  /** @} */

  /**
   * @brief Get a SubstitutionModel pointer toward the model associated to this instance, if possible.
   *
   * Performs a cast operation on the pointer. Return NULL if cast failed.
   * @return A SubstitutionModel pointer toward the model associated to this instance.
   */
  virtual std::shared_ptr<const SubstitutionModelInterface> getSubstitutionModel() const
  {
    return std::dynamic_pointer_cast<const SubstitutionModelInterface>(model_);
  }

  /**
   * @brief Get a SubstitutionModel pointer toward the model associated to this instance, if possible.
   *
   * Performs a cast operation on the pointer. Return NULL if cast failed.
   * @return A SubstitutionModel pointer toward the model associated to this instance.
   *
   * @param nodeId Id of the node
   * @param siteIndex Position of the site
   */
  virtual std::shared_ptr<const SubstitutionModelInterface> getSubstitutionModel(int nodeId, size_t siteIndex) const
  {
    return std::dynamic_pointer_cast<const SubstitutionModelInterface>(model_);
  }

public:
  // Specific methods:

  /**
   * @brief This builds the <i>parameters</i> list from all parametrizable objects,
   * <i>i.e.</i> substitution model, rate distribution and tree.
   */
  virtual void initParameters();

  /**
   * @brief All parameters are stored in a parameter list.
   * This function apply these parameters to the substitution model,
   * to the rate distribution and to the branch lengths.
   */
  virtual void applyParameters();

  virtual void initBranchLengthsParameters(bool verbose = true);

  virtual void setMinimumBranchLength(double minimum)
  {
    if (minimum > maximumBrLen_)
      throw Exception("AbstractHomogeneousTreeLikelihood::setMinimumBranchLength. Minimum branch length sould be lower than the maximum one: " + TextTools::toString(maximumBrLen_));
    minimumBrLen_ = minimum;
    brLenConstraint_ = std::make_shared<IntervalConstraint>(minimumBrLen_, maximumBrLen_, true, true);
    initBranchLengthsParameters();
  }

  virtual void setMaximumBranchLength(double maximum)
  {
    if (maximum < minimumBrLen_)
      throw Exception("AbstractHomogeneousTreeLikelihood::setMaximumBranchLength. Maximum branch length sould be higher than the minimum one: " + TextTools::toString(minimumBrLen_));
    maximumBrLen_ = maximum;
    brLenConstraint_ = std::make_shared<IntervalConstraint>(minimumBrLen_, maximumBrLen_, true, true);
    initBranchLengthsParameters();
  }

  virtual double getMinimumBranchLength() const { return minimumBrLen_; }
  virtual double getMaximumBranchLength() const { return maximumBrLen_; }

protected:
  /**
   * @brief Fill the pxy_, dpxy_ and d2pxy_ arrays for all nodes.
   */
  virtual void computeAllTransitionProbabilities();
  /**
   * @brief Fill the pxy_, dpxy_ and d2pxy_ arrays for one node.
   */
  virtual void computeTransitionProbabilitiesForNode(const Node* node);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H
