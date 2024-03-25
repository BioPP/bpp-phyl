// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_TREELIKELIHOOD_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_TREELIKELIHOOD_H

#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/VectorTools.h>

#include "../../Model/SubstitutionModel.h"
#include "TreeLikelihoodData.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignmentData.h>

namespace bpp
{
/**
 * @brief The TreeLikelihood interface.
 *
 * This interface defines the methods needed for computing the likelihood
 * of a phylogenetic tree, given a dataset.
 */
class TreeLikelihoodInterface :
  public virtual SecondOrderDerivable
{
public:
  /**
   * @brief An iterator over a set of branches, specified by their node ids.
   */
  class BranchIterator
  {
public:
    virtual ~BranchIterator() {}

public:
    /**
     * @return The id of the next node in the set.
     */
    virtual int next() = 0;
    /**
     * @return True if there is at least another node in the set.
     */
    virtual bool hasNext() const = 0;
  };

  /**
   * @brief An iterator over a set of sites, specified by their position.
   *
   * In most cases, the position will reflect the index of an inner array used for likelihood storage.
   */
  class SiteIterator
  {
public:
    virtual ~SiteIterator() {}

public:
    /**
     * @return The position of the next site in the set.
     */
    virtual size_t next() = 0;
    /**
     * @return True is there is at least another site in the set.
     */
    virtual bool hasNext() const = 0;
  };

  /**
   * @brief A pair of SubstitutionModel / SiteIterator.
   */
  class ConstBranchModelDescription
  {
public:
    virtual ~ConstBranchModelDescription() {}

public:
    virtual const TransitionModelInterface& model() const = 0;
    virtual std::shared_ptr<const TransitionModelInterface> getModel() const = 0;

    virtual const SubstitutionModelInterface& substitutionModel() const = 0;
    virtual std::shared_ptr<const SubstitutionModelInterface> getSubstitutionModel() const = 0;

    virtual SiteIterator* getNewSiteIterator() const = 0;
  };

  /**
   * @brief Iterates through all models used for all sites on a given branch.
   */
  class ConstBranchModelIterator
  {
public:
    virtual ~ConstBranchModelIterator() {}

public:
    virtual ConstBranchModelDescription* next() = 0;
    virtual bool hasNext() const = 0;
  };

  /**
   * @brief A pair of SubstitutionModel / BranchIterator.
   */
  class ConstSiteModelDescription
  {
public:
    virtual ~ConstSiteModelDescription() {}

public:
    virtual const TransitionModelInterface& model() const = 0;
    virtual std::shared_ptr<const TransitionModelInterface> getModel() const = 0;

    virtual const SubstitutionModelInterface& substitutionModel() const = 0;
    virtual std::shared_ptr<const SubstitutionModelInterface> getSubstitutionModel() const = 0;

    virtual BranchIterator* getNewBranchIterator() const = 0;
  };

  /**
   * @brief Iterates through all models used for all branches on a given site.
   */
  class ConstSiteModelIterator
  {
public:
    virtual ~ConstSiteModelIterator() {}

public:
    virtual ConstSiteModelDescription* next() = 0;
    virtual bool hasNext() const = 0;
  };

public:
  TreeLikelihoodInterface() {}
  virtual ~TreeLikelihoodInterface() {}

  TreeLikelihoodInterface* clone() const override = 0;

public:
  /**
   * @brief Set the dataset for which the likelihood must be evaluated.
   *
   * @param sites The data set to use.
   */
  virtual void setData(const AlignmentDataInterface& sites) = 0;

  /**
   * @brief Get the dataset for which the likelihood must be evaluated.
   *
   * @return The site container where the sequences are stored.
   */
  virtual const AlignmentDataInterface& data() const = 0;

  /**
   * @brief Init the likelihood object.
   *
   * This method is used to initialize all parameters.
   * It is typically called after the constructor and the setData method.
   * It contains virtual methods that can't be called in the constructor.
   * @throw Exception if something bad happened, for instance if no data are associated to the likelihood function.
   */
  virtual void initialize() = 0;

  /**
   * @return 'true' is the likelihood function has been initialized.
   */
  virtual bool isInitialized() const = 0;

  /**
   * @return True if the TreeLikelihood object has a LikelihoodData.
   */
  virtual bool hasLikelihoodData() const = 0;

  /**
   * @return The underlying likelihood data structure.
   */
  virtual TreeLikelihoodData& likelihoodData() = 0;

  /**
   * @return The underlying likelihood data structure.
   */
  virtual const TreeLikelihoodData& likelihoodData() const = 0;

  /**
   * @brief Get the likelihood for a site.
   *
   * @param site The site index to analyse.
   * @return The likelihood for site <i>site</i>.
   */
  virtual double getLikelihoodForASite(size_t site) const = 0;

  /**
   * @brief Get the logarithm of the likelihood for a site.
   *
   * @param site The site index to analyse.
   * @return The logarithm of the likelihood for site <i>site</i>.
   */
  virtual double getLogLikelihoodForASite(size_t site) const = 0;

  /**
   * @brief Get the likelihood for a site and for a state.
   *
   * @param site The site index to analyse.
   * @param state The state to consider.
   * @return The likelihood for site <i>site</i> and state <i>state</i>.
   */
  virtual double getLikelihoodForASiteForAState(size_t site, int state) const = 0;

  /**
   * @brief Get the logarithm of the likelihood for a site and for a state.
   *
   * @param site The site index to analyse.
   * @param state The state to consider.
   * @return The logarithm of the likelihood for site <i>site</i> and state <i>state</i>.
   */
  virtual double getLogLikelihoodForASiteForAState(size_t site, int state) const = 0;

  /**
   * @brief Get the likelihood for each site.
   *
   * @return A vector with all likelihoods for each site.
   */
  virtual Vdouble getLikelihoodPerSite() const = 0;

  /**
   * @brief Get the logarithm of the likelihood for each site.
   *
   * @return A vector with all log likelihoods for each site.
   */
  virtual Vdouble getLogLikelihoodPerSite() const = 0;

  /**
   * @brief Get the likelihood for each site and for each state.
   *
   * @return A 2d vector with all likelihoods for each site and for each state.
   */
  virtual VVdouble getLikelihoodPerSitePerState() const = 0;

  /**
   * @brief Get the logarithm of the likelihood for each site and for each state.
   *
   * @return A 2d vector with all log likelihoods for each site and for each state.
   */
  virtual VVdouble getLogLikelihoodPerSitePerState() const = 0;

  /**
   * @brief Get the likelihood for the whole dataset.
   *
   * @return The likelihood of the dataset.
   */
  virtual double getLikelihood() const = 0;

  /**
   * @brief Get the logarithm of the likelihood for the whole dataset.
   *
   * @return The logarithm of the likelihood of the dataset.
   */
  virtual double getLogLikelihood() const = 0;

  /**
   * @brief Get the tree (topology and branch lengths).
   *
   * @return The tree of this TreeLikelihood object.
   */
  virtual const Tree& tree() const = 0;

  /**
   * @brief Get the number of sites in the dataset.
   *
   * @return the number of sites in the dataset.
   */
  virtual size_t getNumberOfSites() const = 0;

  /**
   * @return the number of model states of the underlying Markov chain.
   */
  virtual size_t getNumberOfStates() const = 0;

  /**
   * @return the alphabet state corresponding to the given model state.
   */
  virtual int getAlphabetStateAsInt(size_t i) const = 0;

  /**
   * @return the alphabet state corresponding to the given model state.
   */
  virtual std::string getAlphabetStateAsChar(size_t i) const = 0;

  /**
   * @return the alphabet states corresponding to all model states.
   */
  virtual const std::vector<int>& getAlphabetStates() const = 0;

  /**
   * @brief Get the alphabet associated to the dataset.
   *
   * @return the alphabet associated to the dataset.
   */
  virtual std::shared_ptr<const Alphabet> getAlphabet() const = 0;

  /**
   * @name Retrieve some particular parameters subsets.
   *
   * @{
   */

  /**
   * @brief Get the branch lengths parameters.
   *
   * @return A ParameterList with all branch lengths.
   */
  virtual ParameterList getBranchLengthsParameters() const = 0;

  /**
   * @brief Get the parameters associated to substitution model(s).
   *
   * @return A ParameterList.
   */
  virtual ParameterList getSubstitutionModelParameters() const = 0;

  /**
   * @brief Get the substitution model associated to a given node and alignment column.
   *
   * @param nodeId The id of the request node.
   * @param siteIndex The index of the alignment position.
   * @see getSiteIndex
   * @return A pointer toward the corresponding model.
   */
  virtual std::shared_ptr<const TransitionModelInterface> getModelForSite(int nodeId, size_t siteIndex) const = 0;

  /**
   * @brief Get the substitution model associated to a given node and alignment column.
   *
   * @param nodeId The id of the request node.
   * @param siteIndex The index of the alignment position.
   * @see getSiteIndex
   * @return A pointer toward the corresponding model.
   * @throw NodeNotFoundException This exception may be thrown if the node is not found (depending on the implementation).
   */
  virtual std::shared_ptr<TransitionModelInterface> getModelForSite(int nodeId, size_t siteIndex) = 0;

  /**
   * @brief Retrieves all Pij(t) for a particular branch, defined by the upper node and site.
   *
   * These intermediate results may be used by other methods.
   *
   * @param nodeId The node defining the branch of interest.
   * @param siteIndex The index of the alignment position.
   * @see getSiteIndex
   * @return An array of dimension 2, where a[x][y] is the probability of substituting from x to y.
   */
  virtual VVdouble getTransitionProbabilities(int nodeId, size_t siteIndex) const = 0;

  virtual ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const = 0;

  virtual ConstSiteModelIterator* getNewSiteModelIterator(size_t siteIndex) const = 0;

  /**
   * @brief Get the index (used for inner computations) of a given site (original alignment column).
   *
   * @param site An alignment position.
   * @return The site index corresponding to the given input alignment position.
   */
  virtual size_t getSiteIndex(size_t site) const = 0;

  /**
   * @brief Get the values of the frequencies for each state in the alphabet at the root node.
   *
   * For reversible models, these are the equilibrium frequencies.
   * For non-reversible models, these usually are distinct parameters.
   *
   * For models without site partitioning, the set of frequencies is the same for all positions.
   * For partition models, the frequencies may differ from one site to another.
   *
   * @param siteIndex The index of the alignment position.
   * @see getSiteIndex
   * @return A vector with ancestral frequencies for each state in the alphabet;
   */
  virtual const std::vector<double>& getRootFrequencies(size_t siteIndex) const = 0;

  /** @} */

  /**
   * @brief Tell if derivatives must be computed.
   *
   * This methods calls the enableFirstOrderDerivatives and enableSecondOrderDerivatives.
   *
   * @param yn Yes or no.
   */
  virtual void enableDerivatives(bool yn) = 0;

  /**
   * @brief All derivable parameters.
   *
   * Usually, this contains all branch lengths parameters.
   *
   * @return A ParameterList.
   */
  virtual ParameterList getDerivableParameters() const = 0;

  /**
   * @brief All non derivable parameters.
   *
   * Usually, this contains all substitution model parameters and rate distribution.
   *
   * @return A ParameterList.
   */
  virtual ParameterList getNonDerivableParameters() const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_TREELIKELIHOOD_H
