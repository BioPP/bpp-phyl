// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_SUBSTITUTIONPROCESS_H
#define BPP_PHYL_LIKELIHOOD_SUBSTITUTIONPROCESS_H


#include "../Model/SubstitutionModel.h"
#include "ModelScenario.h"
#include "ParametrizablePhyloTree.h"

// From bpp-core:
#include <Bpp/Numeric/ParameterAliasable.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

#include <Bpp/Seq/Container/AlignmentData.h>

// From the STL:
#include <memory>

namespace bpp
{
/**
 * @brief This interface describes the substitution process along the tree and sites of the alignment.
 *
 * It main purpose is to provide the necessary calculus for each branch-site-model class combination,
 * such as Markov generator and transition probabilities.
 * These are typically provided by a BranchModel class, applied in various combination along the
 * tree (eg non-homogeneous models) and alignment (eg partition models).
 * The so-called "model class" refers to mixture models.
 *
 * An instance of the SubstitutionProcess class is always associated to an instance of a ParametrizableTree.
 * The substitution process is in charge of computing the transition probabilities for each branch of the tree,
 * and to update them in case a parameter (including branch length) changes.
 *
 * As several branches and sites can share the same generator/transition probabilities, calling these values
 * and therefore performing the underlying calculation for each branch-site can result in a very unefficient
 * code. Therefore, "model iterators" are provided, to allow smarter loops over the model structure,
 * by minimizing the amount of computation to be done.
 * Such iterator allow to loop over branches and sites in a clever way, through both directions, but do not
 * perform any calculations.
 * These are achieved through calls to the corresponding SubstitutionProcess class.
 */
class SubstitutionProcessInterface :
  public virtual ParameterAliasable
{
public:
  virtual SubstitutionProcessInterface* clone() const = 0;

public:
  /**
   * @return The state map associated with the models of this process
   */

  virtual const StateMapInterface& stateMap() const = 0;

  virtual std::shared_ptr<const StateMapInterface> getStateMap() const = 0;

  virtual bool isCompatibleWith(const AlignmentDataInterface& data) const = 0;

  virtual const ParametrizablePhyloTree& parametrizablePhyloTree() const = 0;

  virtual std::shared_ptr<const ParametrizablePhyloTree> getParametrizablePhyloTree() const = 0;

  virtual size_t getNumberOfClasses() const = 0;

  virtual size_t getNumberOfStates() const = 0;

  /**
   * @return The current number of distinct models.
   */

  virtual size_t getNumberOfModels() const = 0;

  /**
   * @return The current indexes of models used in the process
   */

  virtual std::vector<size_t> getModelNumbers() const = 0;

  /**
   * @return the model with given index.
   */
  virtual const BranchModelInterface& model(size_t i) const = 0;

  /**
   * @return A shared pointer toward the model with given index.
   */
  virtual std::shared_ptr<const BranchModelInterface> getModel(size_t i) const = 0;

  /**
   * @brief Get the substitution model corresponding to a certain
   * branch, site pattern, and model class.
   *
   * @param nodeId The id of the node.
   * @param classIndex The model class index.
   * @return the model with given index.
   */
  virtual const BranchModelInterface& model(
      unsigned int nodeId,
      size_t classIndex) const = 0;

  /**
   * @brief Get the substitution model corresponding to a certain
   * branch, site pattern, and model class.
   *
   * @param nodeId The id of the node.
   * @param classIndex The model class index.
   * @return A shared pointer toward the model with given index.
   */
  virtual std::shared_ptr<const BranchModelInterface> getModel(
      unsigned int nodeId,
      size_t classIndex) const = 0;

  /**
   * @brief Get the Model Scenario associated with this process, in
   * case there are mixture models involved.
   *
   * When a mixture model is not included in the ModelScenario, it
   * is considered as non-mixed (transition probabilities are then
   * computed as mixture of submodel transition probalities).
   *
   * It returns 0 if there is no model path, which means that all
   * mixture models are considered as non-mixed.
   */
  virtual std::shared_ptr<const ModelScenario> getModelScenario() const = 0;

  /**
   * @brief Get a list of nodes id for which the given model is associated.
   *
   * @param i The index of the model in the set.
   * @return A vector with the ids of the node associated to this model.
   * @throw IndexOutOfBoundsException If the index is not valid.
   */

  virtual const std::vector<unsigned int> getNodesWithModel(size_t i) const = 0;

  /**
   * @brief Get the number of the model associated to a particular
   * node id.
   *
   * @param nodeId The id of the query node.
   * @return The number of the model associated to the given node.
   * @throw Exception If no model is found for this node.
   */
  virtual size_t getModelNumberForNode(unsigned int nodeId) const = 0;

  /**
   * @brief Get the model associated to a particular node id.
   *
   * @param nodeId The id of the query node.
   * @return A pointer toward the corresponding model.
   * @throw Exception If no model is found for this node.
   */
  virtual std::shared_ptr<const BranchModelInterface> getModelForNode(unsigned int nodeId) const = 0;

  /**
   * @brief Get the rate distribution
   * @throw NullPointerException if there is no associated rate distribution.
   */
  virtual const DiscreteDistributionInterface& rateDistribution() const = 0;

  /**
   * @brief Get the rate distribution
   * @throw NullPointerException if there is no associated rate distribution.
   */
  virtual DiscreteDistributionInterface& rateDistribution() = 0;

  /**
   * @brief Get a pointer to the rate distribution (or null if there
   * is no rate distribution).
   */
  virtual std::shared_ptr<const DiscreteDistributionInterface> getRateDistribution() const = 0;

  /**
   * @brief Get a pointer to the rate distribution (or null if there
   * is no rate distribution).
   */
  virtual std::shared_ptr<DiscreteDistributionInterface> getRateDistribution() = 0;


  /**
   * @brief Methods to retrieve the parameters of specific objects.
   */
  virtual ParameterList getSubstitutionModelParameters(bool independent) const = 0;

  virtual ParameterList getRateDistributionParameters(bool independent) const = 0;

  virtual ParameterList getRootFrequenciesParameters(bool independent) const = 0;

  virtual ParameterList getBranchLengthParameters(bool independent) const = 0;

  virtual ParameterList getNonDerivableParameters() const = 0;

  /**
   * @brief Get the values of the frequencies for each state in the
   * alphabet at the root node.
   *
   * For reversible models, these are the equilibrium frequencies.
   * For non-reversible models, these usually are distinct parameters.
   *
   * @return A vector with ancestral frequencies for each state in the alphabet;
   */
  virtual const std::vector<double>& getRootFrequencies() const = 0;

  /**
   * @return true if the process has parametrized root frequencies (non-stationary model)
   */
  virtual bool hasRootFrequencySet() const = 0;

  /**
   * @return The set of parametrized root frequencies.
   */
  virtual const FrequencySetInterface& rootFrequencySet() const = 0;

  /**
   * @return A pointer toward rhe set of parametrized root frequencies.
   */
  virtual std::shared_ptr<const FrequencySetInterface> getRootFrequencySet() const = 0;

  /**
   * @return The set of parametrized root frequencies.
   */
  virtual FrequencySetInterface& rootFrequencySet() = 0;

  /**
   * @return A pointer toward rhe set of parametrized root frequencies.
   */
  virtual std::shared_ptr<FrequencySetInterface> getRootFrequencySet() = 0;

  /**
   * @return The probability associated to the given model class.
   *
   * @param classIndex The model class index.
   */
  virtual double getProbabilityForModel(size_t classIndex) const = 0;

  /**
   * @return The vector of the probabilities associated to the
   * classes.
   */
  virtual Vdouble getClassProbabilities() const = 0;

  /**
   * @return The substitution rate associated to the given model class.
   *
   * @param classIndex The model class index.
   */
  virtual double getRateForModel(size_t classIndex) const = 0;

  /**
   * @brief Tell if the transition probabilities have changed after the last call to setParameters().
   * @return True if transition probabilities have changed.
   */

  // virtual bool transitionProbabilitiesHaveChanged() const = 0; Not sure we need that anymore...
};
} // end namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_SUBSTITUTIONPROCESS_H
