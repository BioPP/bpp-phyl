// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_ABSTRACTKRONECKERWORDSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_ABSTRACTKRONECKERWORDSUBSTITUTIONMODEL_H

#include <set>

#include "AbstractWordSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Abstract Kronecker Word Model.
 *
 * Objects of this class are built from several substitution models.
 * Each model corresponds to a position in the word. No model is
 * directly accessible.
 *
 * A substitution rate between words is the product of the rates of
 * changing letters of position generators.
 *
 * Either all substitutions are allowed, or any subset of them is
 * allowed. For example, in case of 3 words model, if @f$Q@f$ is the
 * generator of the position model, either all substitutions are
 * allowed (ie generator is @f$ Q_1 \otimes Q_1 \otimes Q_1 @f$ with
 * @f$Q_1@f$ is @f$Q@f$ with diagonal of 1), or all and only three
 * letters substitutions are allowed (ie generator is @f$ Q_0
 * \otimes Q_0 \otimes Q_0 @f$ with @f$Q_0@f$ is @f$Q@f$ with
 * diagonal of 0), or all and only two letters substitutions are
 * allowed (generator is @f$ I \otimes Q_0 \otimes Q_0 + Q_0 \otimes
 * I \otimes Q_0 + Q_0 \otimes Q_0 \otimes I @f$), or all and only
 * two neighbors letters substitutions are allowed (generator is @f$
 * I \otimes Q_0 \otimes Q_0 + Q_0 \otimes Q_0 \otimes I @f$).
 *
 * Default construction is all substitutions, but a vector of
 * vectors of simultaneously changing positions can be given.
 *
 */
class AbstractKroneckerWordSubstitutionModel :
  public AbstractWordSubstitutionModel
{
private:
  /**
   * @brief vector of sets of simultaneously changing positions.
   *
   * Shortcut : if empty, any number of any positions can change
   * simultaneously
   */
  std::vector<std::set< size_t>> sChangingPos_;

  /**
   * @brief vector of generators for computation purposes
   */
  std::vector< RowMatrix<double>> vGenerators_;

protected:
  void initGenerators_();

  /**
   * @brief checks that the vector of changing positions is valid
   */
  bool checkChangingPositions_();

private:
  /**
   * @brief First fill of the generator, from the position model
   */
  void fillBasicGenerator_();

public:
  /**
   * @brief Build a new AbstractKroneckerWordSubstitutionModel
   * object from a vector of pointers to SubstitutionModels,
   * allowing for all substitutions.
   *
   * @param modelList the list of substitution models to use, in
   *   the order of the positions in the words from left to right. All
   *   the models must be different objects to avoid parameters
   *   redundancy, otherwise only the first model is used. The used models
   *   are owned by the instance.
   * @param prefix the Namespace.
   */
  AbstractKroneckerWordSubstitutionModel(
      ModelList& modelList,
      const std::string& prefix);

  /**
   * @brief Build a new AbstractKroneckerWordSubstitutionModel
   * object from a vector of pointers to SubstitutionModels,
   * allowing for a given vector of simultaneous substitutions.
   *
   * @param modelList the list of substitution models to use, in
   *   the order of the positions in the words from left to right. All
   *   the models must be different objects to avoid parameters
   *   redundancy, otherwise only the first model is used. The used models
   *   are owned by the instance.
   * @param vPos a vector of sets of simultaneously changing
   *   positions.
   * @param prefix the Namespace.
   */
  AbstractKroneckerWordSubstitutionModel(
      ModelList& modelList,
      const std::vector<std::set< size_t>>& vPos,
      const std::string& prefix);

  /**
   * @brief Build a new AbstractWordSubstitutionModel object from a
   * pointer to an SubstitutionModel and a number of desired models,
   * allowing for all substitutions.
   *
   * @param pmodel A pointer to the substitution model to use in all
   * the positions. It will be owned by the instance.
   * @param num The number of models involved.
   * @param prefix the Namespace.
   */
  AbstractKroneckerWordSubstitutionModel(
      std::unique_ptr<SubstitutionModelInterface> pmodel,
      unsigned int num,
      const std::string& prefix);

  /**
   * @brief Build a new AbstractWordSubstitutionModel object from a
   * pointer to an SubstitutionModel and a number of desired models,
   * allowing for a given vector of simultaneous substitutions.
   *
   * @param pmodel A pointer to the substitution model to use in all
   * the positions. It will be owned by the instance.
   * @param num The number of models involved.
   * @param vPos a vector of vectors of simultaneously changing
   *   positions.
   * @param prefix the Namespace.
   */
  AbstractKroneckerWordSubstitutionModel(
      std::unique_ptr<SubstitutionModelInterface> pmodel,
      unsigned int num,
      const std::vector<std::set< size_t>>& vPos,
      const std::string& prefix);


  virtual ~AbstractKroneckerWordSubstitutionModel() {}

  AbstractKroneckerWordSubstitutionModel(const AbstractKroneckerWordSubstitutionModel&);

  AbstractKroneckerWordSubstitutionModel& operator=(const AbstractKroneckerWordSubstitutionModel&);

protected:
  /**
   * @brief Constructor for the derived classes only
   */
  AbstractKroneckerWordSubstitutionModel(
      std::shared_ptr<const Alphabet> alph,
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& prefix);

  void setChangingPositions(const std::vector<std::set<size_t>>& vPos);
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ABSTRACTKRONECKERWORDSUBSTITUTIONMODEL_H
