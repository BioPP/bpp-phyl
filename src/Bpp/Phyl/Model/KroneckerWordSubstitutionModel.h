// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_KRONECKERWORDSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_KRONECKERWORDSUBSTITUTIONMODEL_H


#include "AbstractKroneckerWordSubstitutionModel.h"

// From bpp-core
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/BppVector.h>

namespace bpp
{
/**
 * @brief Basal class for words of  substitution models with multiple
 * substitutions.
 * @author Laurent Guéguen
 *
 * The equilibrium frequency of each word is due to the combination of
 * the multiple substitutions.
 */
class KroneckerWordSubstitutionModel :
  public AbstractKroneckerWordSubstitutionModel
{
public:
  /**
   * @brief Build a new KroneckerWordSubstitutionModel object from a
   * Vector of pointers to SubstitutionModels.
   *
   * @param modelList the list of substitution models to use, in
   *   the order of the positions in the words from left to right. All
   *   the models must be different objects to avoid parameters
   *   redundancy, otherwise only the first model is used. The used models
   *   are owned by the instance.
   * @param prefix the Namespace.
   */
  KroneckerWordSubstitutionModel(ModelList& modelList, const std::string& prefix = "");

  /**
   * @brief Build a new KroneckerWordSubstitutionModel object from a
   * Vector of pointers to SubstitutionModels.
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
  KroneckerWordSubstitutionModel(ModelList& modelList,
                                 const std::vector<std::set< size_t> >& vPos,
                                 const std::string& prefix = "");

  /**
   * @brief Build a new KroneckerWordSubstitutionModel object from a
   * pointer to an SubstitutionModel and a number of
   * desired models.
   *
   * @param pmodel pointer to the substitution model to use in all the
   *  positions. It is owned by the instance.
   * @param num The number of models involved.
   * @param prefix the Namespace.
   */
  KroneckerWordSubstitutionModel(
      std::unique_ptr<SubstitutionModelInterface> pmodel,
      unsigned int num,
      const std::string& prefix = "");

  /**
   * @brief Build a new KroneckerWordSubstitutionModel object from a
   * pointer to an SubstitutionModel and a number of
   * desired models.
   *
   * @param pmodel pointer to the substitution model to use in all the
   *  positions. It is owned by the instance.
   * @param num The number of models involved.
   * @param vPos a vector of sets of simultaneously changing
   *   positions.
   * @param prefix the Namespace.
   */
  KroneckerWordSubstitutionModel(
      std::unique_ptr<SubstitutionModelInterface> pmodel,
      unsigned int num,
      const std::vector<std::set< size_t> >& vPos,
      const std::string& prefix = "");

  virtual ~KroneckerWordSubstitutionModel() {}

  KroneckerWordSubstitutionModel* clone() const override {
    return new KroneckerWordSubstitutionModel(*this);
  }

protected:

  /**
   * @brief Constructor for the derived classes only
   */
  KroneckerWordSubstitutionModel(
      std::shared_ptr<const Alphabet> alph,
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& prefix = "");

  void completeMatrices_() override {}

public:
  virtual std::string getName() const override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_KRONECKERWORDSUBSTITUTIONMODEL_H
