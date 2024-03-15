// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_WORDSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_WORDSUBSTITUTIONMODEL_H


#include "AbstractWordSubstitutionModel.h"

// From bpp-core
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/BppVector.h>

namespace bpp
{
/**
 * @brief Basal class for words of substitution models.
 * @author Laurent Guéguen
 *
 * Only substitutions with one letter changed are accepted. Hence the
 * equilibrium frequency of each word is the product of the
 * equilibrium frequencies of the letters.</p>
 *
 * If there are @f$n@f$ models, @f$\rho_i@f$ is the rate of
 * model i (@f$\sum_{i=1}^{n} \rho_i = 1@f$) and the rates
 * are defined by relative rates parameters @f$r_i@f$
 * (called "relratei") with:
 * @f[
 * 1 <= i < n, \rho_i = (1-r_1).(1-r_2)...(1-r_{i-1}).r_{i}
 * @f]
 * @f[
 * \rho_n = (1-r_1).(1-r_2)...(1-r_{n-1})
 * @f]
 * and
 * @f[
 * \forall 1 <= i < n, r_i = \frac{\rho_i}{1-(\rho_1+...\rho_{i-1})}
 * @f]
 * where @f$\rho_i@f$ stands for the rate of position @f$i@f$.
 */
class WordSubstitutionModel :
  public AbstractWordSubstitutionModel
{
public:
  /**
   * @brief Build a new WordSubstitutionModel object from a
   * Vector of pointers to SubstitutionModels.
   *
   * @param modelList the list of substitution models to use, in
   *   the order of the positions in the words from left to right. All
   *   the models must be different objects to avoid parameters
   *   redundancy, otherwise only the first model is used. The used models
   *   are owned by the instance.
   * @param prefix the Namespace.
   */
  WordSubstitutionModel(ModelList& modelList, const std::string& prefix = "");

  /**
   * @brief Build a new WordSubstitutionModel object from a
   * pointer to an SubstitutionModel and a number of
   * desired models.
   *
   * @param pmodel pointer to the substitution model to use in all the
   *  positions. It is owned by the instance.
   * @param num The number of models involved.
   * @param prefix the Namespace.
   */
  WordSubstitutionModel(
      std::unique_ptr<SubstitutionModelInterface> pmodel,
      unsigned int num,
      const std::string& prefix = "");

  virtual ~WordSubstitutionModel() {}

  WordSubstitutionModel* clone() const override { return new WordSubstitutionModel(*this); }

protected:
  /**
   * @brief Constructor for the derived classes only
   */
  WordSubstitutionModel(
      std::shared_ptr<const Alphabet> alphabet,
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& prefix = "");

  virtual void updateMatrices_() override;

  virtual void completeMatrices_() override;

public:
  virtual const RowMatrix<double>& getPij_t(double d) const override;

  virtual const RowMatrix<double>& getdPij_dt(double d) const override;

  virtual const RowMatrix<double>& getd2Pij_dt2(double d) const override;

  virtual std::string getName() const override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_WORDSUBSTITUTIONMODEL_H
