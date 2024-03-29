// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_EQUIPROBABLESUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_EQUIPROBABLESUBSTITUTIONMODEL_H

#include <Bpp/Seq/Alphabet/Alphabet.h>

#include "AbstractSubstitutionModel.h"
#include "FrequencySet/FrequencySet.h"

namespace bpp
{
/**
 * @brief The EquiprobableSubstitutionModel substitution model for any kind of
 * alphabet. Jukes-Cantor models are specific case of this model,
 * applied to nucleotides or proteins.
 *
 * If the number of states is \f$N\f$.
 *
 * All rates equal:
 * \f[
 * \begin{pmatrix}
 * \ddots & r      & \ldots & r \\
 * r      & \ddots & \ddots & \vdots \\
 * \vdots & \ddots & \ddots & r \\
 * r      & \ldots & r      & \ddots \\
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = diag\left(\frac{1}{N}, \ldots, \frac{1}{N}\right)
 * \f]
 * Normalization: \f$r\f$ is set so that \f$\sum_i Q_{i,i}\pi_i = -1\f$:
 * \f[
 * S = \begin{pmatrix}
 * -N           & \frac{N}{N-1} & \ldots        & \frac{N}{N-1} \\
 * \frac{N}{N-1} &           -N & \ddots        & \vdots \\
 * \vdots        & \ddots        & \ddots        & \frac{N}{N-1} \\
 * \frac{N}{N-1} & \ldots        & \frac{N}{N-1} & -N \\
 * \end{pmatrix}
 * \f]
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$pi\f$:
 * \f[
 * Q = S . \pi = \begin{pmatrix}
 * -1 & \frac{1}{N-1}     & \ldots & \frac{1}{N-1} \\
 * \frac{1}{N-1} & -1     & \ddots & \vdots \\
 * \vdots       & \ddots & \ddots & \frac{1}{N-1} \\
 * \frac{1}{N-1} & \ldots & \frac{1}{N-1} & -1 \\
 * \end{pmatrix}
 * \f]
 *
 * The eigen values are \f$\left(0, -\frac{N}{N-1}, \ldots, -\frac{N}{N-1}\right)\f$,
 * the left eigen vectors are, by row:
 * \f[
 * U = \begin{pmatrix}
 *   \frac{1}{N} &        \ldots &  \frac{1}{N} &   \frac{1}{N} &  \frac{1}{N} \\
 *  -\frac{1}{N} &        \ldots & -\frac{1}{N} &  \frac{N-1}{N} & -\frac{1}{N} \\
 *         \vdots &        \ddots & \frac{N-1}{N} & -\frac{1}{N}  & -\frac{1}{N} \\
 *  -\frac{1}{N} &        \ddots &        \ddots &         \vdots &        \vdots \\
 *  \frac{N-1}{N} & -\frac{1}{N} &        \ldots &  -\frac{1}{N} & -\frac{1}{N} \\
 * \end{pmatrix}
 * \f]
 * and the right eigen vectors are, by column:
 * \f[
 * U^-1 = \begin{pmatrix}
 *      1 &      0 &  \ldots &      0 &  1 \\
 * \vdots & \vdots &  \ddots & \ddots &  0 \\
 *      1 &      0 &       1 & \ddots & \vdots \\
 *      1 &      1 &       0 & \ldots &  0 \\
 *      1 &     -1 &      -1 & \ldots & -1 \\
 * \end{pmatrix}
 * \f]
 *
 * In addition, a rate_ factor defines the mean rate of the model.
 *
 * The probabilities of changes are computed analytically using the formulas:
 * \f[
 * P_{i,j}(t) = \begin{cases}
 * \frac{1}{N} + \frac{N-1}{N}e^{- rate\_ * \frac{N}{N-1}t}& \text{if $i=j$}, \\
 * \frac{1}{N} - \frac{1}{N}e^{- rate\_ * \frac{N}{N-1}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 *
 * First and second order derivatives are also computed analytically using the formulas:
 * \f[
 * \frac{\partial P_{i,j}(t)}{\partial t} = rate\_ * \begin{cases}
 * -e^{- rate\_ * \frac{N}{N-1}t}           & \text{if $i=j$}, \\
 * \frac{1}{N-1}e^{- rate\_ * \frac{N}{N-1}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 * \f[
 * \frac{\partial^2 P_{i,j}(t)}{\partial t^2} =  rate\_^2 * \begin{cases}
 * \frac{N}{N-1}e^{- rate\_ * \frac{N}{N-1}t}  & \text{if $i=j$}, \\
 * -\frac{N}{(N-1)^2}e^{- rate\_ * \frac{N}{N-1}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 *
 */
class EquiprobableSubstitutionModel :
  public AbstractReversibleSubstitutionModel
{
private:
  mutable double exp_;
  mutable RowMatrix<double> p_;
  std::unique_ptr<FrequencySetInterface> freqSet_;

public:
  /**
   * @brief Build a simple equiprobable model, with original equilibrium frequencies.
   *
   * @param alpha An alphabet.
   */
  EquiprobableSubstitutionModel(std::shared_ptr<const Alphabet> alpha);

  /**
   * @brief Build an equiprobable model with special equilibrium frequencies.
   *
   * @param alpha An alphabet.
   * @param freqSet A pointer toward a frequencies set, which will be owned by this instance.
   * @param initFreqs Tell if the frequency set should be initialized with the original JTT92 values.
   * Otherwise, the values of the set will be used.
   */
  EquiprobableSubstitutionModel(
      std::shared_ptr<const Alphabet> alpha,
      std::unique_ptr<FrequencySetInterface> freqSet,
      bool initFreqs = false);

  EquiprobableSubstitutionModel(const EquiprobableSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    AbstractReversibleSubstitutionModel(model),
    exp_(model.exp_),
    p_(model.p_),
    freqSet_(model.freqSet_->clone())
  {}

  EquiprobableSubstitutionModel& operator=(const EquiprobableSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractReversibleSubstitutionModel::operator=(model);
    exp_ = model.exp_;
    p_   = model.p_;
    freqSet_.reset(model.freqSet_->clone());
    return *this;
  }

  virtual ~EquiprobableSubstitutionModel() {}

  EquiprobableSubstitutionModel* clone() const override { return new EquiprobableSubstitutionModel(*this); }

public:
  double Pij_t    (size_t i, size_t j, double d) const override;
  double dPij_dt  (size_t i, size_t j, double d) const override;
  double d2Pij_dt2(size_t i, size_t j, double d) const override;
  const Matrix<double>& getPij_t    (double d) const override;
  const Matrix<double>& getdPij_dt  (double d) const override;
  const Matrix<double>& getd2Pij_dt2(double d) const override;

  std::string getName() const override
  {
    if (freqSet_->getNamespace().find("+F.") != std::string::npos)
      return "Equi+F";
    else
      return "Equi";
  }

  void fireParameterChanged(const ParameterList& parameters) override
  {
    freqSet_->matchParametersValues(parameters);
    freq_ = freqSet_->getFrequencies();
    AbstractReversibleSubstitutionModel::fireParameterChanged(parameters);
  }

  void setFrequencySet(const FrequencySetInterface& freqSet)
  {
    freqSet_.reset(freqSet.clone());
    resetParameters_();
    addParameters_(freqSet_->getParameters());
  }

  const FrequencySetInterface& frequencySet() const override { return *freqSet_; }

  void setFreq(std::map<int, double>& freq) override;

protected:
  /**
   * In the case of the model of Jukes & Cantor, this method is useless since
   * the generator is fixed! No matrice can be changed... This method is only
   * used in the constructor of the class.
   */
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_EQUIPROBABLESUBSTITUTIONMODEL_H
