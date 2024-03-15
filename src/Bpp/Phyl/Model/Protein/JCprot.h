// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_PROTEIN_JCPROT_H
#define BPP_PHYL_MODEL_PROTEIN_JCPROT_H

#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>

#include "../AbstractSubstitutionModel.h"
#include "../FrequencySet/ProteinFrequencySet.h"
#include "ProteinSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The Jukes-Cantor substitution model for proteins.
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
 * \pi = diag\left(\frac{1}{20}, \ldots, \frac{1}{20}\right)
 * \f]
 * Normalization: \f$r\f$ is set so that \f$\sum_i Q_{i,i}\pi_i = -1\f$:
 * \f[
 * S = \begin{pmatrix}
 * -20           & \frac{20}{19} & \ldots        & \frac{20}{19} \\
 * \frac{20}{19} &           -20 & \ddots        & \vdots \\
 * \vdots        & \ddots        & \ddots        & \frac{20}{19} \\
 * \frac{20}{19} & \ldots        & \frac{20}{19} & -20 \\
 * \end{pmatrix}
 * \f]
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$pi\f$:
 * \f[
 * Q = S . \pi = \begin{pmatrix}
 * -1 & \frac{1}{19}     & \ldots & \frac{1}{19} \\
 * \frac{1}{19} & -1     & \ddots & \vdots \\
 * \vdots       & \ddots & \ddots & \frac{1}{19} \\
 * \frac{1}{19} & \ldots & \frac{1}{19} & -1 \\
 * \end{pmatrix}
 * \f]
 *
 * The eigen values are \f$\left(0, -\frac{20}{19}, \ldots, -\frac{20}{19}\right)\f$,
 * the left eigen vectors are, by row:
 * \f[
 * U = \begin{pmatrix}
 *   \frac{1}{20} &        \ldots &  \frac{1}{20} &   \frac{1}{20} &  \frac{1}{20} \\
 *  -\frac{1}{20} &        \ldots & -\frac{1}{20} &  \frac{19}{20} & -\frac{1}{20} \\
 *         \vdots &        \ddots & \frac{19}{20} & -\frac{1}{20}  & -\frac{1}{20} \\
 *  -\frac{1}{20} &        \ddots &        \ddots &         \vdots &        \vdots \\
 *  \frac{19}{20} & -\frac{1}{20} &        \ldots &  -\frac{1}{20} & -\frac{1}{20} \\
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
 * \frac{1}{20} + \frac{19}{20}e^{- rate\_ * \frac{20}{19}t}& \text{if $i=j$}, \\
 * \frac{1}{20} - \frac{1}{20}e^{- rate\_ * \frac{20}{19}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 *
 * First and second order derivatives are also computed analytically using the formulas:
 * \f[
 * \frac{\partial P_{i,j}(t)}{\partial t} = rate\_ * \begin{cases}
 * -e^{- rate\_ * \frac{20}{19}t}           & \text{if $i=j$}, \\
 * \frac{1}{19}e^{- rate\_ * \frac{20}{19}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 * \f[
 * \frac{\partial^2 P_{i,j}(t)}{\partial t^2} =  rate\_^2 * \begin{cases}
 * \frac{20}{19}e^{- rate\_ * \frac{20}{19}t}  & \text{if $i=j$}, \\
 * -\frac{20}{361}e^{- rate\_ * \frac{20}{19}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 *
 * Reference:
 * - Jukes TH and Cantor CR (1969), Evolution_ of proteins molecules_, 121-123, in Mammalian_ protein metabolism_.
 */
class JCprot :
  public AbstractReversibleProteinSubstitutionModel
{
private:
  mutable double exp_;
  mutable RowMatrix<double> p_;
  std::unique_ptr<ProteinFrequencySetInterface> freqSet_;
  bool withFreq_;

public:
  /**
   * @brief Build a simple JC69 model, with original equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   */
  JCprot(std::shared_ptr<const ProteicAlphabet> alpha);

  /**
   * @brief Build a JC69 model with special equilibrium frequencies.
   *
   * @param alpha A proteic alphabet.
   * @param freqSet A pointer toward a protein frequencies set, which will be owned by this instance.
   * @param initFreqs Tell if the frequency set should be initialized with the original JTT92 values.
   * Otherwise, the values of the set will be used.
   */
  JCprot(
      std::shared_ptr<const ProteicAlphabet> alpha,
      std::unique_ptr<ProteinFrequencySetInterface> freqSet,
      bool initFreqs = false);

  JCprot(const JCprot& model) :
    AbstractParameterAliasable(model),
    AbstractReversibleProteinSubstitutionModel(model),
    exp_(model.exp_),
    p_(model.p_),
    freqSet_(model.freqSet_->clone()),
    withFreq_(model.withFreq_)
  {}

  JCprot& operator=(const JCprot& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractReversibleProteinSubstitutionModel::operator=(model);
    exp_ = model.exp_;
    p_   = model.p_;
    freqSet_.reset(model.freqSet_->clone());
    withFreq_ = model.withFreq_;
    return *this;
  }

  virtual ~JCprot() {}

  JCprot* clone() const override { return new JCprot(*this); }

public:
  double Pij_t    (size_t i, size_t j, double d) const override;
  double dPij_dt  (size_t i, size_t j, double d) const override;
  double d2Pij_dt2(size_t i, size_t j, double d) const override;
  const Matrix<double>& getPij_t    (double d) const override;
  const Matrix<double>& getdPij_dt  (double d) const override;
  const Matrix<double>& getd2Pij_dt2(double d) const override;

  std::string getName() const override
  {
    return withFreq_ ? "JC69+F" : "JC69";
  }

  void fireParameterChanged(const ParameterList& parameters) override
  {
    freqSet_->matchParametersValues(parameters);
    freq_ = freqSet_->getFrequencies();
    AbstractReversibleSubstitutionModel::fireParameterChanged(parameters);
  }

  void setFrequencySet(const ProteinFrequencySetInterface& freqSet)
  {
    freqSet_.reset(freqSet.clone());
    resetParameters_();
    addParameters_(freqSet_->getParameters());
  }

  const FrequencySetInterface& frequencySet() const override
  {
    if (freqSet_)
      return *freqSet_;
    throw NullPointerException("JCprot::frequencySet(). No associated FrequencySet.");
  }
    
  void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix);
    freqSet_->setNamespace(prefix + freqSet_->getName() + ".");
  }

  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override;

protected:
  
  /**
   * In the case of the model of Jukes & Cantor, this method is useless since
   * the generator is fixed! No matrice can be changed... This method is only
   * used in the constructor of the class.
   */
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_JCPROT_H
