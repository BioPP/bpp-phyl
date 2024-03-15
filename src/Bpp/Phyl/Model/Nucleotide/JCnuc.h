// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_NUCLEOTIDE_JCNUC_H
#define BPP_PHYL_MODEL_NUCLEOTIDE_JCNUC_H


#include "../AbstractSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The Jukes-Cantor substitution model for nucleotides.
 *
 * All rates equal:
 * \f[
 * S = \begin{pmatrix}
 * \cdots & r & r & r \\
 * r & \cdots & r & r \\
 * r & r & \cdots & r \\
 * r & r & r & \cdots \\
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = diag\left(\frac{1}{4}, \frac{1}{4}, \frac{1}{4}, \frac{1}{4}\right)
 * \f]
 * Normalization: \f$r\f$ is set so that \f$\sum_i Q_{i,i}\pi_i = -1\f$:
 * \f[
 * S = \begin{pmatrix}
 * -4 & \frac{4}{3} & \frac{4}{3} & \frac{4}{3} \\
 * \frac{4}{3} & -4 & \frac{4}{3} & \frac{4}{3} \\
 * \frac{4}{3} & \frac{4}{3} & -4 & \frac{4}{3} \\
 * \frac{4}{3} & \frac{4}{3} & \frac{4}{3} & -4 \\
 * \end{pmatrix}
 * \f]
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$pi\f$:
 * \f[
 * Q = S . \pi = \begin{pmatrix}
 * -1 & \frac{1}{3} & \frac{1}{3} & \frac{1}{3} \\
 * \frac{1}{3} & -1 & \frac{1}{3} & \frac{1}{3} \\
 * \frac{1}{3} & \frac{1}{3} & -1 & \frac{1}{3} \\
 * \frac{1}{3} & \frac{1}{3} & \frac{1}{3} & -1 \\
 * \end{pmatrix}
 * \f]
 *
 * The eigen values are \f$\left(0, -\frac{4}{3}, -\frac{4}{3}, -\frac{4}{3}\right)\f$,
 * the left eigen vectors are, by row:
 * \f[
 * U = \begin{pmatrix}
 *  \frac{1}{4} &  \frac{1}{4} &  \frac{1}{4} &  \frac{1}{4} \\
 * -\frac{1}{4} & -\frac{1}{4} &  \frac{3}{4} & -\frac{1}{4} \\
 * -\frac{1}{4} &  \frac{3}{4} & -\frac{1}{4} & -\frac{1}{4} \\
 *  \frac{3}{4} & -\frac{1}{4} & -\frac{1}{4} & -\frac{1}{4} \\
 * \end{pmatrix}
 * \f]
 * and the right eigen vectors are, by column:
 * \f[
 * U^-1 = \begin{pmatrix}
 * 1 &  0 &  0 &  1 \\
 * 1 &  0 &  1 &  0 \\
 * 1 &  1 &  0 &  0 \\
 * 1 & -1 & -1 & -1 \\
 * \end{pmatrix}
 * \f]
 *
 * In addition, a rate_ factor defines the mean rate of the model.
 *
 * The probabilities of changes are computed analytically using the formulas:
 * \f[
 * P_{i,j}(t) = \begin{cases}
 * \frac{1}{4} + \frac{3}{4}e^{-rate\_*\frac{4}{3}t} & \text{if $i=j$}, \\
 * \frac{1}{4} - \frac{1}{4}e^{-rate\_*\frac{4}{3}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 *
 * First and second order derivatives are also computed analytically using the formulas:
 * \f[
 * \frac{\partial P_{i,j}(t)}{\partial t} = rate\_ * \begin{cases}
 * -e^{-rate\_*\frac{4}{3}t}           & \text{if $i=j$}, \\
 * \frac{1}{3}e^{-rate*\frac{4}{3}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 * \f[
 * \frac{\partial^2 P_{i,j}(t)}{\partial t^2} = rate\_^2 *\begin{cases}
 * \frac{4}{3}e^{-rate\_*\frac{4}{3}t}  & \text{if $i=j$}, \\
 * -\frac{4}{9}e^{-rate\_*\frac{4}{3}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 *
 * Reference:
 * - Jukes TH and Cantor CR (1969), Evolution_ of proteins molecules_, 121-123, in Mammalian_ protein metabolism_.
 */
class JCnuc :
  public AbstractReversibleNucleotideSubstitutionModel
{
private:
  mutable double exp_;
  mutable RowMatrix<double> p_;

public:
  JCnuc(std::shared_ptr<const NucleicAlphabet> alpha);

  virtual ~JCnuc() {}

  JCnuc* clone() const override { return new JCnuc(*this); }

public:
  double Pij_t    (size_t i, size_t j, double d) const override;
  double dPij_dt  (size_t i, size_t j, double d) const override;
  double d2Pij_dt2(size_t i, size_t j, double d) const override;
  const Matrix<double>& getPij_t    (double d) const override;
  const Matrix<double>& getdPij_dt  (double d) const override;
  const Matrix<double>& getd2Pij_dt2(double d) const override;

  std::string getName() const override { return "JC69"; }

  /**
   * @brief This method is disabled in this model since frequencies are not free parameters.
   *
   * Consider using the HKY85 model for instance if you want to set frequencies as parameters.
   *
   * @param data Unused parameter.
   * @param pseudoCount Unused parameter.
   */
  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override {}

protected:
  
  /**
   * In the case of the model of Jukes & Cantor, this method is not usefull since
   * the generator is fully determined. No matrice can be changed... This method is only
   * used in the constructor of the class.
   */
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_JCNUC_H
