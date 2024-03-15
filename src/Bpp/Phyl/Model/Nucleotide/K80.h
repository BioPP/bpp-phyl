// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_NUCLEOTIDE_K80_H
#define BPP_PHYL_MODEL_NUCLEOTIDE_K80_H


#include "../AbstractSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief The Kimura 2-rates substitution model for nucleotides.
 *
 * All two rates: one for transitions and one for transversions.
 * This models include one parameter, the transition / transversion
 * relative rate \f$\kappa\f$.
 * \f[
 * S = \begin{pmatrix}
 * \cdots & r & \kappa r & r \\
 * r & \cdots & r & \kappa r \\
 * \kappa r & r & \cdots & r \\
 * r & \kappa r & r & \cdots \\
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = \left(\frac{1}{4}, \frac{1}{4}, \frac{1}{4}, \frac{1}{4}\right)
 * \f]
 * Normalization: \f$r\f$ is set so that \f$\sum_i Q_{i,i}\pi_i = -1\f$:
 * \f[
 * S = \begin{pmatrix}
 * -4 & \frac{4}{\kappa+2} & \frac{4\kappa}{\kappa+2} & \frac{4}{\kappa+2} \\
 * \frac{4}{\kappa+2} & -4 & \frac{4}{\kappa+2} & \frac{4\kappa}{\kappa+2} \\
 * \frac{4\kappa}{\kappa+2} & \frac{4}{\kappa+2} & -4 & \frac{4}{\kappa+2} \\
 * \frac{4}{\kappa+2} & \frac{4\kappa}{\kappa+2} & \frac{4}{\kappa+2} & -4 \\
 * \end{pmatrix}
 * \f]
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$\pi\f$:
 * \f[
 * Q = S . \pi = \begin{pmatrix}
 * -1 & \frac{1}{\kappa+2} & \frac{\kappa}{\kappa+2} & \frac{1}{\kappa+2} \\
 * \frac{1}{\kappa+2} & -1 & \frac{1}{\kappa+2} & \frac{\kappa}{\kappa+2} \\
 * \frac{\kappa}{\kappa+2} & \frac{1}{\kappa+2} & -1 & \frac{1}{\kappa+2} \\
 * \frac{1}{\kappa+2} & \frac{\kappa}{\kappa+2} & \frac{1}{\kappa+2} & -1 \\
 * \end{pmatrix}
 * \f]
 *
 * The eigen values are \f$\left(0, -\frac{2\kappa+2}{\kappa+2}, -\frac{2\kappa+2}{\kappa+2}, -\frac{4}{\kappa+2}\right)\f$,
 * the left eigen vectors are, by row:
 * \f[
 * U = \begin{pmatrix}
 * \frac{1}{4} &  \frac{1}{4} &  \frac{1}{4} &  \frac{1}{4} \\
 *           0 &  \frac{1}{2} &            0 & -\frac{1}{2} \\
 * \frac{1}{2} &            0 & -\frac{1}{2} &            0 \\
 * \frac{1}{4} & -\frac{1}{4} &  \frac{1}{4} & -\frac{1}{4} \\
 * \end{pmatrix}
 * \f]
 * and the right eigen vectors are, by column:
 * \f[
 * U^-1 = \begin{pmatrix}
 * 1 &  0 &  1 &  1 \\
 * 1 &  1 &  0 & -1 \\
 * 1 &  0 & -1 &  1 \\
 * 1 & -1 &  0 & -1 \\
 * \end{pmatrix}
 * \f]
 *
 * In addition, a rate_ factor defines the mean rate of the model.
 *
 * The probabilities of changes are computed analytically using the formulas:
 * \f[
 * P_{i,j}(t) = \begin{pmatrix}
 * \frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} & \frac{1}{4} - \frac{1}{4}B & -\frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} & \frac{1}{4} - \frac{1}{4}B \\
 * \frac{1}{4} - \frac{1}{4}B & \frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} & \frac{1}{4} - \frac{1}{4}B & -\frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} \\
 * -\frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} & \frac{1}{4} - \frac{1}{4}B & \frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} & \frac{1}{4} - \frac{1}{4}B \\
 * \frac{1}{4} - \frac{1}{4}B & -\frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} & \frac{1}{4} - \frac{1}{4}B & \frac{1}{2}A + \frac{1}{4}B + \frac{1}{4} \\
 * \end{pmatrix}
 * \f]
 * with \f$A=e^{-\frac{rate\_ * (2\kappa+2)t}{\kappa+2}}\f$ and \f$B = e^{-\frac{rate\_ * 4t}{\kappa+2}}\f$.
 *
 * First and second order derivatives are also computed analytically using the formulas:
 * \f[
 * \frac{\partial P_{i,j}(t)}{\partial t} = rate\_ * \begin{pmatrix}
 * -\frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B & \frac{1}{\kappa+2}B & \frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B & \frac{1}{\kappa+2}B \\
 * \frac{1}{\kappa+2}B & -\frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B & \frac{1}{\kappa+2}B & \frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B \\
 * \frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B & \frac{1}{\kappa+2}B & -\frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B & \frac{1}{\kappa+2}B \\
 * \frac{1}{\kappa+2}B & \frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B & \frac{1}{\kappa+2}B & -\frac{2\kappa+2}{2(\kappa+2)}A - \frac{1}{\kappa+2}B \\
 * \end{pmatrix}
 * \f]
 * \f{multline*}
 * \frac{\partial^2 P_{i,j}(t)}{\partial t^2} = rate\_^2 * \\
 * \begin{pmatrix}
 * \frac{{(2\kappa+2)}^2}{2{(\kappa+2)}^2}A - \frac{4}{{(\kappa+2)}^2}B & -\frac{4}{{(\kappa+2)}^2}B & -\frac{{(2\kappa+2)}^2}{2{(\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B & -\frac{4}{{(\kappa+2)}^2}B \\
 * -\frac{4}{{(\kappa+2)}^2}B & \frac{{(2\kappa+2)}^2}{2{(\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B & -\frac{4}{{(\kappa+2)}^2}B & -\frac{{(2\kappa+2)}^2}{2{(\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B \\
 * -\frac{{(2\kappa+2)}^2}{2{(\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B & -\frac{4}{{(\kappa+2)}^2}B & \frac{{(2\kappa+2)}^2}{2{(\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B & -\frac{4}{{(\kappa+2)}^2}B \\
 * -\frac{4}{{(\kappa+2)}^2}B & -\frac{2{(\kappa+2)}^2}{2{(\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B & -\frac{4}{{(\kappa+2)}^2}B & \frac{2{(\kappa+2)}^2}{{(2\kappa+2)}^2}A + \frac{4}{{(\kappa+2)}^2}B \\
 * \end{pmatrix}
 * \f}
 *
 * The parameter is named \c "kappa"
 * and its value may be retrieve with the command
 * \code
 * getParameterValue("kappa")
 * \endcode
 *
 * Reference:
 * - Kimura M (1980), Journal_ Of Molecular Evolution_ 16(2) 111-20.
 */
class K80 :
  public AbstractReversibleNucleotideSubstitutionModel
{
private:
  double kappa_, r_;
  mutable double l_, k_, exp1_, exp2_;
  mutable RowMatrix<double> p_;

public:
  K80(std::shared_ptr<const NucleicAlphabet> alpha, double kappa = 1.);

  virtual ~K80() {}

  K80* clone() const override { return new K80(*this); }

public:
  double Pij_t    (size_t i, size_t j, double d) const override;
  double dPij_dt  (size_t i, size_t j, double d) const override;
  double d2Pij_dt2(size_t i, size_t j, double d) const override;
  const Matrix<double>& getPij_t    (double d) const override;
  const Matrix<double>& getdPij_dt  (double d) const override;
  const Matrix<double>& getd2Pij_dt2(double d) const override;

  std::string getName() const override { return "K80"; }

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
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_K80_H
