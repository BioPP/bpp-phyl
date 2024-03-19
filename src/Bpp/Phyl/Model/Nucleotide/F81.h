// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_NUCLEOTIDE_F81_H
#define BPP_PHYL_MODEL_NUCLEOTIDE_F81_H

#include <Bpp/Numeric/Constraints.h>

#include "../AbstractSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief The Felsenstein  (1981) substitution model for nucleotides.
 *
 * This model is similar to J69 model,
 * but allows distinct equilibrium frequencies between A, C, G and T.
 * \f[
 * \begin{pmatrix}
 * \cdots & r & r & r \\
 * r & \cdots & r & r \\
 *  r & r & \cdots & r \\
 * r &  r & r & \cdots \\
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = \left(\pi_A, \pi_C, \pi_G, \pi_T\right)
 * \f]
 * This models hence includes four frequencies parameters \f$\pi_A, \pi_C, \pi_G, \pi_T\f$.
 * These four frequencies are not independent parameters, since they have the constraint to
 * sum to 1.
 * We use instead a different parametrization to remove this constraint:
 * \f[
 * \begin{cases}
 * \theta = \pi_C + \pi_G\\
 * \theta_1 = \frac{\pi_A}{1 - \theta} = \frac{\pi_A}{\pi_A + \pi_T}\\
 * \theta_2 = \frac{\pi_G}{\theta} = \frac{\pi_G}{\pi_C + \pi_G}\\
 * \end{cases}
 * \Longleftrightarrow
 * \begin{cases}
 * \pi_A = \theta_1 (1 - \theta)\\
 * \pi_C = (1 - \theta_2) \theta\\
 * \pi_G = \theta_2 \theta\\
 * \pi_T = (1 - \theta_1)(1 - \theta).
 * \end{cases}
 * \f]
 * These parameters can also be measured from the data and not optimized.
 *
 * Normalization: \f$r\f$ is set so that \f$\sum_i Q_{i,i}\pi_i = -1\f$:
 * \f[
 * S = \frac{1}{P}\begin{pmatrix}
 * \frac{-\pi_T-\pi_G-\pi_C}{\pi_A} & 1 & 1 & 1 \\
 * 1 & \frac{-\pi_T-\pi_G-\pi_A}{\pi_C} & 1 & 1 \\
 * 1 & 1 & \frac{-\pi_T-\pi_C-\pi_A}{\pi_G} & 1 \\
 * 1 & 1 & 1 & \frac{-\pi_G-\pi_C-\pi_A}{\pi_T} \\
 * \end{pmatrix}
 * \f]
 * with \f$P=2\left(\pi_A \pi_C + \pi_C \pi_G + \pi_A \pi_T + \pi_G \pi_T +  \left(\pi_C \pi_T + \pi_A \pi_G\right)\right)\f$.
 *
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$\pi\f$:
 * \f[
 * Q = S . \pi = \frac{1}{P}\begin{pmatrix}
 * -\pi_T-\pi_G-\pi_C & \pi_C & \pi_G & \pi_T \\
 * \pi_A & -\pi_T-\pi_G-\pi_A & \pi_G & \pi_T \\
 * \pi_A & \pi_C & -\pi_T-\pi_C-\pi_A & \pi_T \\
 * \pi_A & \pi_C & \pi_G & -\pi_G-\pi_C-\pi_A \\
 * \end{pmatrix}
 * \f]
 *
 * The eigen values are \f$\left(0, -\frac{\pi_R + \pi_Y}{P}, -\frac{\pi_Y + \pi_R}{P}, -\frac{1}{P}\right)\f$,
 * with \f$\pi_R=\pi_A+\pi_G\f$ and \f$\pi_Y=\pi_C+\pi_T\f$.
 * The left eigen vectors are, by row:
 * \f[
 * U = \begin{pmatrix}
 *                    \pi_A &               \pi_C &                    \pi_G & \pi_T \\
 *                        0 & \frac{\pi_T}{\pi_Y} &                        0 & -\frac{\pi_T}{\pi_Y} \\
 *      \frac{\pi_G}{\pi_R} &                   0 &     -\frac{\pi_G}{\pi_R} & 0 \\
 * \frac{\pi_A\pi_Y}{\pi_R} &              -\pi_C & \frac{\pi_G\pi_Y}{\pi_R} & -\pi_T \\
 * \end{pmatrix}
 * \f]
 * and the right eigen vectors are, by column:
 * \f[
 * U^-1 = \begin{pmatrix}
 * 1 &  0 &  1 &  1 \\
 * 1 &  1 &  0 & -\frac{\pi_R}{\pi_Y} \\
 * 1 &  0 & \frac{\pi_A}{\pi_G} &  1 \\
 * 1 & -\frac{\pi_C}{\pi_T} &  0 & -\frac{\pi_R}{\pi_Y} \\
 * \end{pmatrix}
 * \f]
 *
 * In addition, a rate_ factor defines the mean rate of the model.
 *
 * The probabilities of changes are computed analytically using the formulas:
 * \f{multline*}
 * P_{i,j}(t) = \\
 * \begin{pmatrix}
 * \frac{\pi_G}{\pi_R}A + \frac{\pi_A\pi_Y}{\pi_R}B + \pi_A & \pi_C - \pi_CB & -\frac{\pi_G}{\pi_R}A + \frac{\pi_G\pi_Y}{\pi_R}B + \pi_G & \pi_T - \pi_TB \\
 * \pi_A - \pi_AB & \frac{\pi_T}{\pi_Y}A + \frac{\pi_C\pi_R}{\pi_Y}B + \pi_C & \pi_G - \pi_GB & -\frac{\pi_T}{\pi_Y}A + \frac{\pi_T\pi_R}{\pi_Y}B + \pi_T \\
 * -\frac{\pi_A}{\pi_R}A + \frac{\pi_A\pi_Y}{\pi_R}B + \pi_A & \pi_C - \pi_CB & \frac{\pi_A}{\pi_R}A + \frac{\pi_G\pi_Y}{\pi_R}B + \pi_G & \pi_T - \pi_TB \\
 * \pi_A - \pi_AB & -\frac{\pi_C}{\pi_Y}A + \frac{\pi_C\pi_R}{\pi_Y}B + \pi_C & \pi_G - \pi_GB & \frac{\pi_C}{\pi_Y}A + \frac{\pi_R\pi_T}{\pi_Y}B + \pi_T \\
 * \end{pmatrix}
 * \f}
 * with \f$A=e^{-\frac{rate\_*(\pi_Y+\pi_R)t}{P}}\f$ and \f$B = e^{-\frac{rate\_*t}{P}}\f$.
 *
 * First and second order derivatives are also computed analytically using the formulas:
 * \f{multline*}
 * \frac{\partial P_{i,j}(t)}{\partial t} = rate\_ * \\
 * \frac{1}{P}
 * \begin{pmatrix}
 * -\frac{\pi_G(\pi_Y+\pi_R)}{\pi_R}A - \frac{\pi_A\pi_Y}{\pi_R}B & \pi_CB & \frac{\pi_G(\pi_Y+\pi_R)}{\pi_R}A - \frac{\pi_G\pi_Y}{\pi_R}B & \pi_TB \\
 * \pi_AB & -\frac{\pi_T(\pi_R+\pi_Y)}{\pi_Y}A - \frac{\pi_C\pi_R}{\pi_Y}B & \pi_GB & \frac{\pi_T(\pi_R+\pi_Y)}{\pi_Y}A - \frac{\pi_T\pi_R}{\pi_Y}B \\
 * \frac{\pi_A(\pi_Y+\pi_R)}{\pi_R}A - \frac{\pi_A\pi_Y}{\pi_R}B & \pi_CB & -\frac{\pi_A(\pi_Y+\pi_R)}{\pi_R}A - \frac{\pi_G\pi_Y}{\pi_R}B & \pi_TB \\
 * \pi_AB & \frac{\pi_C(\pi_R+\pi_Y)}{\pi_Y}A - \frac{\pi_C\pi_R}{\pi_Y}B & \pi_GB & -\frac{\pi_C(\pi_R+\pi_Y)}{\pi_Y}A - \frac{\pi_R\pi_T}{\pi_Y}B \\
 * \end{pmatrix}
 * \f}
 * \f{multline*}
 * \frac{\partial^2 P_{i,j}(t)}{\partial t^2} = rate\_^2 * \\
 * \frac{1}{P^2}
 * \begin{pmatrix}
 * \frac{\pi_G{(\pi_Y+\pi_R)}^2}{\pi_R}A + \frac{\pi_A\pi_Y}{\pi_R}B & -\pi_CB & -\frac{\pi_G{(\pi_Y+\pi_R)}^2}{\pi_R}A + \frac{\pi_G\pi_Y}{\pi_R}B & -\pi_TB \\
 * -\pi_AB & \frac{\pi_T{(\pi_R+\pi_Y)}^2}{\pi_Y}A + \frac{\pi_C\pi_R}{\pi_Y}B & -\pi_GB & -\frac{\pi_T{(\pi_R+\pi_Y)}^2}{\pi_Y}A + \frac{\pi_T\pi_R}{\pi_Y}B \\
 * -\frac{\pi_A{(\pi_Y+\pi_R)}^2}{\pi_R}A + \frac{\pi_A\pi_Y}{\pi_R}B & -\pi_CB & \frac{\pi_A{(\pi_Y+\pi_R)}^2}{\pi_R}A + \frac{\pi_G\pi_Y}{\pi_R}B & -\pi_TB \\
 * -\pi_AB & -\frac{\pi_C{(\pi_R+\pi_Y)}^2}{\pi_Y}A + \frac{\pi_C\pi_R}{\pi_Y}B & -\pi_GB & \frac{\pi_C{(\pi_R+\pi_Y)}^2}{\pi_Y}A + \frac{\pi_R\pi_T}{\pi_Y}B \\
 * \end{pmatrix}
 * \f}
 *
 * The parameters are named \c "theta", \c "theta1", and \c "theta2"
 * and their values may be retrieve with the command
 * \code
 * getParameterValue("theta")
 * \endcode
 * for instance.
 *
 * Reference:
 * - Felsenstein J. (1981), J Mol Evol. 1981;17(6):368-76.
 */

class F81 :
  public AbstractReversibleNucleotideSubstitutionModel
{
private:
  double r_, piA_, piC_, piG_, piT_, piY_, piR_, theta_, theta1_, theta2_;
  mutable double exp1_, l_;
  mutable RowMatrix<double> p_;

public:
  F81(
      std::shared_ptr<const NucleicAlphabet> alpha,
      double piA = 0.25,
      double piC = 0.25,
      double piG = 0.25,
      double piT = 0.25);

  virtual ~F81() {}

  F81* clone() const override { return new F81(*this); }

public:
  double Pij_t    (size_t i, size_t j, double d) const override;
  double dPij_dt  (size_t i, size_t j, double d) const override;
  double d2Pij_dt2(size_t i, size_t j, double d) const override;
  const Matrix<double>& getPij_t    (double d) const override;
  const Matrix<double>& getdPij_dt  (double d) const override;
  const Matrix<double>& getd2Pij_dt2(double d) const override;

  std::string getName() const override { return "F81"; }

  /**
   * @brief This method is redefined to actualize the corresponding parameters piA, piT, piG and piC too.
   */
  void setFreq(std::map<int, double>& freqs) override;

protected:
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_F81_H
