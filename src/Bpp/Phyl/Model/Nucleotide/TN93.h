// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_NUCLEOTIDE_TN93_H
#define BPP_PHYL_MODEL_NUCLEOTIDE_TN93_H

#include <Bpp/Numeric/Constraints.h>

#include "../AbstractSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief The Tamura and Nei (1993) substitution model for nucleotides.
 *
 * This model has two rate of transitions and one rate of transversion.
 * It also allows distinct equilibrium frequencies between A, C, G and T.
 * This models hence includes six parameters, two transition / transversion
 * relative rates \f$\kappa_1\f$ and \f$\kappa_2\f$, and four frequencies \f$\pi_A, \pi_C, \pi_G, \pi_T\f$.
 * These four frequencies are not independent parameters, since they have the constraint to
 * sum to 1. Usually, these parameters are measured from the data and not optimized.
 * \f[
 * S = \begin{pmatrix}
 * \cdots & r & \kappa_1 r & r \\
 * r & \cdots & r & \kappa_2 r \\
 * \kappa_1 r & r & \cdots & r \\
 * r & \kappa_2 r & r & \cdots \\
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = \left(\pi_A, \pi_C, \pi_G, \pi_T\right)
 * \f]
 * This models hence includes five parameters, two transition / transversion
 * relative rates \f$\kappa_1, \kappa_2\f$ and four frequencies \f$\pi_A, \pi_C, \pi_G, \pi_T\f$.
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
 * \frac{-\pi_T-\kappa_1\pi_G-\pi_C}{\pi_A} & 1 & \kappa_1 & 1 \\
 * 1 & \frac{-\kappa_2\pi_T-\pi_G-\pi_A}{\pi_C} & 1 & \kappa_2 \\
 * \kappa_1 & 1 & \frac{-\pi_T-\pi_C-\kappa_1\pi_A}{\pi_G} & 1 \\
 * 1 & \kappa_2 & 1 & \frac{-\pi_G-\kappa_2\pi_C-\pi_A}{\pi_T} \\
 * \end{pmatrix}
 * \f]
 * with \f$P=2\left(\pi_A \pi_C + \pi_C \pi_G + \pi_A \pi_T + \pi_G \pi_T + \kappa_2 \pi_C \pi_T + \kappa_1 \pi_A \pi_G\right)\f$.
 *
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$\pi\f$:
 * \f[
 * Q = S . \pi = \frac{1}{P}\begin{pmatrix}
 * -\pi_T-\kappa_1\pi_G-\pi_C & \pi_C & \kappa_1\pi_G & \pi_T \\
 * \pi_A & -\kappa_2\pi_T-\pi_G-\pi_A & \pi_G & \kappa_2\pi_T \\
 * \kappa_1\pi_A & \pi_C & -\pi_T-\pi_C-\kappa_1\pi_A & \pi_T \\
 * \pi_A & \kappa_2\pi_C & \pi_G & -\pi_G-\kappa_2\pi_C-\pi_A \\
 * \end{pmatrix}
 * \f]
 *
 * For now, the generator of this model is diagonalized numericaly.
 * See AbstractSubstitutionModel for details of how the porbabilities are computed.
 *
 * The parameters are named \c "kappa1", \c "kappa2", \c "theta", \c "theta1" and \c "theta2"
 * and their values may be retrieve with the command
 * \code
 * getParameterValue("kappa1")
 * \endcode
 * for instance.
 *
 * Reference:
 * - Tamura N and Nei K (1993), Molecular_ Biology And Evolution_ 10(3) 512-26.
 */
class TN93 :
  public AbstractReversibleNucleotideSubstitutionModel
{
private:
  double kappa1_, kappa2_, piA_, piC_, piG_, piT_, piY_, piR_, r_, k1_, k2_, theta_, theta1_, theta2_;
  mutable double exp1_, exp21_, exp22_, l_;
  mutable RowMatrix<double> p_;

public:
  TN93(
    std::shared_ptr<const NucleicAlphabet> alpha,
    double kappa1 = 1.,
    double kappa2 = 1.,
    double piA = 0.25,
    double piC = 0.25,
    double piG = 0.25,
    double piT = 0.25);

  virtual ~TN93() {}

  TN93* clone() const override { return new TN93(*this); }

public:
  double Pij_t    (size_t i, size_t j, double d) const override;
  double dPij_dt  (size_t i, size_t j, double d) const override;
  double d2Pij_dt2(size_t i, size_t j, double d) const override;
  const Matrix<double>& getPij_t    (double d) const override;
  const Matrix<double>& getdPij_dt  (double d) const override;
  const Matrix<double>& getd2Pij_dt2(double d) const override;

  std::string getName() const override { return "TN93"; }

  /**
   * @brief This method is over-defined to actualize the corresponding parameters piA, piT, piG and piC too.
   */
  void setFreq(std::map<int, double>& freqs) override;

protected:

  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_TN93_H
