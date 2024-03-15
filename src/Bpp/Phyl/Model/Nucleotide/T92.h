// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_NUCLEOTIDE_T92_H
#define BPP_PHYL_MODEL_NUCLEOTIDE_T92_H

#include <Bpp/Numeric/Constraints.h>

#include "../AbstractSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief The Tamura (1992) substitution model for nucleotides.
 *
 * This model is similar to the K80 model,
 * but allows distinct equilibrium frequencies between GC and AT.
 * This models hence includes two parameters, the transition / transversion
 * relative rate \f$\kappa\f$ and the frequency of GC, \f$\theta\f$.
 * \f[
 * S = \begin{pmatrix}
 * \cdots & r & \kappa r & r \\
 * r & \cdots & r & \kappa r \\
 * \kappa r & r & \cdots & r \\
 * r & \kappa r & r & \cdots \\
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = \left(\frac{1-\theta}{2}, \frac{\theta}{2}, \frac{\theta}{2}, \frac{1 - \theta}{2}\right)
 * \f]
 * Normalization: \f$r\f$ is set so that \f$\sum_i Q_{i,i}\pi_i = -1\f$:
 * \f[
 * S = \frac{1}{P}\begin{pmatrix}
 * \frac{2\theta\kappa + 2}{\theta - 1} & 2 & 2\kappa & 2 \\
 * 2 & \frac{2(\theta-1)\kappa - 2}{\theta} & 2 & 2\kappa \\
 * 2\kappa & 2 & \frac{2(\theta-1)\kappa - 2}{\theta} & 2 \\
 * 2 & 2\kappa & 2 & \frac{2\theta\kappa + 2}{\theta - 1} \\
 * \end{pmatrix}
 * \f]
 * with \f$P=1+2\theta\kappa-2\theta^2\kappa\f$.
 *
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$\pi\f$:
 * \f[
 * Q = S . \pi = \frac{1}{P}\begin{pmatrix}
 * -(\theta\kappa+1) & \theta & \theta\kappa & 1-\theta \\
 * 1-\theta & -((1-\theta)\kappa+1) & \theta & (1-\theta)\kappa \\
 * (1-\theta)\kappa & \theta & -((1-\theta)\kappa+1) & 1-\theta \\
 * 1-\theta & \theta\kappa & \theta & -(\theta\kappa+1) \\
 * \end{pmatrix}
 * \f]
 *
 * The eigen values are \f$\left(0, -\frac{\kappa+1}{P}, -\frac{\kappa+1}{P}, -\frac{2}{P}\right)\f$,
 * the left eigen vectors are, by row:
 * \f[
 * U = \begin{pmatrix}
 * \frac{1-\theta}{2} &  \frac{\theta}{2} &  \frac{\theta}{2} & \frac{1-\theta}{2} \\
 *                  0 &          1-\theta &                 0 & \theta-1 \\
 *             \theta &                   0 &         -\theta & 0 \\
 * \frac{1-\theta}{2} & -\frac{\theta}{2} &  \frac{\theta}{2} & -\frac{1-\theta}{2} \\
 * \end{pmatrix}
 * \f]
 * and the right eigen vectors are, by column:
 * \f[
 * U^-1 = \begin{pmatrix}
 * 1 &  0 &  1 &  1 \\
 * 1 &  1 &  0 & -1 \\
 * 1 &  0 & -\frac{1-\theta}{\theta} &  1 \\
 * 1 & -\frac{\theta}{1-\theta} &  0 & -1 \\
 * \end{pmatrix}
 * \f]
 *
 * In addition, a rate_ factor defines the mean rate of the model.
 *
 * The probabilities of changes are computed analytically using the formulas:
 * \f{multline*}
 * P_{i,j}(t) = \\
 * \begin{pmatrix}
 * \theta A + \frac{1-\theta}{2}B + \frac{1-\theta}{2} & \frac{\theta}{2} - \frac{\theta}{2}B & -\theta A + \frac{\theta}{2}B + \frac{\theta}{2} & \frac{1-\theta}{2} - \frac{1-\theta}{2}B \\
 * \frac{1-\theta}{2} - \frac{1-\theta}{2}B & (1-\theta)A + \frac{\theta}{2}B + \frac{\theta}{2} & \frac{\theta}{2} - \frac{\theta}{2}B & -(1-\theta)A + \frac{1-\theta}{2}B + \frac{1-\theta}{2} \\
 * -(1-\theta)A + \frac{1-\theta}{2}B + \frac{1-\theta}{2} & \frac{\theta}{2} - \frac{\theta}{2}B & (1-\theta)A + \frac{\theta}{2}B + \frac{\theta}{2} & \frac{1-\theta}{2} - \frac{1-\theta}{2}B \\
 * \frac{1-\theta}{2} - \frac{1-\theta}{2}B & -\theta A + \frac{\theta}{2}B + \frac{\theta}{2} & \frac{\theta}{2} - \frac{\theta}{2}B & \theta A + \frac{1-\theta}{2}B + \frac{1-\theta}{2} \\
 * \end{pmatrix}
 * \f}
 * with \f$A=e^{-\frac{rate\_ * (\kappa+1)t}{P}}\f$ and \f$B = e^{-\frac{rate\_ * 2t}{P}}\f$.
 *
 * First and second order derivatives are also computed analytically using the formulas:
 * \f{multline*}
 * \frac{\partial P_{i,j}(t)}{\partial t} = rate\_ * \\
 * \frac{1}{P}
 * \begin{pmatrix}
 * -\theta(\kappa+1)A - (1-\theta)B & \theta B & \theta(\kappa+1)A - \theta B & (1-\theta)B \\
 * (1-\theta)B & -(1-\theta)(\kappa+1)A - \theta B & \theta B & (1-\theta)(\kappa+1)A - (1-\theta)B \\
 * (1-\theta)(\kappa+1)A - (1-\theta)B & \theta B & -(1-\theta)(\kappa+1)A - \theta B & (1-\theta)B \\
 * (1-\theta)B & \theta(\kappa+1)A - \theta B & \theta B & -\theta(\kappa+1)A - (1-\theta)B \\
 * \end{pmatrix}
 * \f}
 * \f{multline*}
 * \frac{\partial^2 P_{i,j}(t)}{\partial t^2} = rate\_^2 * \\
 * \frac{1}{P^2}
 * \begin{pmatrix}
 * \theta{(\kappa+1)}^2A + 2(1-\theta)B & -2\theta B & -\theta{(\kappa+1)}^2A + 2\theta B & -2(1-\theta)B \\
 * -2(1-\theta)B & (1-\theta){(\kappa+1)}^2A + 2\theta B & -2\theta B & -(1-\theta){(\kappa+1)}^2A + 2(1-\theta)B \\
 * -(1-\theta){(\kappa+1)}^2A + 2(1-\theta)B & -2\theta B & (1-\theta){(\kappa+1)}^2A + 2\theta B & -2(1-\theta)B \\
 * -2(1-\theta)B & -\theta{(\kappa+1)}^2A + 2\theta B & -2\theta B & \theta{(\kappa+1)}^2A + 2(1-\theta)B \\
 * \end{pmatrix}
 * \f}
 *
 * The parameters are named \c "kappa" and \c "theta"
 * and their values may be retrieve with the commands
 * \code
 * getParameterValue("kappa")
 * getParameterValue("theta")
 * \endcode
 *
 * Reference:
 * - Tamura K (1992), Molecular_ Biology And Evolution_ 9(5) 814-25.
 */
class T92 :
  public AbstractReversibleNucleotideSubstitutionModel
{
private:
	
  double kappa_, theta_, k_, r_, piA_, piC_, piG_, piT_;
  mutable double exp1_, exp2_, l_;
  mutable RowMatrix<double> p_;

public:

  T92(std::shared_ptr<const NucleicAlphabet> alpha, double kappa = 1., double theta = 0.5);

  virtual ~T92() {}

  T92* clone() const override { return new T92(*this); }

public:
  double Pij_t    (size_t i, size_t j, double d) const override;
  double dPij_dt  (size_t i, size_t j, double d) const override;
  double d2Pij_dt2(size_t i, size_t j, double d) const override;
  const Matrix<double>& getPij_t(double d) const override;
  const Matrix<double>& getdPij_dt(double d) const override;
  const Matrix<double>& getd2Pij_dt2(double d) const override;

  std::string getName() const override { return "T92"; }

  /**
   * @brief This method is over-defined to actualize the 'theta' parameter too.
   */
  void setFreq(std::map<int, double>& freqs) override;

protected:

  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_T92_H
