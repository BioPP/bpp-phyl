// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_TWOPARAMETERBINARYSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_TWOPARAMETERBINARYSUBSTITUTIONMODEL_H

#include <Bpp/Seq/Alphabet/BinaryAlphabet.h>

#include "AbstractSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The Model on two states
 *
 * \f[
 * Q = r.\begin{pmatrix}
 * -\kappa & \kappa  \\
 * 1 & -1  \\
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = diag\left(\frac{1}{\kappa+1}, \frac{\kappa}{\kappa+1}\right)
 * \f]
 * Normalization: \f$r\f$ is set so that \f$\sum_i Q_{i,i}\pi_i = -1\f$:
 * \f[
 * Q = \begin{pmatrix}
 * -\frac{\kappa + 1}2 & \frac{\kappa + 1}2 \\
 * \frac{\kappa+1}{2\kappa} & -\frac{\kappa+1}{2\kappa}\\
 * \end{pmatrix}
 * \f]
 *
 * The eigen values are \f$\left(0, - \frac{(\kappa+1)^2}{2\kappa}\right)\f$,
 * and IF \f$\kappa \neq 1\f$, the left eigen vectors are, by row:
 * \f[
 * U = \begin{pmatrix}
 *  \frac{1}{1+\kappa} &  \frac{\kappa}{1+\kappa} \\
 *  \frac{\kappa-1}{\kappa+1} & -\frac{\kappa-1}{\kappa+1} \\
 * \end{pmatrix}
 * \f]
 * and the right eigen vectors are by column:
 * \f[
 * U^{-1} = \begin{pmatrix}
 *  1 &  \frac \kappa{\kappa-1} \\
 *  1 &  - \frac 1{\kappa-1} \\
 * \end{pmatrix}
 * \f]
 *
 * The probabilities of changes are computed analytically using the formulas, with \f$\lambda= \frac{(\kappa+1)^2}{2\kappa}\f$ :
 * \f[
 * P_{i,j}(t) = \begin{pmatrix}
 * \frac{1}{\kappa+1} + \frac{\kappa}{\kappa+1}e^{-\lambda t} & \frac{\kappa}{\kappa+1} - \frac{\kappa}{\kappa+1}e^{-\lambda t} \\
 * \frac{1}{\kappa+1} - \frac{1}{\kappa+1}e^{-\lambda t} & \frac{\kappa}{\kappa+1} + \frac{1}{\kappa+1}e^{-\lambda t} \\
 * \end{pmatrix}
 * \f]
 *
 * \f[
 * \frac{\partial P_{i,j}(t)}{\partial t} = \begin{pmatrix}
 * -\frac {\kappa+1} 2 e^{-\lambda t}  & \frac {\kappa+1} 2 e^{-\lambda t} \\
 * \frac {\kappa+1} {2\kappa} e^{-\lambda t}  & - \frac {\kappa+1} {2\kappa} e^{-\lambda t} \\
 * \end{pmatrix}
 * \f]
 * \f{multline*}
 * \frac{\partial^2 P_{i,j}(t)}{\partial t^2} = \\
 * \begin{pmatrix}
 * \frac {\lambda (\kappa+1)} 2 e^{-\lambda t}  & -\ frac {\lambda (\kappa+1)} 2 e^{-\lambda t} \\
 * \frac {\lambda (\kappa+1)} {2\kappa} e^{-\lambda t}  & - \frac {\lambda (\kappa+1)} {2\kappa} e^{-\lambda t} \\
 * \end{pmatrix}
 * \f}
 *
 * The parameter is named \c "kappa"
 * and its value may be retrieve with the command
 * \code
 * getParameterValue("kappa")
 * \endcode
 *
 */

class TwoParameterBinarySubstitutionModel :
  public AbstractReversibleSubstitutionModel
{
private:
  double mu_;
  double pi0_;

protected:
  mutable double lambda_, exp_;
  mutable RowMatrix<double> p_;

public:
  TwoParameterBinarySubstitutionModel(std::shared_ptr<const BinaryAlphabet> alpha, double mu = 1., double pi0 = 0.5);

  virtual ~TwoParameterBinarySubstitutionModel() {}

  TwoParameterBinarySubstitutionModel* clone() const { return new TwoParameterBinarySubstitutionModel(*this); }

public:
  // the inherited functions don't do the work - need to override them with the correct computation
  double Pij_t    (size_t i, size_t j, double d) const;
  double dPij_dt  (size_t i, size_t j, double d) const;
  double d2Pij_dt2(size_t i, size_t j, double d) const;
  const Matrix<double>& getPij_t    (double d) const;
  const Matrix<double>& getdPij_dt  (double d) const;
  const Matrix<double>& getd2Pij_dt2(double d) const;

  std::string getName() const { return "TwoParameterBinary"; }

  size_t getNumberOfStates() const { return 2; }

  void setMuBounds(double lb, double ub);

protected:
  void updateMatrices_();
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_TWOPARAMETERBINARYSUBSTITUTIONMODEL_H
