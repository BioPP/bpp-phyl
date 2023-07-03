//
// File: TwoParameterBinarySubstitutionModel.h
// Authors:
//   Laurent Gueguen
// Created: 2009-10-08 00:00:00
//

/*
  Copyright or Ã¯Â¿Â½ or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

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
