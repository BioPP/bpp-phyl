// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_NUCLEOTIDE_SSR_H
#define BPP_PHYL_MODEL_NUCLEOTIDE_SSR_H

#include <Bpp/Numeric/Constraints.h>

#include "../AbstractSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief The Strand Symmetric Reversible substitution model for
 * nucleotides.
 *
 * We use a parametrization derived from Hobolth et al 2007
 * \f[
 * S = \begin{pmatrix}
 * \cdots & \beta & 1 & \gamma \\
 * \beta & \cdots & \delta & 1 \\
 * 1 & \delta & \cdots & \beta \\
 * \gamma & 1 & \beta & \cdots \\
 * \end{pmatrix}
 * \f]
 * The equilibrium frequencies
 * \f[
 * \pi = \left(1-\frac{\theta}{2}, \frac{\theta}{2}, \frac{\theta}{2}, 1-\frac{\theta}{2}\right)
 * \f]
 * This models hence includes four parameters, three relative rates \f$\beta, \gamma, \delta\f$ and the GC content \f$\theta\f$.
 *
 * Normalization: we set \f$f\f$ to 1, and scale the matrix so that \f$\sum_i Q_{i,i}\pi_i = -1\f$.
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$\pi\f$:
 * \f[
 * Q = S . \pi = \frac{1}{P}\begin{pmatrix}
 * -\gamma\pi_T-\pi_G-\beta\pi_C & \beta\pi_C & \pi_G & \gamma\pi_T \\
 * \beta\pi_A & -\pi_T-\delta\pi_G-\beta\pi_A & \delta\pi_G & \pi_T \\
 * \pi_A & \delta\pi_C & -\beta\pi_T-\delta\pi_C-\pi_A & \beta\pi_T \\
 * \gamma\pi_A & \pi_C & \beta\pi_G & -\beta\pi_G-\pi_C-\gamma\pi_A \\
 * \end{pmatrix}
 * \f]
 * where P is the normalization constant.
 * For now, the generator of this model is diagonalized numericaly.
 * See AbstractSubstitutionModel for details of how the probabilities are computed.
 *
 * The parameters are named \c "beta", \c "gamma", \c "delta", and \c "theta"
 * and their values may be retrieved with the command
 * \code
 * getParameterValue("beta")
 * \endcode
 * for instance.
 *
 * References:
 * - Hobolth A, Christensen O Fm Mailund T, Schierup M H (2007), PLoS Genetics 3(2) e7.
 * - Yap VB, Speed TP (2004), Journal Of Molecular Evolution 58(1) 12-18
 */
class SSR :
  public AbstractReversibleNucleotideSubstitutionModel
{
private:
  double beta_, gamma_, delta_, theta_, piA_, piC_, piG_, piT_;

public:
  SSR(std::shared_ptr<const NucleicAlphabet> alpha,
      double beta = 1.,
      double gamma = 1.,
      double delta = 1.,
      double theta = 0.5);

  virtual ~SSR() {}

  SSR* clone() const override { return new SSR(*this); }

public:
  std::string getName() const override { return "SSR"; }

  /**
   * @brief This method is redefined to actualize the corresponding parameters theta too.
   */
  void setFreq(std::map<int, double>&) override;

protected:
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_SSR_H
