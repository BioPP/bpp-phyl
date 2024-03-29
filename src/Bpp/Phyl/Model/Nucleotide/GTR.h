// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_NUCLEOTIDE_GTR_H
#define BPP_PHYL_MODEL_NUCLEOTIDE_GTR_H

#include <Bpp/Numeric/Constraints.h>

#include "../AbstractSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief The General Time-Reversible substitution model for nucleotides.
 *
 * This model is sometimes called the REV model, following Yang (1994), see references.
 * It was used in Lanave et al (1984), described in Tavare et al. 1886 and Rodriguez et al. 1990.
 * It is the most general reversible one, it has 6 substitution rates and 4 frequency
 * parameters.
 * We used the parametrization proposed by Yang (1994):
 * \f[
 * S = \begin{pmatrix}
 * \cdots & d & f & b \\
 * d & \cdots & e & a \\
 * f & e & \cdots & c \\
 * b & a & c & \cdots \\
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = \left(\pi_A, \pi_C, \pi_G, \pi_T\right)
 * \f]
 * This models hence includes nine parameters, five relative rates \f$a, b, c, d, e\f$
 * and four frequencies \f$\pi_A, \pi_C, \pi_G, \pi_T\f$.
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
 * Normalization: we set \f$f\f$ to 1, and scale the matrix so that \f$\sum_i Q_{i,i}\pi_i = -1\f$.
 * Parameters \f$a,b,c,d,e\f$ are hence relative rates.
 * \f[
 * S = \frac{1}{P}\begin{pmatrix}
 * \frac{-b\pi_T-\pi_G-d\pi_C}{\pi_A} & d & 1 & b \\
 * d & \frac{-a\pi_T-e\pi_G-d\pi_A}{\pi_C} & e & a \\
 * 1 & e & \frac{-c\pi_T-e\pi_C-\pi_A}{\pi_G} & c \\
 * b & a & c & \frac{-c\pi_G-a\pi_C-b\pi_A}{\pi_T} \\
 * \end{pmatrix}
 * \f]
 * with \f{eqnarray*}
 * P &=& \pi_A(d\pi_C+ \pi_G+b\pi_T)\\
 *   &+& \pi_C(d\pi_A+e\pi_G+a\pi_T)\\
 *   &+& \pi_G( \pi_A+e\pi_C+c\pi_T)\\
 *   &+& \pi_T(b\pi_A+a\pi_C+c\pi_G)\\
 *   &=& 2*(a*\pi_C*\pi_T+b*\pi_A*\pi_T+c*\pi_G*\pi_T+d*\pi_A*\pi_C+e*\pi_C*\pi_G+\pi_A*\pi_G)
 * \f}
 *
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$\pi\f$:
 * \f[
 * Q = S . \pi = \frac{1}{P}\begin{pmatrix}
 * -b\pi_T-\pi_G-d\pi_C & d\pi_C & \pi_G & b\pi_T \\
 * d\pi_A & -a\pi_T-e\pi_G-d\pi_A & e\pi_G & a\pi_T \\
 * \pi_A & e\pi_C & -c\pi_T-e\pi_C-\pi_A & c\pi_T \\
 * b\pi_A & a\pi_C & c\pi_G & -c\pi_G-a\pi_C-b\pi_A \\
 * \end{pmatrix}
 * \f]
 *
 * The generator of this model is diagonalized numericaly.
 * See AbstractSubstitutionModel for details of how the probabilities are computed.
 *
 * The parameters are named \c "a", \c "b", \c "c", \c "d", \c "e", \c "theta", \c "theta1", and \c "theta2"
 * and their values may be retrieve with the command
 * \code
 * getParameterValue("a")
 * \endcode
 * for instance.
 *
 * Reference:
 * - Yang Z (1994), Journal_ Of Molecular Evolution_ 39(1) 105-11.
 * - Lanave C, Preparata G, Saccone C and Serio G (1984), Journal_ Of Molecular Evolution_ 20 86-93.
 * - Tavaré S (1986), Lect_. Math. Life Sci._ 17 57-86.
 * - Rodriguez F (1990, Journal_ Of Theoretical Biology_ 142(4) 485-501.
 */

class GTR :
  public AbstractReversibleNucleotideSubstitutionModel
{
protected:
  double a_, b_, c_, d_, e_, piA_, piC_, piG_, piT_, theta_, theta1_, theta2_, p_;

public:
  GTR(
      std::shared_ptr<const NucleicAlphabet> alpha,
      double a = 1.,
      double b = 1.,
      double c = 1.,
      double d = 1.,
      double e = 1.,
      double piA = 0.25,
      double piC = 0.25,
      double piG = 0.25,
      double piT = 0.25);

  virtual ~GTR() {}

  GTR* clone() const override { return new GTR(*this); }

public:
  std::string getName() const override { return "GTR"; }

  /**
   * @brief This method is redefined to actualize the corresponding parameters piA, piT, piG and piC too.
   */
  void setFreq(std::map<int, double>& freqs) override;

protected:
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_GTR_H
