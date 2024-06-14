// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_NUCLEOTIDE_RN95S_H
#define BPP_PHYL_MODEL_NUCLEOTIDE_RN95S_H

#include <Bpp/Numeric/Constraints.h>

#include "../AbstractSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief Intersection of models RN95 and L95.
 *
 * The two hypotheses are that the transversion rates are only
 * dependent of the target nucleotide, and strand symmetry.
 *
 * After normalization this model has 3 parameters:
 * \f[
 * Q= \frac 1P
 * \begin{pmatrix}
 * . & \gamma & \alpha & \delta \\
 * \delta & . &  \gamma & \beta \\
 * \beta & \gamma & . & \delta \\
 * \delta & \alpha & \gamma & .\\
 * \end{pmatrix}\f]
 *
 * with
 *
 * \f[ P = \frac {
 * \left((\alpha+\gamma)(\delta+\beta+\gamma) +
 * (\delta+\beta)(\alpha+\gamma+\delta)\right)}{\alpha+\beta+\gamma+\delta}\f]
 *
 * After normalization, this model has 3 free parameters so without
 * loss of generality we set: \f[\alpha+\beta+\gamma+\delta=1\f].
 *
 * The stationnary distribution is:
 * \f[
 * \pi = \left(\frac{\delta + \beta}2, \frac{\alpha+ \gamma}2, \frac{\alpha+ \gamma}2, , \frac{\alpha+ \gamma}2\right)
 * \f]
 *
 * We use as  parameters:
 *
 *\f[
 * \begin{cases}
 * \theta =\pi_C + \pi_G = \alpha + \gamma\\
 * \alpha' =  \frac{\alpha}{\theta} \in ]0;1[\\
 * \beta' =  \frac{\beta}{1-\theta} \in ]0;1[\\
 * \end{cases}
 * \f]
 *
 * The generator is then computed as:
 *
 *\f[
 * \begin{cases}
 * \alpha=\theta \alpha'\\
 * \beta=(1-\theta) \beta'\\
 * \gamma = \theta - \alpha\\
 * \delta = 1 - \theta - \beta\\
 * \end{cases}
 * \f]
 *
 * Using
 * \f[
 * \begin{cases}
 * c_1 = \gamma+\delta-\beta-\alpha
 * c_2 = \delta \gamma - \alpha \beta
 * c_3 = \beta \gamma - \alpha \delta
 * \end{cases}
 * \f]
 *
 * The eigen values are \f$ \left(- \frac{2 (\gamma + \delta)}P, -\frac 1 P, -\frac 1 P, 0 \right)\f$,
 *
 * the right eigen vectors are, by column:
 * \f[
 * U^{-1} = \begin{pmatrix}
 * 1 &  c_2 & 0  &  1 \\
 * -1 & 0 &  c_2 &  1 \\
 * 1 &  (\beta - \delta) (\delta + \beta)  & c_3 &  1 \\
 * -1 & - c_3  & (\alpha - \gamma) (\alpha + \gamma) &1 \\
 * \end{pmatrix}
 * \f]
 *
 * and the left eigen vectors are, by row:
 * \f[
 * U = \begin{pmatrix}
 * \frac{\delta-\beta}{2 c_1} &  \frac{\alpha - \gamma}{2 c_1} &  \frac{\gamma  - \alpha}{2 c_1} & \frac{\beta-\delta}{2 c_1} \\
 * \frac{\gamma (\gamma + \delta) - \alpha (\alpha + \beta)}{c_1 c_2} &  \frac{c_3}{c_1 c_2} &  \frac{\alpha (\alpha + \beta) - \gamma (\gamma + \delta)}{c_1 c_2} & - \frac{c_3}{c_1 c_2} \\
 *  - \frac {c_3}{c_1 c_2}  & \frac{\delta(\gamma + \delta) - \beta(\alpha + \beta)}{c_1 c_2} &  \frac{c_3}{c_1 c_2} &  \frac{\beta(\alpha + \beta) - \delta(\gamma + \delta)}{c_1 c_2} \\
 * \frac{\delta+\beta}{2} &  \frac{\alpha + \gamma}{2} &  \frac{\gamma  + \alpha}{2} & \frac{\beta+\delta}{2} \\
 * \end{pmatrix}
 * \f]
 *
 *
 * The parameters are named \c "theta", \c "alphaP", \c "betaP".
 *
 * References:
 * - Rhetsky A. \& Ney M. (1995) MBE 12(1) 131-151.
 * - Lobry J R (1995), Journal_ Of Molecular Evolution_ 40 326-330.
 * - Schadt, Sinsheimer \& Lange (1998) Genome Research 8 222-233.
 */

class RN95s :
  public AbstractNucleotideSubstitutionModel
{
private:
  double alpha_, beta_, gamma_, delta_;

public:
  RN95s(std::shared_ptr<const NucleicAlphabet> alphabet,
      double alpha = 0.25,
      double beta = 0.25,
      double gamma = 0.25,
      double delta = 0.25);

  virtual ~RN95s() {}

  RN95s* clone() const override { return new RN95s(*this); }

public:
  std::string getName() const override { return "RN95s"; }

  /**
   * @brief This method takes the average value between observed @f$\pi_A@f$ and @f$\pi_T@f$.
   */
  void setFreq(std::map<int, double>&) override;

protected:
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_RN95S_H
