// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_NUCLEOTIDE_L95_H
#define BPP_PHYL_MODEL_NUCLEOTIDE_L95_H

#include <Bpp/Numeric/Constraints.h>

#include "../AbstractSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief The no-strand bias substitution model for nucleotides, from
 * Lobry 1995. The point of this model is that the substitution rate
 * from a nucleotide N towards another M is the same as the rate from
 * the complement of N towards the complement of M. Note that this
 * model is not reversible.
 *
 * After normalization, this model contains 5 parameters:
 * \f[
 * Q = \frac 1{2*\kappa*\theta*(1-\theta)+\gamma+\theta-2*\theta*\gamma} \begin{pmatrix}
 * \cdots & \kappa.\beta.\theta & \kappa.(1-\beta).\theta & \gamma \\
 * \kappa.\alpha.(1-\theta) & \cdots & 1-\gamma & \kappa.(1-\alpha).(1-\theta) \\
 * \kappa.(1-\alpha).(1-\theta) & 1-\gamma & \cdots & \kappa.\alpha.(1-\theta) \\
 * \gamma & \kappa.(1-\beta).\theta & \kappa.\beta.\theta & \cdots \\
 * \end{pmatrix}
 * \f]
 * The equilibrium frequencies are
 * \f[
 * \pi = \left(\frac{1-\theta}{2}, \frac{\theta}{2}, \frac{\theta}{2}, \frac{1-\theta}{2}\right)
 * \f]
 *
 * and then \f$\sum_i Q_{i,i}\pi_i = -1\f$.
 *
 * The generator of this model is diagonalized numerically.
 * See AbstractSubstitutionModel for details of how the probabilities are computed.
 *
 * The parameters are named \c "alpha", \c "beta", \c "gamma", \c
 * "kappa" and \c "theta". The values of \c "alpha", \c "beta", \c
 * "gamma" are between 0 and 1, The values of \c "gamma" are between 0
 * and 1 excluded, the values of "kappa" which are positive. Their
 * values may be retrieved with the command:
 *
 * \code
 * getParameterValue("alpha")
 * \endcode for instance.
 *
 * Reference:
 * - Lobry J R (1995), _Journal Of Molecular Evolution_ 40 326-330.
 */
class L95 :
  public AbstractNucleotideSubstitutionModel
{
private:
  double alpha_, beta_, gamma_, kappa_, theta_;

public:
  L95(
    std::shared_ptr<const NucleicAlphabet> alphabet,
    double alpha = 0.5,
    double beta = 0.5,
    double gamma = 0.5,
    double kappa = 1.,
    double theta = 0.5);

  virtual ~L95() {}

  L95* clone() const override { return new L95(*this); }

public:

  std::string getName() const override { return "L95"; }

  /**
   * @brief This method is redefined to actualize the corresponding parameters theta too.
   */
  void setFreq(std::map<int, double>&) override;
  
protected:

  void updateMatrices_() override;

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_L95_H
