//
// File: RN95.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: jeudi 24 fÃÂ©vrier 2011, ÃÂ  20h 43
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_MODEL_NUCLEOTIDE_RN95_H
#define BPP_PHYL_MODEL_NUCLEOTIDE_RN95_H

#include <Bpp/Numeric/Constraints.h>

#include "../AbstractSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief The model described by Rhetsky \& Nei, where the only
 * hypothesis is that the transversion rates are only dependent of
 * the target nucleotide. This model is not reversible.
 *
 * This model has been thoroughly studied by Schadt & al, and we
 * follow their notations.
 *
 * \f[
 * Q= \frac 1P  \begin{pmatrix}
 * . & \gamma & \alpha & \lambda \\
 * \delta & . &  \kappa & \beta \\
 * \epsilon & \gamma & . & \lambda \\
 * \delta & \sigma & \kappa & .\\
 * \end{pmatrix}
 * \f]
 * with
 *
 * \f[ P = \frac 1{\lambda+\gamma+\kappa+\delta}.\left(
 * \frac{\left(\kappa+\delta+\beta\right)\,\left(\sigma\,\lambda+\sigma
 * \,\gamma+\kappa\,\gamma+\delta\,\gamma\right)+\left(\sigma+\kappa+
 * \delta\right)\,\left(\kappa\,\lambda+\delta\,\lambda+\beta\,\lambda+
 * \beta\,\gamma\right)}{\sigma+\kappa+\delta+\beta}+\frac{\left(
 * \lambda+\gamma+\epsilon\right)\,\left(\kappa\,\lambda+\kappa\,
 * \gamma+\alpha\,\kappa+\alpha\,\delta\right)+\left(\lambda+\gamma+
 * \alpha\right)\,\left(\delta\,\lambda+\delta\,\gamma+\epsilon\,
 * \kappa+\delta\,\epsilon\right)}{\lambda+\gamma+\epsilon+ \alpha}
 * \right)\f]
 *
 * After normalization, this model has 7 free parameters so without
 * loss of generality we set: \f[\gamma+\lambda+\delta+\kappa=1\f]
 *
 * We use
 * \f[\begin{cases}
 * \gamma + \lambda + \delta + \kappa=1\\
 * c_1 = \kappa + \delta + \sigma + \beta\\
 * c_2 = \lambda + \gamma - (\sigma + \beta)  = 1 - c_1\\
 * c_3 = \lambda + \gamma + \alpha + \epsilon \\ 
 * c_4 = \kappa + \delta - (\alpha + \epsilon) = 1 - c_3 \\
 * \end{cases}
 *\f]
 
 * The stationnary distribution is then:
 * \f[
 * \pi = \frac{\delta\,(\lambda+\gamma)+\epsilon\,(\kappa+\delta)}{c_3} ,
 * \frac{\sigma\,(\lambda+\gamma)+\gamma\,(\kappa+\delta)}{c_1} ,
 * \frac{\kappa\,(\lambda+\gamma)+\alpha\,(\kappa+\delta)}{c_3} ,
 * \frac{\beta\,(\lambda+\gamma)+\lambda\,(\kappa+\delta)}{c_1}
 * \f]
 *
 * which means: \f[\pi_R = \delta + \kappa \f] and \f[\pi_Y = \gamma + \lambda \f].
 *
 * So we set as parameters with values in ]0;1[:
 *
 *\f[
 * \begin{cases}
 * \theta_R = \pi_A + \pi_G\\
 * \kappa'=\frac{\kappa}{\theta_R}\\
 * \gamma'=\frac{\gamma}{1-\theta_R}\\
 * \end{cases}
 * \f]
 *
 * from which:
 *
 * \f[
 * \begin{cases}
 * \kappa=\kappa' \theta_R\\
 * \delta=(1 - \kappa') \theta_R \\
 * \gamma=\gamma' (1-\theta_R)\\
 * \lambda=(1 -\gamma') (1-\theta_R)\\
 * \end{cases}
 * \f]
 *
 *
 * and 4 other positive parameters: \f[ \alpha, \beta, \epsilon, \sigma \f] 
 *
 * The eigen values are \f$\left(-\frac{1}{P}, - \frac{c_3}{P}, -\frac{c_1}{P}, 0\right)\f$,
 *
 * the right eigen vectors are, by column:
 * \f[
 * U^{-1} = \begin{pmatrix}
 * \lambda + \gamma &  (\kappa - \alpha) (\lambda + \gamma) + \alpha c_4 & \sigma \lambda - \beta \gamma  &  1 \\
 * -(\kappa + \delta) &  \alpha \delta - \epsilon \kappa &  (\kappa + \delta) (\beta - \lambda) - \beta  c_2 & 1 \\
 * \lambda + \gamma &  (\epsilon - \delta) (\lambda + \gamma) - \epsilon c_4  & \sigma \lambda - \beta \gamma &  1 \\
 * -(\kappa + \delta)& \alpha \delta - \epsilon \kappa & (\kappa + \delta) (\gamma - \sigma) + \sigma c_2 &1 \\
 * \end{pmatrix}
 * \f]
 *
 * 
 * and the left eigen vectors are, by row:
 * \f[
 * U = \begin{pmatrix}
 * \frac{\delta-\epsilon}{c_4} &  \frac{\sigma  - \gamma}{c_2} &  \frac{\kappa - \alpha}{c_4} &  \frac{\beta - \lambda}{c_2} \\
 *  \frac 1{c_4 c_3}    &   0 &  -\frac{1}{c_4 c_3} &  0 \\
 *  - \frac{1}{c_1 c_2} &   0 & \frac{1}{c_1 c_2}   &  0 \\
 * \frac{\epsilon + (\delta - \epsilon ) (\lambda +  \gamma)}{c_3} & \frac{\gamma + (\sigma - \gamma) (\lambda + \gamma) }{c_1} &  \frac{\alpha + (\kappa - \alpha) (\lambda +  \gamma)}{c_3} & \frac{\lambda + (\beta - \lambda)(\lambda + \gamma)}{c_1} \\
 * \end{pmatrix}
 * \f]
 *
 *
 *
 * The parameters are named \c "thetaR", \c "kappaP", \c "gammaP", \c
 * "alpha". \c "beta", \c "epsilon",\c "sigma".
 *
 * References:
 * - Rhetsky A. \& Nei M. (1995) MBE 12(1) 131-151.
 * - Schadt, Sinsheimer \& Lange (1998) Genome Research 8 222-233.
 *
 */

class RN95 :
  public AbstractNucleotideSubstitutionModel
{
private:
  double alpha_, beta_, gamma_, delta_, epsilon_, kappa_, lambda_, sigma_;

public:
  RN95(
    const NucleicAlphabet* alphabet,
    double alpha = 1,
    double beta = 1,
    double gamma = 0.25,
    double delta = 0.25,
    double epsilon = 1,
    double kappa = 0.25,
    double lambda = 0.25,
    double sigma = 1);

  virtual ~RN95() {}

  RN95* clone() const { return new RN95(*this); }

public:
  std::string getName() const { return "RN95"; }

  void updateMatrices();

  void setFreq(std::map<int, double>&);
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_RN95_H
