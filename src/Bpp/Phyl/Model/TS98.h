// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_TS98_H
#define BPP_PHYL_MODEL_TS98_H


#include "MarkovModulatedSubstitutionModel.h"

// From NumCalc:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Parameter.h>

namespace bpp
{
/**
 * @brief Tuffley and Steel's 1998 covarion model.
 *
 * This model is a subclass of the so-called Markov-modulated substitution models,
 * with a rate matrix
 * @f[
 * G = \begin{pmatrix}
 * -s_1 & s_1\\
 * s_2 & -s_2
 * \end{pmatrix}
 * @f]
 * and
 * @f[
 * D_R = \begin{pmatrix}
 * 0 & 0\\
 * 0 & \dfrac{s_1+s_2}{s_1}
 * \end{pmatrix}.
 * @f]
 * This model was originally designed for nucleotides sequences, but it can be used with other alphabets.
 *
 * @see MarkovModulatedSubstitutionModel
 *
 * Tuffley C. and Steel M. A., Modelling the covarion hypothesis of nucleotide substitution (1998),
 * _Math. Biosci._, 147:63-91.
 */
class TS98 :
  public MarkovModulatedSubstitutionModel
{
public:
  /**
   * @brief Build a new TS98 substitution model.
   *
   * @param model The substitution model to use. May be of any alphabet type.
   * @param s1    First rate parameter.
   * @param s2    Second rate parameter.
   * @param normalizeRateChanges Tell if the rate transition matrix should be normalized.
   */
  TS98(
      std::unique_ptr<ReversibleSubstitutionModelInterface> model,
      double s1 = 1.,
      double s2 = 1.,
      bool normalizeRateChanges = false) :
    MarkovModulatedSubstitutionModel(std::move(model), 2, normalizeRateChanges, "TS98.")
  {
    addParameter_(new Parameter("TS98.s1", s1, Parameter::R_PLUS_STAR));
    addParameter_(new Parameter("TS98.s2", s2, Parameter::R_PLUS_STAR));
    updateRatesModel_();
    updateMatrices_();
  }

  virtual ~TS98() {}

  TS98* clone() const override { return new TS98(*this); }

public:
  std::string getName() const override { return "TS98"; }

  double getRate() const override { return 1.; }

  void setRate(double rate) override {}

  void addRateParameter() override {}

protected:
  
  void updateRatesModel_() override
  {
    double s1 = getParameterValue("s1");
    double s2 = getParameterValue("s2");
    ratesFreq_[0] = s2 / (s1 + s2);
    ratesFreq_[1] = s1 / (s1 + s2);
    rates_(1, 1) = (s1 + s2) / s1;
    ratesExchangeability_(0, 1) = ratesExchangeability_(1, 0) = s1 + s2;
    ratesExchangeability_(0, 0) = -s1 * (s1 + s2) / s2;
    ratesExchangeability_(1, 1) = -s2 * (s1 + s2) / s1;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_TS98_H
