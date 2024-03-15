// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_PSEUDONEWTONOPTIMIZER_H
#define BPP_PHYL_PSEUDONEWTONOPTIMIZER_H

#include <Bpp/Numeric/Function/AbstractOptimizer.h>


namespace bpp
{
/**
 * @brief This Optimizer implements Newton's algorithm for finding a minimum of a function.
 * This is in fact a modified version of the algorithm, as suggested by Nicolas Galtier, for
 * the purpose of optimizing phylogenetic likelihoods.
 *
 * Only second simple order derivative are computed, no cross derivative, following Galtier's
 * algorithm.
 * Felsenstein and Churchill's (1996) correction is applied when new trial as a likelihood
 * lower than the starting point.
 */
class PseudoNewtonOptimizer :
  public AbstractOptimizer
{
public:
  class PNStopCondition :
    public AbstractOptimizationStopCondition
  {
public:
    PNStopCondition(PseudoNewtonOptimizer* pno) :
      AbstractOptimizationStopCondition(pno) {}
    virtual ~PNStopCondition() {}

    PNStopCondition* clone() const { return new PNStopCondition(*this); }

public:
    bool isToleranceReached() const;
    double getCurrentTolerance() const;
  };

  friend class PNStopCondition;

private:
  ParameterList previousPoint_; // Current point is in parameters_

  double previousValue_;

  size_t n_; // Number of parameters

  std::vector<std::string> params_; // All parameter names

  unsigned int maxCorrection_;

  bool useCG_;

public:
  PseudoNewtonOptimizer(std::shared_ptr<SecondOrderDerivable> function);

  virtual ~PseudoNewtonOptimizer() {}

  PseudoNewtonOptimizer* clone() const { return new PseudoNewtonOptimizer(*this); }

public:
  const SecondOrderDerivable& secondOrderDerivableFunction() const
  {
    return *dynamic_pointer_cast<const SecondOrderDerivable>(function_);
  }

  SecondOrderDerivable& secondOrderDerivableFunction()
  {
    return *dynamic_pointer_cast<SecondOrderDerivable>(function_);
  }

  std::shared_ptr<const SecondOrderDerivable> getSecondOrderDerivableFunction() const
  {
    return dynamic_pointer_cast<const SecondOrderDerivable>(function_);
  }

  std::shared_ptr<SecondOrderDerivable> getSecondOrderDerivableFunction()
  {
    return dynamic_pointer_cast<SecondOrderDerivable>(function_);
  }

  /**
   * @name The Optimizer interface.
   *
   * @{
   */
  double getFunctionValue() const { return currentValue_; }
  /** @} */

  void doInit(const ParameterList& params);

  double doStep();

  void setMaximumNumberOfCorrections(unsigned int mx) { maxCorrection_ = mx; }

  void disableCG() { useCG_ = false; }

};
} // end of namespace bpp.
#endif // BPP_PHYL_PSEUDONEWTONOPTIMIZER_H
