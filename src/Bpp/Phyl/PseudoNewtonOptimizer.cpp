// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

/**************************************************************************/

#include "PseudoNewtonOptimizer.h"

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>

#include "Likelihood/PhyloLikelihoods/PhyloLikelihood.h"
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>

using namespace bpp;

/**************************************************************************/

double PseudoNewtonOptimizer::PNStopCondition::getCurrentTolerance() const
{
  return NumTools::abs<double>(
        dynamic_cast<const PseudoNewtonOptimizer*>(optimizer_)->currentValue_ -
        dynamic_cast<const PseudoNewtonOptimizer*>(optimizer_)->previousValue_);
}

/**************************************************************************/

bool PseudoNewtonOptimizer::PNStopCondition::isToleranceReached() const
{
  return getCurrentTolerance() < tolerance_;
}

/**************************************************************************/

PseudoNewtonOptimizer::PseudoNewtonOptimizer(shared_ptr<SecondOrderDerivable> function) :
  AbstractOptimizer(function),
  previousPoint_(),
  previousValue_(0),
  n_(0),
  params_(),
  maxCorrection_(10),
  useCG_(true)
{
  setDefaultStopCondition_(make_shared<FunctionStopCondition>(this));
  setStopCondition(getDefaultStopCondition());
}

/**************************************************************************/

void PseudoNewtonOptimizer::doInit(const ParameterList& params)
{
  n_ = getParameters().size();
  params_ = getParameters().getParameterNames();
  secondOrderDerivableFunction().enableSecondOrderDerivatives(true);
  secondOrderDerivableFunction().setParameters(getParameters());
}

/**************************************************************************/

double PseudoNewtonOptimizer::doStep()
{
  ParameterList* bckPoint = 0;
  if (updateParameters())
    bckPoint = new ParameterList(getFunction()->getParameters());
  double newValue = 0;
  // Compute derivative at current point:
  std::vector<double> movements(n_);
  ParameterList newPoint = getParameters();

  for (size_t i = 0; i < n_; i++)
  {
    double firstOrderDerivative = secondOrderDerivableFunction().getFirstOrderDerivative(params_[i]);
    double secondOrderDerivative = secondOrderDerivableFunction().getSecondOrderDerivative(params_[i]);
    if (secondOrderDerivative == 0)
    {
      movements[i] = 0;
    }
    else if (secondOrderDerivative < 0)
    {
      printMessage("!!! Second order derivative is negative for parameter " + params_[i] + "(" + TextTools::toString(getParameters()[i].getValue()) + "). Moving in the other direction.");
      // movements[i] = 0;  // We want to reach a minimum, not a maximum!
      // My personal improvement:
      movements[i] = -firstOrderDerivative / secondOrderDerivative;
    }
    else
      movements[i] = firstOrderDerivative / secondOrderDerivative;
    if (std::isnan(movements[i]))
    {
      printMessage("!!! Non derivable point at " + params_[i] + ". No move performed. (f=" + TextTools::toString(currentValue_) + ", d1=" + TextTools::toString(firstOrderDerivative) + ", d2=" + TextTools::toString(secondOrderDerivative) + ").");
      movements[i] = 0; // Either first or second order derivative is infinity. This may happen when the function == inf at this point.
    }
    // DEBUG:
    // cerr << "PN[" << params_[i] << "]=" << getParameters().parameter(params_[i]).getValue() << "\t" << movements[i] << "\t " << firstOrderDerivative << "\t" << secondOrderDerivative << endl;
    newPoint[i].setValue(getParameters()[i].getValue() - movements[i]);
    // Correct the movement in case of constraint (this is used in case of Felsenstein-Churchill correction:
    movements[i] = getParameters()[i].getValue() - newPoint[i].getValue();
  }
  newValue = getFunction()->f(newPoint);

  // Check newValue:
  unsigned int count = 0;
  while ((count < maxCorrection_) && ((newValue > currentValue_ + getStopCondition()->getTolerance()) || std::isnan(newValue)))
  {
    // Restore previous point (all parameters in case of global constraint):
    if ((count == 0) && updateParameters())
      getFunction()->setParameters(*bckPoint);

    if (!(useCG_ && (count == 3)))
    {
      printMessage("!!! Function at new point is greater than at current point: " + TextTools::toString(newValue) + ">" + TextTools::toString(currentValue_) + ". Applying Felsenstein-Churchill correction: " + TextTools::toString(count));

      for (size_t i = 0; i < movements.size(); i++)
      {
        movements[i] = movements[i] / 2;
        newPoint[i].setValue(getParameters()[i].getValue() - movements[i]);
      }
      newValue = getFunction()->f(newPoint);
    }
    else
    {
      printMessage("!!! Felsenstein-Churchill correction applied too many times.");
      printMessage("Use conjugate gradients optimization.");
      secondOrderDerivableFunction().enableSecondOrderDerivatives(false);
      ConjugateGradientMultiDimensions opt(dynamic_pointer_cast<FirstOrderDerivable>(function_));
      opt.setConstraintPolicy(getConstraintPolicy());
      opt.setProfiler(getProfiler());
      opt.setMessageHandler(getMessageHandler());
      opt.setVerbose(0);
      double tol = std::max(getStopCondition()->getCurrentTolerance() / 2., getStopCondition()->getTolerance());
      opt.getStopCondition()->setTolerance(tol);
      opt.setMaximumNumberOfEvaluations(nbEvalMax_);
      getFunction()->setParameters(getParameters());
      opt.init(getParameters());
      opt.optimize();
      newPoint = opt.getParameters();
      newValue = opt.getFunctionValue();

      if (newValue > currentValue_ + tol)
      {
        printMessage("!!! Conjugate gradient method failed to improve likelihood.");
        printMessage("Back to Felsenstein-Churchill method.");
      }
    }
    count++;
  }

  if (newValue > currentValue_ + getStopCondition()->getTolerance())
  {
    printMessage("PseudoNewtonOptimizer::doStep. Value could not be ameliorated!");
    newValue = currentValue_;
  }
  else
  {
    secondOrderDerivableFunction().enableSecondOrderDerivatives(true);
    secondOrderDerivableFunction().setParameters(newPoint); // Compute derivatives for this point

    previousPoint_ = getParameters();
    previousValue_ = currentValue_;
    getParameters_() = newPoint;
  }

  if (updateParameters())
    delete bckPoint;
  return newValue;
}

/**************************************************************************/
