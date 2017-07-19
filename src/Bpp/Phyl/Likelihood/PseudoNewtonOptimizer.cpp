//
// File: PseudoNewtonOptimizer.cpp
// Created by: Julien Dutheil
// Created on: Tue Nov 16 12:33 2004
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

/**************************************************************************/

#include "../Tree/TreeTemplateTools.h"
#include "PseudoNewtonOptimizer.h"
#include "DRHomogeneousTreeLikelihood.h"

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>

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
  
PseudoNewtonOptimizer::PseudoNewtonOptimizer(DerivableSecondOrder* function) :
  AbstractOptimizer(function),
  previousPoint_(),
  previousValue_(0),
  n_(0),
  params_(),
  maxCorrection_(10),
  useCG_(true)
{
  setDefaultStopCondition_(new FunctionStopCondition(this));
  setStopCondition(*getDefaultStopCondition());
}

/**************************************************************************/

void PseudoNewtonOptimizer::doInit(const ParameterList& params) throw (Exception)
{
  n_ = getParameters().size();
  params_ = getParameters().getParameterNames();
  getFunction()->enableSecondOrderDerivatives(true);
  getFunction()->setParameters(getParameters());
}

/**************************************************************************/

double PseudoNewtonOptimizer::doStep() throw (Exception)
{
  ParameterList* bckPoint = 0;
  if (updateParameters()) bckPoint = new ParameterList(getFunction()->getParameters());
  double newValue = 0;
  // Compute derivative at current point:
  std::vector<double> movements(n_);
  ParameterList newPoint = getParameters();

  for (size_t i = 0; i < n_; i++)
  {
    double firstOrderDerivative = getFunction()->getFirstOrderDerivative(params_[i]);
    double secondOrderDerivative = getFunction()->getSecondOrderDerivative(params_[i]);
    
    if (secondOrderDerivative == 0)
    {
      movements[i] = 0;
    }
    else if (secondOrderDerivative < 0)
    {
      printMessage("!!! Second order derivative is negative for parameter " + params_[i] + "(" + TextTools::toString(getParameters()[i].getValue()) + "). Moving in the other direction.");
      //movements[i] = 0;  // We want to reach a minimum, not a maximum!
      // My personnal improvement:
      movements[i] = -firstOrderDerivative / secondOrderDerivative;
    }
    else movements[i] = firstOrderDerivative / secondOrderDerivative;
    if (std::isnan(movements[i]))
    {
      printMessage("!!! Non derivable point at " + params_[i] + ". No move performed. (f=" + TextTools::toString(currentValue_) + ", d1=" + TextTools::toString(firstOrderDerivative) + ", d2=" + TextTools::toString(secondOrderDerivative) + ").");
      movements[i] = 0; // Either first or second order derivative is infinity. This may happen when the function == inf at this point.
    }
    //DEBUG:
    //cerr << "PN[" << params_[i] << "]=" << getParameters().getParameter(params_[i]).getValue() << "\t" << movements[i] << "\t " << firstOrderDerivative << "\t" << secondOrderDerivative << endl;
    newPoint[i].setValue(getParameters()[i].getValue() - movements[i]);
    //Correct the movement in case of constraint (this is used in case of Felsenstein-Churchill correction:
    movements[i] = getParameters()[i].getValue() - newPoint[i].getValue(); 
  }
  newValue = getFunction()->f(newPoint);

  // Check newValue:
  unsigned int count = 0;
  while ((count < maxCorrection_) && ((newValue > currentValue_ + getStopCondition()->getTolerance()) || std::isnan(newValue))) {
    //Restore previous point (all parameters in case of global constraint):
    if ((count==0) && updateParameters()) getFunction()->setParameters(*bckPoint);
    
    if (!(useCG_ && (count==3))){
      printMessage("!!! Function at new point is greater than at current point: " + TextTools::toString(newValue) + ">" + TextTools::toString(currentValue_) + ". Applying Felsenstein-Churchill correction: " + TextTools::toString(count));
        
      for (size_t i = 0; i < movements.size(); i++) {
        movements[i] = movements[i] / 2;
        newPoint[i].setValue(getParameters()[i].getValue() - movements[i]);
      }
      newValue = getFunction()->f(newPoint);
    }
    else {
      printMessage("!!! Felsenstein-Churchill correction applied too many times.");
      printMessage("Use conjugate gradients optimization.");
      getFunction()->enableSecondOrderDerivatives(false);
      ConjugateGradientMultiDimensions opt(getFunction());
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
      
      if (newValue > currentValue_ + tol) {
        printMessage("!!! Conjugate gradient method failed to improve likelihood.");
        printMessage("Back to Felsenstein-Churchill method.");
      }
    }
    count++;
  }
  
  if (newValue > currentValue_ + getStopCondition()->getTolerance()){
    printMessage("PseudoNewtonOptimizer::doStep. Value could not be ameliorated!");
    newValue = currentValue_;
  }
  else  {
    //getFunction()->enableFirstOrderDerivatives(true);
    getFunction()->enableSecondOrderDerivatives(true);
    getFunction()->setParameters(newPoint); //Compute derivatives for this point
    
    previousPoint_ = getParameters();
    previousValue_ = currentValue_;
    getParameters_() = newPoint;
  }

  if (updateParameters()) delete bckPoint;
  return newValue;
  
}

/**************************************************************************/

