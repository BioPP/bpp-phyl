//
// File: OptimizationTools.cpp
// Authors:
//   Julien Dutheil
// Created: 2003-12-14 09:43:32
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

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/TwoPointsNumericalDerivative.h>
#include <Bpp/Numeric/ParameterList.h>

#include "../Io/Newick.h"
#include "../PseudoNewtonOptimizer.h"
#include "../Tree/NNISearchable.h"
#include "Likelihood/GlobalClockTreeLikelihoodFunctionWrapper.h"
#include "OptimizationTools.h"
#include "Tree/NNITopologySearch.h"

// From bpp-seq:
#include <Bpp/Seq/Io/Fasta.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

LegacyOptimizationTools::LegacyOptimizationTools() {}

LegacyOptimizationTools::~LegacyOptimizationTools() {}

/******************************************************************************/

LegacyOptimizationTools::ScaleFunction::ScaleFunction(
    shared_ptr<TreeLikelihoodInterface> tl) :
  tl_(tl),
  brLen_(),
  lambda_()
{
  // We work only on the branch lengths:
  brLen_ = tl->getBranchLengthsParameters();
  if (brLen_.hasParameter("RootPosition"))
    brLen_.deleteParameter("RootPosition");
  lambda_.addParameter(Parameter("scale factor", 0));
}

LegacyOptimizationTools::ScaleFunction::~ScaleFunction() {}

void LegacyOptimizationTools::ScaleFunction::setParameters(const ParameterList& lambda)
{
  if (lambda.size() != 1)
    throw Exception("LegacyOptimizationTools::ScaleFunction::f(). This is a one parameter function!");
  lambda_.setParametersValues(lambda);
}

double LegacyOptimizationTools::ScaleFunction::getValue() const
{
  // Scale the tree:
  ParameterList brLen = brLen_;
  double s = exp(lambda_[0].getValue());
  for (unsigned int i = 0; i < brLen.size(); i++)
  {
    try
    {
      brLen[i].setValue(brLen[i].getValue() * s);
    }
    catch (ConstraintException& cex)
    {
      // Do nothing. Branch value is already at bound...
    }
  }
  return tl_->f(brLen);
}

/******************************************************************************/

unsigned int LegacyOptimizationTools::optimizeTreeScale(
  shared_ptr<TreeLikelihoodInterface> tl,
  double tolerance,
  unsigned int tlEvalMax,
  shared_ptr<OutputStream> messageHandler,
  shared_ptr<OutputStream> profiler,
  unsigned int verbose)
{
  auto sf = make_shared<ScaleFunction>(tl);
  BrentOneDimension bod(sf);
  bod.setMessageHandler(messageHandler);
  bod.setProfiler(profiler);
  ParameterList singleParameter;
  singleParameter.addParameter(Parameter("scale factor", 0));
  bod.setInitialInterval(-0.5, 0.5);
  bod.init(singleParameter);
  auto PS =make_shared<ParametersStopCondition>(&bod, tolerance);
  bod.setStopCondition(PS);
  bod.setMaximumNumberOfEvaluations(tlEvalMax);
  bod.optimize();
  ApplicationTools::displayTaskDone();
  if (verbose > 0)
    ApplicationTools::displayResult("Tree scaled by", exp(sf->getParameters()[0].getValue()));
  return bod.getNumberOfEvaluations();
}

/******************************************************************************/

unsigned int LegacyOptimizationTools::optimizeNumericalParameters(
  shared_ptr<DiscreteRatesAcrossSitesTreeLikelihoodInterface> tl,
  const ParameterList& parameters,
  shared_ptr<OptimizationListener> listener,
  unsigned int nstep,
  double tolerance,
  unsigned int tlEvalMax,
  shared_ptr<OutputStream> messageHandler,
  shared_ptr<OutputStream> profiler,
  bool reparametrization,
  unsigned int verbose,
  const string& optMethodDeriv,
  const string& optMethodModel)
{
  shared_ptr<SecondOrderDerivable> f = tl;
  ParameterList pl = parameters;

  // Shall we reparametrize the function to remove constraints?
  if (reparametrization)
  {
    f = make_shared<ReparametrizationDerivableSecondOrderWrapper>(f, parameters);

    // Reset parameters to remove constraints:
    pl = f->getParameters().createSubList(parameters.getParameterNames());
  }

  // ///////////////
  // Build optimizer:

  // Branch lengths

  auto desc = make_unique<MetaOptimizerInfos>();
  unique_ptr<MetaOptimizer> poptimizer;
  shared_ptr<AbstractNumericalDerivative> fnum(new ThreePointsNumericalDerivative(f));

  if (optMethodDeriv == OPTIMIZATION_GRADIENT)
    desc->addOptimizer("Branch length parameters", make_shared<ConjugateGradientMultiDimensions>(f), tl->getBranchLengthsParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
  else if (optMethodDeriv == OPTIMIZATION_NEWTON)
    desc->addOptimizer("Branch length parameters", make_shared<PseudoNewtonOptimizer>(f), tl->getBranchLengthsParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
  else if (optMethodDeriv == OPTIMIZATION_BFGS)
    desc->addOptimizer("Branch length parameters", make_shared<BfgsMultiDimensions>(f), tl->getBranchLengthsParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
  else
    throw Exception("LegacyOptimizationTools::optimizeNumericalParameters. Unknown optimization method: " + optMethodDeriv);

  // Other parameters

  if (optMethodModel == OPTIMIZATION_BRENT)
  {
    ParameterList plsm = parameters.getCommonParametersWith(tl->getSubstitutionModelParameters());
    desc->addOptimizer("Substitution model parameter", make_shared<SimpleMultiDimensions>(f), plsm.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);


    ParameterList plrd = parameters.getCommonParametersWith(tl->getRateDistributionParameters());
    desc->addOptimizer("Rate distribution parameter", make_shared<SimpleMultiDimensions>(f), plrd.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);
    poptimizer.reset(new MetaOptimizer(f, move(desc), nstep));
  }
  else if (optMethodModel == OPTIMIZATION_BFGS)
  {
    vector<string> vNameDer;

    ParameterList plsm = parameters.getCommonParametersWith(tl->getSubstitutionModelParameters());
    vNameDer = plsm.getParameterNames();

    ParameterList plrd = parameters.getCommonParametersWith(tl->getRateDistributionParameters());

    vector<string> vNameDer2 = plrd.getParameterNames();

    vNameDer.insert(vNameDer.begin(), vNameDer2.begin(), vNameDer2.end());
    fnum->setParametersToDerivate(vNameDer);

    desc->addOptimizer("Rate & model distribution parameters", make_shared<BfgsMultiDimensions>(fnum), vNameDer, 1, MetaOptimizerInfos::IT_TYPE_FULL);
    poptimizer.reset(new MetaOptimizer(fnum, move(desc), nstep));
  }
  else
    throw Exception("LegacyOptimizationTools::optimizeNumericalParameters. Unknown optimization method: " + optMethodModel);

  poptimizer->setVerbose(verbose);
  poptimizer->setProfiler(profiler);
  poptimizer->setMessageHandler(messageHandler);
  poptimizer->setMaximumNumberOfEvaluations(tlEvalMax);
  poptimizer->getStopCondition()->setTolerance(tolerance);

  // Optimize TreeLikelihood function:
  poptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  auto nanListener = make_shared<NaNListener>(poptimizer.get(), tl.get());
  poptimizer->addOptimizationListener(nanListener);
  if (listener)
    poptimizer->addOptimizationListener(listener);
  poptimizer->init(pl);
  poptimizer->optimize();

  if (verbose > 0)
    ApplicationTools::displayMessage("\n");

  // We're done.
  unsigned int nb = poptimizer->getNumberOfEvaluations();
  return nb;
}


/******************************************************************************/

unsigned int LegacyOptimizationTools::optimizeNumericalParameters2(
  shared_ptr<DiscreteRatesAcrossSitesTreeLikelihoodInterface> tl,
  const ParameterList& parameters,
  shared_ptr<OptimizationListener> listener,
  double tolerance,
  unsigned int tlEvalMax,
  shared_ptr<OutputStream> messageHandler,
  shared_ptr<OutputStream> profiler,
  bool reparametrization,
  bool useClock,
  unsigned int verbose,
  const std::string& optMethodDeriv)
{
  shared_ptr<SecondOrderDerivable> f = tl;
  shared_ptr<GlobalClockTreeLikelihoodFunctionWrapper> fclock;

  ParameterList pl = parameters;
  // Shall we use a molecular clock constraint on branch lengths?
  if (useClock)
  {
    fclock.reset(new GlobalClockTreeLikelihoodFunctionWrapper(tl));
    f = fclock;
    if (verbose > 0)
      ApplicationTools::displayResult("Log-likelihood after adding clock", -(tl->getLogLikelihood()));

    // Reset parameters to use new branch lengths. WARNING! 'old' branch parameters do not exist anymore and have been replaced by heights
    pl = fclock->getParameters().getCommonParametersWith(parameters);
    pl.addParameters(fclock->getHeightParameters());
  }
  // Shall we reparametrize the function to remove constraints?
  if (reparametrization)
  {
    f = make_shared<ReparametrizationDerivableSecondOrderWrapper>(f, pl);

    // Reset parameters to remove constraints:
    pl = f->getParameters().createSubList(pl.getParameterNames());
  }

  shared_ptr<AbstractNumericalDerivative> fnum;
  
  // Build optimizer:
  unique_ptr<OptimizerInterface> optimizer;
  if (optMethodDeriv == OPTIMIZATION_GRADIENT)
  {
    fnum.reset(new TwoPointsNumericalDerivative(dynamic_pointer_cast<FirstOrderDerivable>(f)));
    fnum->setInterval(0.0000001);
    optimizer.reset(new ConjugateGradientMultiDimensions(fnum));
  }
  else if (optMethodDeriv == OPTIMIZATION_NEWTON)
  {
    fnum.reset(new ThreePointsNumericalDerivative(f));
    fnum->setInterval(0.0001);
    optimizer.reset(new PseudoNewtonOptimizer(fnum));
  }
  else if (optMethodDeriv == OPTIMIZATION_BFGS)
  {
    fnum.reset(new TwoPointsNumericalDerivative(dynamic_pointer_cast<FirstOrderDerivable>(f)));
    fnum->setInterval(0.0001);
    optimizer.reset(new BfgsMultiDimensions(fnum));
  }
  else
    throw Exception("LegacyOptimizationTools::optimizeNumericalParameters2. Unknown optimization method: " + optMethodDeriv);

  // Numerical derivatives:
  ParameterList tmp = tl->getNonDerivableParameters();

  if (useClock)
    tmp.addParameters(fclock->getHeightParameters());
  fnum->setParametersToDerivate(tmp.getParameterNames());
  optimizer->setVerbose(verbose);
  optimizer->setProfiler(profiler);
  optimizer->setMessageHandler(messageHandler);
  optimizer->setMaximumNumberOfEvaluations(tlEvalMax);
  optimizer->getStopCondition()->setTolerance(tolerance);

  // Optimize TreeLikelihood function:
  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  auto nanListener = make_shared<NaNListener>(optimizer.get(), tl.get());
  optimizer->addOptimizationListener(nanListener);
  if (listener)
    optimizer->addOptimizationListener(listener);
  optimizer->init(pl);
  optimizer->optimize();

  if (verbose > 0)
    ApplicationTools::displayMessage("\n");

  // We're done.
  return optimizer->getNumberOfEvaluations();
}


/******************************************************************************/

unsigned int LegacyOptimizationTools::optimizeBranchLengthsParameters(
  shared_ptr<DiscreteRatesAcrossSitesTreeLikelihoodInterface> tl,
  const ParameterList& parameters,
  shared_ptr<OptimizationListener> listener,
  double tolerance,
  unsigned int tlEvalMax,
  shared_ptr<OutputStream> messageHandler,
  shared_ptr<OutputStream> profiler,
  unsigned int verbose,
  const string& optMethodDeriv)
{
  // Build optimizer:
  unique_ptr<OptimizerInterface> optimizer;
  if (optMethodDeriv == OPTIMIZATION_GRADIENT)
  {
    tl->enableFirstOrderDerivatives(true);
    tl->enableSecondOrderDerivatives(false);
    optimizer.reset(new ConjugateGradientMultiDimensions(tl));
  }
  else if (optMethodDeriv == OPTIMIZATION_NEWTON)
  {
    tl->enableFirstOrderDerivatives(true);
    tl->enableSecondOrderDerivatives(true);
    optimizer.reset(new PseudoNewtonOptimizer(tl));
  }
  else if (optMethodDeriv == OPTIMIZATION_BFGS)
  {
    tl->enableFirstOrderDerivatives(true);
    tl->enableSecondOrderDerivatives(false);
    optimizer.reset(new BfgsMultiDimensions(tl));
  }
  else
    throw Exception("LegacyOptimizationTools::optimizeBranchLengthsParameters. Unknown optimization method: " + optMethodDeriv);
  optimizer->setVerbose(verbose);
  optimizer->setProfiler(profiler);
  optimizer->setMessageHandler(messageHandler);
  optimizer->setMaximumNumberOfEvaluations(tlEvalMax);
  optimizer->getStopCondition()->setTolerance(tolerance);

  // Optimize TreeLikelihood function:
  ParameterList pl = parameters.getCommonParametersWith(tl->getBranchLengthsParameters());
  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  if (listener)
    optimizer->addOptimizationListener(listener);
  optimizer->init(pl);
  optimizer->optimize();
  if (verbose > 0)
    ApplicationTools::displayMessage("\n");

  // We're done.
  unsigned int n = optimizer->getNumberOfEvaluations();
  return n;
}

/******************************************************************************/

void NNITopologyListener::topologyChangeSuccessful(const TopologyChangeEvent& event)
{
  optimizeCounter_++;
  if (optimizeCounter_ == optimizeNumerical_)
  {
    auto likelihood = dynamic_pointer_cast<DiscreteRatesAcrossSitesTreeLikelihoodInterface>(topoSearch_->getSearchableObject());
    parameters_.matchParametersValues(likelihood->getParameters());
    LegacyOptimizationTools::optimizeNumericalParameters(likelihood, parameters_, 0, nStep_, tolerance_, 1000000, messenger_, profiler_, reparametrization_, verbose_, optMethod_);
    optimizeCounter_ = 0;
  }
}

/******************************************************************************/

void NNITopologyListener2::topologyChangeSuccessful(const TopologyChangeEvent& event)
{
  optimizeCounter_++;
  if (optimizeCounter_ == optimizeNumerical_)
  {
    auto likelihood = dynamic_pointer_cast<DiscreteRatesAcrossSitesTreeLikelihoodInterface>(topoSearch_->getSearchableObject());
    parameters_.matchParametersValues(likelihood->getParameters());
    LegacyOptimizationTools::optimizeNumericalParameters2(likelihood, parameters_, 0, tolerance_, 1000000, messenger_, profiler_, reparametrization_, false, verbose_, optMethod_);
    optimizeCounter_ = 0;
  }
}

// ******************************************************************************/

shared_ptr<NNIHomogeneousTreeLikelihood> LegacyOptimizationTools::optimizeTreeNNI(
    shared_ptr<NNIHomogeneousTreeLikelihood> tl,
    const ParameterList& parameters,
    bool optimizeNumFirst,
    double tolBefore,
    double tolDuring,
    unsigned int tlEvalMax,
    unsigned int numStep,
    shared_ptr<OutputStream> messageHandler,
    shared_ptr<OutputStream> profiler,
    bool reparametrization,
    unsigned int verbose,
    const string& optMethodDeriv,
    unsigned int nStep,
    const string& nniMethod)
{
  // Roughly optimize parameter
  if (optimizeNumFirst)
  {
    LegacyOptimizationTools::optimizeNumericalParameters(tl, parameters, NULL, nStep, tolBefore, 1000000, messageHandler, profiler, reparametrization, verbose, optMethodDeriv);
  }
  // Begin topo search:
  auto topoSearch = make_shared<NNITopologySearch>(tl, nniMethod, verbose > 2 ? verbose - 2 : 0);
  auto topoListener = make_shared<NNITopologyListener>(topoSearch, parameters, tolDuring, messageHandler, profiler, verbose, optMethodDeriv, nStep, reparametrization);
  topoListener->setNumericalOptimizationCounter(numStep);
  topoSearch->addTopologyListener(topoListener);
  topoSearch->search();
  return dynamic_pointer_cast<NNIHomogeneousTreeLikelihood>(topoSearch->getSearchableObject());
}

/******************************************************************************/

shared_ptr<NNIHomogeneousTreeLikelihood> LegacyOptimizationTools::optimizeTreeNNI2(
    shared_ptr<NNIHomogeneousTreeLikelihood> tl,
    const ParameterList& parameters,
    bool optimizeNumFirst,
    double tolBefore,
    double tolDuring,
    unsigned int tlEvalMax,
    unsigned int numStep,
    shared_ptr<OutputStream> messageHandler,
    shared_ptr<OutputStream> profiler,
    bool reparametrization,
    unsigned int verbose,
    const string& optMethodDeriv,
    const string& nniMethod)
{
  // Roughly optimize parameter
  if (optimizeNumFirst)
  {
    LegacyOptimizationTools::optimizeNumericalParameters2(tl, parameters, nullptr, tolBefore, 1000000, messageHandler, profiler, reparametrization, false, verbose, optMethodDeriv);
  }
  // Begin topo search:
  auto topoSearch = make_shared<NNITopologySearch> (tl, nniMethod, verbose > 2 ? verbose - 2 : 0);
  auto topoListener = make_shared<NNITopologyListener2>(topoSearch, parameters, tolDuring, messageHandler, profiler, verbose, optMethodDeriv, reparametrization);
  topoListener->setNumericalOptimizationCounter(numStep);
  topoSearch->addTopologyListener(topoListener);
  topoSearch->search();
  return dynamic_pointer_cast<NNIHomogeneousTreeLikelihood>(topoSearch->getSearchableObject());
}

/******************************************************************************/

shared_ptr<DRTreeParsimonyScore> LegacyOptimizationTools::optimizeTreeNNI(
    shared_ptr<DRTreeParsimonyScore> tp,
    unsigned int verbose)
{
  auto topo = dynamic_pointer_cast<NNISearchable>(tp);
  NNITopologySearch topoSearch(topo, NNITopologySearch::PHYML, verbose);
  topoSearch.search();
  return dynamic_pointer_cast<DRTreeParsimonyScore>(topoSearch.getSearchableObject());
}


