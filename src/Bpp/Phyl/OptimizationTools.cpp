// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Function/MetaOptimizer.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/SimpleMultiDimensions.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/TwoPointsNumericalDerivative.h>
#include <Bpp/Numeric/ParameterList.h>

#include "Io/Newick.h"
#include "OptimizationTools.h"
#include "PseudoNewtonOptimizer.h"
#include "Tree/PhyloTreeTools.h"

// From bpp-seq:
#include <Bpp/Seq/Io/Fasta.h>

#include <memory>

using namespace bpp;
using namespace std;

/******************************************************************************/

OptimizationTools::OptimizationTools() {}

OptimizationTools::~OptimizationTools() {}

/******************************************************************************/

std::string OptimizationTools::OPTIMIZATION_NEWTON = "newton";
std::string OptimizationTools::OPTIMIZATION_GRADIENT = "gradient";
std::string OptimizationTools::OPTIMIZATION_BRENT = "Brent";
std::string OptimizationTools::OPTIMIZATION_BFGS = "BFGS";

/******************************************************************************/

unsigned int OptimizationTools::optimizeNumericalParameters(
  shared_ptr<PhyloLikelihoodInterface> lik,
  const ParameterList& parameters,
  shared_ptr<OptimizationListener> listener,
  unsigned int nstep,
  double tolerance,
  unsigned int tlEvalMax,
  shared_ptr<OutputStream> messageHandler,
  shared_ptr<OutputStream> profiler,
  bool reparametrization,
  unsigned int verbose,
  const std::string& optMethodDeriv,
  const std::string& optMethodModel)
{
  shared_ptr<SecondOrderDerivable> f = lik;
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

  if (optMethodDeriv == OPTIMIZATION_GRADIENT)
    desc->addOptimizer("Branch length parameters", make_shared<ConjugateGradientMultiDimensions>(f), lik->getBranchLengthParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
  else if (optMethodDeriv == OPTIMIZATION_NEWTON)
    desc->addOptimizer("Branch length parameters", make_shared<PseudoNewtonOptimizer>(f), lik->getBranchLengthParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
  else if (optMethodDeriv == OPTIMIZATION_BFGS)
    desc->addOptimizer("Branch length parameters", make_shared<BfgsMultiDimensions>(f), lik->getBranchLengthParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
  else
    throw Exception("OptimizationTools::optimizeNumericalParameters. Unknown optimization method: " + optMethodDeriv);

  // Other parameters

  if (optMethodModel == OPTIMIZATION_BRENT)
  {
    ParameterList plTmp = lik->getSubstitutionModelParameters();
    plTmp.addParameters(lik->getRootFrequenciesParameters());
    ParameterList plsm = parameters.getCommonParametersWith(plTmp);
    desc->addOptimizer("Substitution model parameters", make_shared<SimpleMultiDimensions>(f), plsm.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);


    ParameterList plrd = parameters.getCommonParametersWith(lik->getRateDistributionParameters());
    desc->addOptimizer("Rate distribution parameters", make_shared<SimpleMultiDimensions>(f), plrd.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);
    poptimizer = make_unique<MetaOptimizer>(f, std::move(desc), nstep);
  }
  else if (optMethodModel == OPTIMIZATION_BFGS)
  {
    vector<string> vNameDer;
    auto fnum = make_shared<ThreePointsNumericalDerivative>(f);

    ParameterList plsm = parameters.getCommonParametersWith(lik->getSubstitutionModelParameters());
    vNameDer = plsm.getParameterNames();

    ParameterList plrd = parameters.getCommonParametersWith(lik->getRateDistributionParameters());

    vector<string> vNameDer2 = plrd.getParameterNames();

    vNameDer.insert(vNameDer.begin(), vNameDer2.begin(), vNameDer2.end());
    fnum->setParametersToDerivate(vNameDer);

    desc->addOptimizer("Rate & model distribution parameters", make_shared<BfgsMultiDimensions>(fnum), vNameDer, 1, MetaOptimizerInfos::IT_TYPE_FULL);
    poptimizer = make_unique<MetaOptimizer>(fnum, std::move(desc), nstep);
  }
  else
    throw Exception("OptimizationTools::optimizeNumericalParameters. Unknown optimization method: " + optMethodModel);

  poptimizer->setVerbose(verbose);
  poptimizer->setProfiler(profiler);
  poptimizer->setMessageHandler(messageHandler);
  poptimizer->setMaximumNumberOfEvaluations(tlEvalMax);
  poptimizer->getStopCondition()->setTolerance(tolerance);

  // Optimize TreeLikelihood function:
  poptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  auto nanListener = make_shared<NaNListener>(poptimizer.get(), lik.get());
  poptimizer->addOptimizationListener(nanListener);
  if (listener)
    poptimizer->addOptimizationListener(listener);
  poptimizer->init(pl);
  poptimizer->optimize();

  if (verbose > 0)
    ApplicationTools::displayMessage("\n");

  // We're done.
  uint nb = poptimizer->getNumberOfEvaluations();
  return nb;
}


/************************************************************/

unsigned int OptimizationTools::optimizeNumericalParameters2(
  shared_ptr<PhyloLikelihoodInterface> lik,
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
  shared_ptr<SecondOrderDerivable> f = lik;
  ParameterList pl = parameters;

  // Shall we use a molecular clock constraint on branch lengths?
  // unique_ptr<GlobalClockTreeLikelihoodFunctionWrapper> fclock;
  // if (useClock)
  //   {
  //     fclock.reset(new GlobalClockTreeLikelihoodFunctionWrapper(lik));
  //     f = fclock.get();
  //     if (verbose > 0)
  //       ApplicationTools::displayResult("Log-likelihood after adding clock", -lik->getLogLikelihood());

  //     // Reset parameters to use new branch lengths. WARNING! 'old' branch parameters do not exist anymore and have been replaced by heights
  //     pl = fclock->getParameters().getCommonParametersWith(parameters);
  //     pl.addParameters(fclock->getHeightParameters());
  //   }

  // Shall we reparametrize the function to remove constraints?
  shared_ptr<SecondOrderDerivable> frep;
  if (reparametrization)
  {
    frep.reset(new ReparametrizationDerivableSecondOrderWrapper(f, pl));
    f = frep;

    // Reset parameters to remove constraints:
    pl = f->getParameters().createSubList(pl.getParameterNames());
  }

  // Build optimizer:
  unique_ptr<OptimizerInterface> optimizer;
  shared_ptr<AbstractNumericalDerivative> fnum;

  if (optMethodDeriv == OPTIMIZATION_GRADIENT)
  {
    fnum = make_shared<TwoPointsNumericalDerivative>(dynamic_pointer_cast<FirstOrderDerivable>(f));
    fnum->setInterval(0.0000001);
    optimizer = make_unique<ConjugateGradientMultiDimensions>(fnum);
  }
  else if (optMethodDeriv == OPTIMIZATION_NEWTON)
  {
    fnum = make_shared<ThreePointsNumericalDerivative>(f);
    fnum->setInterval(0.0001);
    optimizer = make_unique<PseudoNewtonOptimizer>(fnum);
  }
  else if (optMethodDeriv == OPTIMIZATION_BFGS)
  {
    fnum = make_shared<TwoPointsNumericalDerivative>(dynamic_pointer_cast<FirstOrderDerivable>(f));
    fnum->setInterval(0.0001);
    optimizer = make_unique<BfgsMultiDimensions>(fnum);
  }
  else
    throw Exception("OptimizationTools::optimizeNumericalParameters2. Unknown optimization method: " + optMethodDeriv);

  // Numerical derivatives:
  // Variables not derivatived in Likelihood DF but in numerical way
  ParameterList tmp = lik->getParameters();

  // if (useClock)
  //   tmp.addParameters(fclock->getHeightParameters());

  fnum->setParametersToDerivate(tmp.getParameterNames());

  optimizer->setVerbose(verbose);
  optimizer->setProfiler(profiler);
  optimizer->setMessageHandler(messageHandler);
  optimizer->setMaximumNumberOfEvaluations(tlEvalMax);
  optimizer->getStopCondition()->setTolerance(tolerance);

  // Optimize TreeLikelihood function:
  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  auto nanListener = make_shared<NaNListener>(optimizer.get(), lik.get());
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

/************************************************************/

unsigned int OptimizationTools::optimizeNumericalParameters2(
  shared_ptr<SingleProcessPhyloLikelihood> lik,
  const ParameterList& parameters,
  shared_ptr<OptimizationListener> listener,
  double tolerance,
  unsigned int tlEvalMax,
  shared_ptr<OutputStream> messageHandler,
  shared_ptr<OutputStream> profiler,
  bool reparametrization,
  bool useClock,
  unsigned int verbose,
  const string& optMethodDeriv)
{
  shared_ptr<SecondOrderDerivable> f = lik;
  ParameterList pl = parameters;
  if (reparametrization)
  {
    // Shall we reparametrize the function to remove constraints?
    if (reparametrization)
    {
      f = make_shared<ReparametrizationDerivableSecondOrderWrapper>(f, parameters);

      // Reset parameters to remove constraints:
      pl = f->getParameters().createSubList(parameters.getParameterNames());
    }
  }

  // Build optimizer:
  shared_ptr<AbstractNumericalDerivative> fnum;
  unique_ptr<OptimizerInterface> optimizer;

  if (optMethodDeriv == OPTIMIZATION_GRADIENT)
  {
    lik->likelihoodCalculationSingleProcess().setNumericalDerivateConfiguration(0.00001, NumericalDerivativeType::ThreePoints);

    fnum = make_shared<TwoPointsNumericalDerivative>(dynamic_pointer_cast<FirstOrderDerivable>(f));
    fnum->setInterval(0.0000001);
    optimizer.reset(new ConjugateGradientMultiDimensions(fnum));
  }
  else if (optMethodDeriv == OPTIMIZATION_NEWTON)
  {
    lik->likelihoodCalculationSingleProcess().setNumericalDerivateConfiguration(0.0001, NumericalDerivativeType::FivePoints);
    fnum = make_shared<ThreePointsNumericalDerivative>(f);
    fnum->setInterval(0.0001);
    optimizer.reset(new PseudoNewtonOptimizer(fnum));
  }
  else if (optMethodDeriv == OPTIMIZATION_BFGS)
  {
    lik->likelihoodCalculationSingleProcess().setNumericalDerivateConfiguration(0.0001, NumericalDerivativeType::ThreePoints);
    fnum = make_shared<TwoPointsNumericalDerivative>(dynamic_pointer_cast<FirstOrderDerivable>(f));
    fnum->setInterval(0.0001);
    optimizer.reset(new BfgsMultiDimensions(fnum));
  }
  else
    throw Exception("OptimizationTools::optimizeNumericalParameters2. Unknown optimization method: " + optMethodDeriv);


  // Variables not derivatived in Likelihood DF but in numerical way
  ParameterList tmp = f->getParameters();

  fnum->setParametersToDerivate(tmp.getParameterNames());

  optimizer->setVerbose(verbose);
  optimizer->setProfiler(profiler);
  optimizer->setMessageHandler(messageHandler);
  optimizer->setMaximumNumberOfEvaluations(tlEvalMax);
  optimizer->getStopCondition()->setTolerance(tolerance);

  // Optimize TreeLikelihood function:
  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  auto nanListener = make_shared<NaNListener>(optimizer.get(), lik.get());
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

std::string OptimizationTools::DISTANCEMETHOD_INIT       = "init";
std::string OptimizationTools::DISTANCEMETHOD_PAIRWISE   = "pairwise";
std::string OptimizationTools::DISTANCEMETHOD_ITERATIONS = "iterations";

/******************************************************************************/

unique_ptr<DistanceMatrix> OptimizationTools::estimateDistanceMatrix(
  DistanceEstimation& estimationMethod,
  const ParameterList& parametersToIgnore,
  const std::string& param,
  unsigned int verbose)
{
  if (param != DISTANCEMETHOD_PAIRWISE && param != DISTANCEMETHOD_INIT)
    throw Exception("OptimizationTools::estimateDistanceMatrix. Invalid option param=" + param + ".");
  estimationMethod.resetAdditionalParameters();
  estimationMethod.setVerbose(verbose);
  if (param == DISTANCEMETHOD_PAIRWISE)
  {
    ParameterList tmp = estimationMethod.model().getIndependentParameters();
    tmp.addParameters(estimationMethod.rateDistribution().getIndependentParameters());
    tmp.deleteParameters(parametersToIgnore.getParameterNames());
    estimationMethod.setAdditionalParameters(tmp);
  }
  // Compute matrice:
  if (verbose > 0)
    ApplicationTools::displayTask("Estimating distance matrix", true);
  estimationMethod.computeMatrix();
  auto matrix = estimationMethod.getMatrix();
  
  if (verbose > 0)
    ApplicationTools::displayTaskDone();

  return matrix;
}

/******************************************************************************/

unique_ptr<TreeTemplate<Node>> OptimizationTools::buildDistanceTree(
  DistanceEstimation& estimationMethod,
  AgglomerativeDistanceMethodInterface& reconstructionMethod,
  const ParameterList& parametersToIgnore,
  bool optimizeBrLen,
  const std::string& param,
  double tolerance,
  unsigned int tlEvalMax,
  shared_ptr<OutputStream> profiler,
  shared_ptr<OutputStream> messenger,
  unsigned int verbose)
{
  estimationMethod.resetAdditionalParameters();
  estimationMethod.setVerbose(verbose);
  if (param == DISTANCEMETHOD_PAIRWISE)
  {
    ParameterList tmp = estimationMethod.model().getIndependentParameters();
    tmp.addParameters(estimationMethod.rateDistribution().getIndependentParameters());
    tmp.deleteParameters(parametersToIgnore.getParameterNames());
    estimationMethod.setAdditionalParameters(tmp);
  }
  unique_ptr<TreeTemplate<Node>> tree = nullptr;
  unique_ptr<TreeTemplate<Node>> previousTree = nullptr;
  bool test = true;
  while (test)
  {
    // Compute matrice:
    if (verbose > 0)
      ApplicationTools::displayTask("Estimating distance matrix", true);
    estimationMethod.computeMatrix();
    auto matrix = estimationMethod.getMatrix();
    if (verbose > 0)
      ApplicationTools::displayTaskDone();

    // Compute tree:
    if (matrix->size() == 2)
    {
      // Special case, there is only one possible tree:
      Node* n1 = new Node(0);
      Node* n2 = new Node(1, matrix->getName(0));
      n2->setDistanceToFather((*matrix)(0, 0) / 2.);
      Node* n3 = new Node(2, matrix->getName(1));
      n3->setDistanceToFather((*matrix)(0, 0) / 2.);
      n1->addSon(n2);
      n1->addSon(n3);
      tree.reset(new TreeTemplate<Node>(n1));
      break;
    }
    /* For future integration
       shared_ptr<PhyloTree> tree;
       shared_ptr<PhyloTree> previousTree;
       bool test = true;
       while (test)
       {
       // Compute matrice:
       if (verbose > 0) 
       ApplicationTools::displayTask("Estimating distance matrix", true);
       estimationMethod.computeMatrix();
       DistanceMatrix* matrix = estimationMethod.getMatrix();
       if (verbose > 0)
       ApplicationTools::displayTaskDone();
       
       // Compute tree:
       if (matrix->size() == 2)
       {
       // Special case, there is only one possible tree:
       auto root=std::make_shared<PhyloNode>();
       tree.createNode(root);
       tree.setNodeIndex(root, 0);
       auto n1=std::make_shared<PhyloNode>(matrix->getName(0));
       auto n2=std::make_shared<PhyloNode>(matrix->getName(1));
       auto branch1=shared_ptr<PhyloBranch> ((*matrix)(0, 0) / 2.);
       auto branch2=shared_ptr<PhyloBranch> ((*matrix)(0, 0) / 2.);
       tree.createNode(root, n1, branch1);
       tree.createNode(root, n2, branch2);
       tree.setNodeIndex(n1, 1);
       tree.setNodeIndex(n2, 2);
       break;
       }
    */
    if (verbose > 0)
      ApplicationTools::displayTask("Building tree");
    reconstructionMethod.setDistanceMatrix(*matrix);
    reconstructionMethod.computeTree();
    previousTree = std::move(tree);

    tree = make_unique<TreeTemplate<Node>>(reconstructionMethod.tree());
    if (verbose > 0)
      ApplicationTools::displayTaskDone();
    if (previousTree && verbose > 0)
    {
      int rf = TreeTools::robinsonFouldsDistance(*previousTree, *tree, false);
      ApplicationTools::displayResult("Topo. distance with previous iteration", TextTools::toString(rf));
      test = (rf == 0);
    }
    if (param != DISTANCEMETHOD_ITERATIONS)
      break;                              // Ends here.

    // Now, re-estimate parameters:
    Context context;
  
    shared_ptr<BranchModelInterface> model(estimationMethod.model().clone());
    shared_ptr<DiscreteDistributionInterface> rdist(estimationMethod.rateDistribution().clone());
    auto phyloT  = PhyloTreeTools::buildFromTreeTemplate(*tree);
    auto process = make_shared<RateAcrossSitesSubstitutionProcess>(model, rdist, phyloT);
    auto lik     = make_shared<LikelihoodCalculationSingleProcess>(context, estimationMethod.getData(), process);
    auto tl      = make_shared<SingleProcessPhyloLikelihood>(context, lik);

    ParameterList parameters = tl->getParameters();
    if (!optimizeBrLen)
    {
      vector<string> vs = tl->getBranchLengthParameters().getParameterNames();
      parameters.deleteParameters(vs);
    }
    parameters.deleteParameters(parametersToIgnore.getParameterNames());
    optimizeNumericalParameters(tl, parameters, NULL, 0, tolerance, tlEvalMax, messenger, profiler, verbose > 0 ? verbose - 1 : 0);
    if (verbose > 0)
    {
      ParameterList tmp = tl->getSubstitutionModelParameters();
      for (unsigned int i = 0; i < tmp.size(); ++i)
      {
        ApplicationTools::displayResult(tmp[i].getName(), TextTools::toString(tmp[i].getValue()));
      }
      tmp = tl->getRateDistributionParameters();
      for (unsigned int i = 0; i < tmp.size(); ++i)
      {
        ApplicationTools::displayResult(tmp[i].getName(), TextTools::toString(tmp[i].getValue()));
      }
    }
  }
  return tree;
}

/******************************************************************************/
