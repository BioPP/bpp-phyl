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
  PhyloLikelihood* lik,
  const ParameterList& parameters,
  OptimizationListener* listener,
  unsigned int nstep,
  double tolerance,
  unsigned int tlEvalMax,
  OutputStream* messageHandler,
  OutputStream* profiler,
  bool reparametrization,
  unsigned int verbose,
  const std::string& optMethodDeriv,
  const std::string& optMethodModel)
{
  DerivableSecondOrder* f = lik;
  ParameterList pl = parameters;

  // Shall we reparametrize the function to remove constraints?
  unique_ptr<DerivableSecondOrder> frep;
  if (reparametrization)
  {
    frep.reset(new ReparametrizationDerivableSecondOrderWrapper(f, parameters));
    f = frep.get();

    // Reset parameters to remove constraints:
    pl = f->getParameters().createSubList(parameters.getParameterNames());
  }

  // ///////////////
  // Build optimizer:

  // Branch lengths

  MetaOptimizerInfos* desc = new MetaOptimizerInfos();
  MetaOptimizer* poptimizer = 0;
  AbstractNumericalDerivative* fnum = new ThreePointsNumericalDerivative(f);

  if (optMethodDeriv == OPTIMIZATION_GRADIENT)
    desc->addOptimizer("Branch length parameters", new ConjugateGradientMultiDimensions(f), lik->getBranchLengthParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
  else if (optMethodDeriv == OPTIMIZATION_NEWTON)
    desc->addOptimizer("Branch length parameters", new PseudoNewtonOptimizer(f), lik->getBranchLengthParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
  else if (optMethodDeriv == OPTIMIZATION_BFGS)
    desc->addOptimizer("Branch length parameters", new BfgsMultiDimensions(f), lik->getBranchLengthParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
  else
    throw Exception("OptimizationTools::optimizeNumericalParameters. Unknown optimization method: " + optMethodDeriv);

  // Other parameters

  if (optMethodModel == OPTIMIZATION_BRENT)
  {
    ParameterList plsm = parameters.getCommonParametersWith(lik->getSubstitutionModelParameters());
    desc->addOptimizer("Substitution model parameter", new SimpleMultiDimensions(f), plsm.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);


    ParameterList plrd = parameters.getCommonParametersWith(lik->getRateDistributionParameters());
    desc->addOptimizer("Rate distribution parameter", new SimpleMultiDimensions(f), plrd.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);
    poptimizer = new MetaOptimizer(f, desc, nstep);
  }
  else if (optMethodModel == OPTIMIZATION_BFGS)
  {
    vector<string> vNameDer;

    ParameterList plsm = parameters.getCommonParametersWith(lik->getSubstitutionModelParameters());
    vNameDer = plsm.getParameterNames();

    ParameterList plrd = parameters.getCommonParametersWith(lik->getRateDistributionParameters());

    vector<string> vNameDer2 = plrd.getParameterNames();

    vNameDer.insert(vNameDer.begin(), vNameDer2.begin(), vNameDer2.end());
    fnum->setParametersToDerivate(vNameDer);

    desc->addOptimizer("Rate & model distribution parameters", new BfgsMultiDimensions(fnum), vNameDer, 1, MetaOptimizerInfos::IT_TYPE_FULL);
    poptimizer = new MetaOptimizer(fnum, desc, nstep);
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
  NaNListener* nanListener = new NaNListener(poptimizer, lik);
  poptimizer->addOptimizationListener(nanListener);
  if (listener)
    poptimizer->addOptimizationListener(listener);
  poptimizer->init(pl);
  poptimizer->optimize();

  if (verbose > 0)
    ApplicationTools::displayMessage("\n");

  // We're done.
  uint nb = poptimizer->getNumberOfEvaluations();
  delete poptimizer;
  return nb;
}


/************************************************************/

unsigned int OptimizationTools::optimizeNumericalParameters2(
  PhyloLikelihood* lik,
  const ParameterList& parameters,
  OptimizationListener* listener,
  double tolerance,
  unsigned int tlEvalMax,
  OutputStream* messageHandler,
  OutputStream* profiler,
  bool reparametrization,
  bool useClock,
  unsigned int verbose,
  const std::string& optMethodDeriv)
{
  DerivableSecondOrder* f = lik;
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
  // unique_ptr<DerivableSecondOrder> frep;
  // if (reparametrization)
  // {
  //   frep.reset(new ReparametrizationDerivableSecondOrderWrapper(f, pl));
  //   f = frep.get();

  //   // Reset parameters to remove constraints:
  //   pl = f->getParameters().createSubList(pl.getParameterNames());
  // }

  unique_ptr<AbstractNumericalDerivative> fnum;
  // Build optimizer:
  unique_ptr<Optimizer> optimizer;

  if (optMethodDeriv == OPTIMIZATION_GRADIENT)
  {
    fnum.reset(new TwoPointsNumericalDerivative(f));
    fnum->setInterval(0.0000001);
    optimizer.reset(new ConjugateGradientMultiDimensions(fnum.get())); // Removes strict-aliasing warning with gcc 4.4
  }
  else if (optMethodDeriv == OPTIMIZATION_NEWTON)
  {
    fnum.reset(new ThreePointsNumericalDerivative(f));
    fnum->setInterval(0.0001);
    optimizer.reset(new PseudoNewtonOptimizer(fnum.get()));
  }
  else if (optMethodDeriv == OPTIMIZATION_BFGS)
  {
    fnum.reset(new TwoPointsNumericalDerivative(f));
    fnum->setInterval(0.0001);
    optimizer.reset(new BfgsMultiDimensions(fnum.get()));
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
  NaNListener* nanListener = new NaNListener(optimizer.get(), lik);
  optimizer->addOptimizationListener(nanListener);
  if (listener)
    optimizer->addOptimizationListener(listener);

  optimizer->init(pl);

  optimizer->optimize();
  delete nanListener;

  if (verbose > 0)
    ApplicationTools::displayMessage("\n");

  // We're done.
  return optimizer->getNumberOfEvaluations();
}

/************************************************************/

unsigned int OptimizationTools::optimizeNumericalParameters2(
  SingleProcessPhyloLikelihood& lik,
  const ParameterList& parameters,
  OptimizationListener* listener,
  double tolerance,
  unsigned int tlEvalMax,
  OutputStream* messageHandler,
  OutputStream* profiler,
  bool reparametrization,
  bool useClock,
  unsigned int verbose,
  const std::string& optMethodDeriv)
{
  ParameterList pl = parameters;
  if (reparametrization)
  {
    throw Exception("OptimizationTools::optimizeNumericalParameters2 reparametrization not checked for dataflow likelihood calculation");
  }

  // Build optimizer:
  unique_ptr<AbstractNumericalDerivative> fnum;
  unique_ptr<Optimizer> optimizer;

  if (optMethodDeriv == OPTIMIZATION_GRADIENT)
  {
    lik.getLikelihoodCalculationSingleProcess()->setNumericalDerivateConfiguration(0.00001, NumericalDerivativeType::ThreePoints);

    fnum.reset(new TwoPointsNumericalDerivative(&lik));
    fnum->setInterval(0.0000001);
    optimizer.reset(new ConjugateGradientMultiDimensions(fnum.get())); // Removes strict-aliasing warning with gcc 4.4
  }
  else if (optMethodDeriv == OPTIMIZATION_NEWTON)
  {
    lik.getLikelihoodCalculationSingleProcess()->setNumericalDerivateConfiguration(0.0001, NumericalDerivativeType::FivePoints);
    fnum.reset(new ThreePointsNumericalDerivative(&lik));
    fnum->setInterval(0.0001);
    optimizer.reset(new PseudoNewtonOptimizer(fnum.get()));
  }
  else if (optMethodDeriv == OPTIMIZATION_BFGS)
  {
    lik.getLikelihoodCalculationSingleProcess()->setNumericalDerivateConfiguration(0.0001, NumericalDerivativeType::ThreePoints);
    fnum.reset(new TwoPointsNumericalDerivative(&lik));
    fnum->setInterval(0.0001);
    optimizer.reset(new BfgsMultiDimensions(fnum.get()));
  }
  else
    throw Exception("OptimizationTools::optimizeNumericalParameters2. Unknown optimization method: " + optMethodDeriv);


  // Variables not derivatived in Likelihood DF but in numerical way
  ParameterList tmp = lik.getParameters();

  fnum->setParametersToDerivate(tmp.getParameterNames());

  optimizer->setVerbose(verbose);
  optimizer->setProfiler(profiler);
  optimizer->setMessageHandler(messageHandler);
  optimizer->setMaximumNumberOfEvaluations(tlEvalMax);
  optimizer->getStopCondition()->setTolerance(tolerance);

  // Optimize TreeLikelihood function:
  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  NaNListener* nanListener = new NaNListener(optimizer.get(), &lik);
  optimizer->addOptimizationListener(nanListener);
  if (listener)
    optimizer->addOptimizationListener(listener);

  optimizer->init(pl);
  optimizer->optimize();

  delete nanListener;

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

DistanceMatrix* OptimizationTools::estimateDistanceMatrix(
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
    ParameterList tmp = estimationMethod.getModel().getIndependentParameters();
    tmp.addParameters(estimationMethod.getRateDistribution().getIndependentParameters());
    tmp.deleteParameters(parametersToIgnore.getParameterNames());
    estimationMethod.setAdditionalParameters(tmp);
  }
  // Compute matrice:
  if (verbose > 0)
    ApplicationTools::displayTask("Estimating distance matrix", true);
  estimationMethod.computeMatrix();
  unique_ptr<DistanceMatrix> matrix(estimationMethod.getMatrix());
  if (verbose > 0)
    ApplicationTools::displayTaskDone();

  return matrix.release();
}

/******************************************************************************/

TreeTemplate<Node>* OptimizationTools::buildDistanceTree(
  DistanceEstimation& estimationMethod,
  AgglomerativeDistanceMethod& reconstructionMethod,
  const ParameterList& parametersToIgnore,
  bool optimizeBrLen,
  const std::string& param,
  double tolerance,
  unsigned int tlEvalMax,
  OutputStream* profiler,
  OutputStream* messenger,
  unsigned int verbose)
{
  estimationMethod.resetAdditionalParameters();
  estimationMethod.setVerbose(verbose);
  if (param == DISTANCEMETHOD_PAIRWISE)
  {
    ParameterList tmp = estimationMethod.getModel().getIndependentParameters();
    tmp.addParameters(estimationMethod.getRateDistribution().getIndependentParameters());
    tmp.deleteParameters(parametersToIgnore.getParameterNames());
    estimationMethod.setAdditionalParameters(tmp);
  }
  TreeTemplate<Node>* tree = NULL;
  TreeTemplate<Node>* previousTree = NULL;
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
      Node* n1 = new Node(0);
      Node* n2 = new Node(1, matrix->getName(0));
      n2->setDistanceToFather((*matrix)(0, 0) / 2.);
      Node* n3 = new Node(2, matrix->getName(1));
      n3->setDistanceToFather((*matrix)(0, 0) / 2.);
      n1->addSon(n2);
      n1->addSon(n3);
      tree = new TreeTemplate<Node>(n1);
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
    previousTree = tree;
    delete matrix;

    tree = dynamic_cast<TreeTemplate<Node>*>(reconstructionMethod.getTree());
    if (verbose > 0)
      ApplicationTools::displayTaskDone();
    if (previousTree && verbose > 0)
    {
      int rf = TreeTools::robinsonFouldsDistance(*previousTree, *tree, false);
      ApplicationTools::displayResult("Topo. distance with previous iteration", TextTools::toString(rf));
      test = (rf == 0);
      delete previousTree;
    }
    if (param != DISTANCEMETHOD_ITERATIONS)
      break;                              // Ends here.

    // Now, re-estimate parameters:
    Context context;
  
    std::shared_ptr<BranchModel> model(estimationMethod.getModel().clone());
    std::shared_ptr<DiscreteDistribution> rdist(estimationMethod.getRateDistribution().clone());
    auto phyloT = PhyloTreeTools::buildFromTreeTemplate(*tree);
    ParametrizablePhyloTree partree(phyloT.get());

    auto process=std::make_shared<RateAcrossSitesSubstitutionProcess>(model, rdist, &partree);
      
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context,
                                                                    *estimationMethod.getData(), *process);

    SingleProcessPhyloLikelihood tl(context, lik);

    ParameterList parameters = tl.getParameters();
    if (!optimizeBrLen)
    {
      vector<string> vs = tl.getBranchLengthParameters().getParameterNames();
      parameters.deleteParameters(vs);
    }
    parameters.deleteParameters(parametersToIgnore.getParameterNames());
    optimizeNumericalParameters(&tl, parameters, NULL, 0, tolerance, tlEvalMax, messenger, profiler, verbose > 0 ? verbose - 1 : 0);
    if (verbose > 0)
    {
      ParameterList tmp = tl.getSubstitutionModelParameters();
      for (unsigned int i = 0; i < tmp.size(); i++)
      {
        ApplicationTools::displayResult(tmp[i].getName(), TextTools::toString(tmp[i].getValue()));
      }
      tmp = tl.getRateDistributionParameters();
      for (unsigned int i = 0; i < tmp.size(); i++)
      {
        ApplicationTools::displayResult(tmp[i].getName(), TextTools::toString(tmp[i].getValue()));
      }
    }
  }
  return tree;
}

/******************************************************************************/
