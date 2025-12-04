// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/KeyvalTools.h>

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
#include <Bpp/Numeric/Matrix/MatrixTools.h>

// bpp-phyl
#include "Likelihood/AutonomousSubstitutionProcess.h"
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

OptimizationTools::OptimizationOptions::OptimizationOptions(
    std::shared_ptr<PhyloLikelihoodInterface> lik,
    const std::map<std::string, std::string>& params,
    const std::string& suffix,
    bool suffixIsOptional,
    bool verb,
    int warn) :
  parameters(),
  listener(nullptr),
  backupFile(),
  nstep(1),
  tolerance(0.000001),
  nbEvalMax(1000000),
  messenger(ApplicationTools::message),
  profiler(ApplicationTools::message),
  reparametrization(false),
  useClock(0),
  verbose(1),
  optMethodDeriv(OPTIMIZATION_NEWTON),
  optMethodModel(OPTIMIZATION_BRENT)
{
  /// Get the method
  string optimization = ApplicationTools::getStringParameter("optimization", params, "FullD(derivatives=Newton)", suffix, suffixIsOptional, warn);

  map<string, string> optArgs;
  KeyvalTools::parseProcedure(optimization, optMethodModel, optArgs);

  if (optMethodModel == "D-Brent")
    // Uses Newton-Brent method or Newton-BFGS method
    optMethodModel = OptimizationTools::OPTIMIZATION_BRENT;
  else if (optMethodModel == "D-BFGS")
    optMethodModel = OptimizationTools::OPTIMIZATION_BFGS;
  else if (optMethodModel == "FullD")
    optMethodModel = OptimizationTools::OPTIMIZATION_NEWTON;
  else
    throw Exception("Unknown optimization method " + optMethodModel);

  nstep = ApplicationTools::getParameter<unsigned int>("nstep", optArgs, 1, "", true, warn + 1);


  // VERBOSITY

  verbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional, warn + 1);

  // MESSAGE

  string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
  messenger =
      (mhPath == "none") ? nullptr :
      (mhPath == "std") ? ApplicationTools::message :
      make_shared<StlOutputStream>(make_unique<ofstream>(mhPath.c_str(), ios::out));
  if (verb)
    ApplicationTools::displayResult("Message handler", mhPath);

  // PROFILE
  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
  profiler =
      (prPath == "none") ? nullptr :
      (prPath == "std") ? ApplicationTools::message :
      make_shared<StlOutputStream>(make_unique<ofstream>(prPath.c_str(), ios::out));
  if (profiler)
    profiler->setPrecision(20);
  if (verb)
    ApplicationTools::displayResult("Profiler", prPath);


  bool scaleFirst = ApplicationTools::getBooleanParameter("optimization.scale_first", params, false, suffix, suffixIsOptional, warn + 1);
  if (scaleFirst)
  {
    ApplicationTools::displayError("Sorry, optimization.scale_first not implemented yet for process.");
    exit(-1);
  }

  // Should I ignore some parameters?

  parameters = lik->getParameters();
  vector<string> parNames = parameters.getParameterNames();

  if (params.find("optimization.ignore_parameter") != params.end())
    throw Exception("optimization.ignore_parameter is deprecated, use optimization.ignore_parameters instead!");

  string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameters", params, "", suffix, suffixIsOptional, warn + 1);
  StringTokenizer st(paramListDesc, ",");
  while (st.hasMoreToken())
  {
    try
    {
      string param = st.nextToken();
      if (param == "BrLen")
      {
        vector<string> vs = lik->getBranchLengthParameters().getParameterNames();
        parameters.deleteParameters(vs);
        if (verb)
          ApplicationTools::displayResult("Parameter ignored", string("Branch lengths"));
      }
      else if (param == "Ancient")
      {
        vector<string> vs = lik->getRootFrequenciesParameters().getParameterNames();
        parameters.deleteParameters(vs);
        if (verb)
          ApplicationTools::displayResult("Parameter ignored", string("Root frequencies"));
      }
      else if (param == "Model")
      {
        vector<string> vs = lik->getSubstitutionModelParameters().getParameterNames();
        parameters.deleteParameters(vs);
        if (verb)
          ApplicationTools::displayResult("Parameter ignored", string("Model"));
      }
      else if (param == "*")
      {
        parameters.reset();
        if (verb)
          ApplicationTools::displayResult("Parameter ignored", string("All"));
      }
      else if (param.find("*") != string::npos)
      {
        vector<string> vs = ApplicationTools::matchingParameters(param, parNames);
        bool verbhere = verb;

        if (vs.size() >= 20)
        {
          if (verb)
          {
            ApplicationTools::displayResult("Number of parameters ignored", vs.size());
            ApplicationTools::displayMessage(" from " + param);
          }
          verbhere = false;
        }

        for (auto& it :  vs)
        {
          parameters.deleteParameter(it);
          if (verbhere)
            ApplicationTools::displayResult("Parameter ignored", it);
        }
      }
      else
      {
        parameters.deleteParameter(param);
        if (verb)
          ApplicationTools::displayResult("Parameter ignored", param);
      }
    }
    catch (ParameterNotFoundException& pnfe)
    {
      ApplicationTools::displayWarning("Parameter '" + pnfe.parameter() + "' not found, and so can't be ignored!");
    }
  }

  // Should I constrain some parameters?
  vector<string> parToEstNames = parameters.getParameterNames();

  if (params.find("optimization.constrain_parameter") != params.end())
    throw Exception("optimization.constrain_parameter is deprecated, use optimization.constrain_parameters instead!");

  paramListDesc = ApplicationTools::getStringParameter("optimization.constrain_parameters", params, "", suffix, suffixIsOptional, warn + 1);

  string constraint = "";
  string pc, param = "";

  StringTokenizer st2(paramListDesc, ",");
  while (st2.hasMoreToken())
  {
    try
    {
      pc = st2.nextToken();
      string::size_type index = pc.find("=");
      if (index == string::npos)
        throw Exception("PhylogeneticsApplicationTools::optimizeParamaters. Bad constrain syntax, should contain `=' symbol: " + pc);
      param = pc.substr(0, index);
      constraint = pc.substr(index + 1);
      std::shared_ptr<IntervalConstraint> ic(new IntervalConstraint(constraint));

      vector<string> parNames2;

      if (param == "BrLen")
        parNames2  = lik->getBranchLengthParameters().getParameterNames();
      else if (param == "Ancient")
        parNames2 = lik->getRootFrequenciesParameters().getParameterNames();
      else if (param == "Model")
      {
        vector<string> vs = lik->getSubstitutionModelParameters().getParameterNames();
      }
      else if (param.find("*") != string::npos)
        parNames2 = ApplicationTools::matchingParameters(param, parToEstNames);
      else
        parNames2.push_back(param);


      for (size_t i = 0; i < parNames2.size(); i++)
      {
        Parameter& par = parameters.parameter(parNames2[i]);
        if (par.hasConstraint())
        {
          par.setConstraint(std::shared_ptr<ConstraintInterface>(*ic & (*par.getConstraint())));
          if (par.getConstraint()->isEmpty())
            throw Exception("Empty interval for parameter " + parNames[i] + par.getConstraint()->getDescription());
        }
        else
          par.setConstraint(ic);

        if (verb)
          ApplicationTools::displayResult("Parameter constrained " + par.getName(), par.getConstraint()->getDescription());
      }
    }
    catch (ParameterNotFoundException& pnfe)
    {
      ApplicationTools::displayWarning("Parameter '" + pnfe.parameter() + "' not found, and so can't be constrained!");
    }
    catch (ConstraintException& pnfe)
    {
      throw Exception("Parameter '" + param + "' does not fit the constraint " + constraint);
    }
  }

  nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional, warn + 1);
  if (verb)
    ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));

  tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001, suffix, suffixIsOptional, warn + 1);
  if (verb)
    ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));

  // Backing up or restoring?
  backupFile = ApplicationTools::getAFilePath("optimization.backup.file", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
  if (backupFile != "none")
  {
    ApplicationTools::displayResult("Parameters will be backup to", backupFile);
    listener.reset(new BackupListener(backupFile));
    if (FileTools::fileExists(backupFile))
    {
      ApplicationTools::displayMessage("A backup file was found! Try to restore parameters from previous run...");
      ifstream bck(backupFile.c_str(), ios::in);
      vector<string> lines = FileTools::putStreamIntoVectorOfStrings(bck);
      double fval = TextTools::toDouble(lines[0].substr(5));
      ParameterList pl = lik->getParameters();
      for (size_t l = 1; l < lines.size(); ++l)
      {
        if (!TextTools::isEmpty(lines[l]))
        {
          StringTokenizer stp(lines[l], "=");
          if (stp.numberOfRemainingTokens() != 2)
          {
            cerr << "Corrupted backup file!!!" << endl;
            cerr << "at line " << l << ": " << lines[l] << endl;
          }
          string pname  = stp.nextToken();
          string pvalue = stp.nextToken();
          if (pl.hasParameter(pname))
          {
            size_t p = pl.whichParameterHasName(pname);
            pl.setParameter(p, AutoParameter(pl[p]));
            pl[p].setValue(TextTools::toDouble(pvalue));
          }
          else
            ApplicationTools::displayWarning("Unknown parameter in backup file : " + pname);
        }
      }
      bck.close();
      lik->setParameters(pl);
      if (convert(abs(lik->getValue() - fval)) > 0.000001)
        ApplicationTools::displayMessage("Changed likelihood from backup file.");
      ApplicationTools::displayResult("Restoring log-likelihood", -lik->getValue());
    }
  }

  string order = ApplicationTools::getStringParameter("derivatives", optArgs, "Newton", "", true, warn + 1);
  if (order == "Gradient")
  {
    optMethodDeriv = OptimizationTools::OPTIMIZATION_GRADIENT;
  }
  else if (order == "Newton")
  {
    optMethodDeriv = OptimizationTools::OPTIMIZATION_NEWTON;
  }
  else if (order == "BFGS")
  {
    optMethodDeriv = OptimizationTools::OPTIMIZATION_BFGS;
  }
  else
    throw Exception("Unknown derivatives algorithm: '" + order + "'.");

  if (verb)
    ApplicationTools::displayResult("Optimization method", optMethodModel);
  if (verb)
    ApplicationTools::displayResult("Algorithm used for derivable parameters", order);

  // See if we should reparametrize:
  reparametrization = ApplicationTools::getBooleanParameter("optimization.reparametrization", params, false, suffix, suffixIsOptional, warn + 1);
  if (verb)
    ApplicationTools::displayResult("Reparametrization", (reparametrization ? "yes" : "no"));

  // See if we should use a molecular clock constraint:
  string clock = ApplicationTools::getStringParameter("optimization.clock", params, "None", suffix, suffixIsOptional, warn + 1);
  if (clock != "None" && clock != "Global")
    throw Exception("Molecular clock option not recognized, should be one of 'Global' or 'None'.");
  useClock = (clock == "Global");
  if (verb)
    ApplicationTools::displayResult("Molecular clock", clock);
}

/******************************************************************************/

std::string OptimizationTools::OPTIMIZATION_NEWTON = "newton";
std::string OptimizationTools::OPTIMIZATION_GRADIENT = "gradient";
std::string OptimizationTools::OPTIMIZATION_BRENT = "Brent";
std::string OptimizationTools::OPTIMIZATION_BFGS = "BFGS";

/******************************************************************************/

unsigned int OptimizationTools::optimizeNumericalParameters(
    std::shared_ptr<PhyloLikelihoodInterface> lik,
    const OptimizationOptions& optopt)
{
  shared_ptr<SecondOrderDerivable> f = lik;
  ParameterList pl = optopt.parameters;

  // Shall we reparametrize the function to remove constraints?
  if (optopt.reparametrization)
  {
    f = make_shared<ReparametrizationDerivableSecondOrderWrapper>(f, optopt.parameters);

    // Reset parameters to remove constraints:
    pl = f->getParameters().createSubList(optopt.parameters.getParameterNames());
  }

  // ///////////////
  // Build optimizer:

  // Branch lengths

  auto desc = make_unique<MetaOptimizerInfos>();
  unique_ptr<MetaOptimizer> poptimizer;

  if (optopt.optMethodDeriv == OPTIMIZATION_GRADIENT)
    desc->addOptimizer("Branch length parameters", make_shared<ConjugateGradientMultiDimensions>(f), lik->getBranchLengthParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
  else if (optopt.optMethodDeriv == OPTIMIZATION_NEWTON)
    desc->addOptimizer("Branch length parameters", make_shared<PseudoNewtonOptimizer>(f), lik->getBranchLengthParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
  else if (optopt.optMethodDeriv == OPTIMIZATION_BFGS)
    desc->addOptimizer("Branch length parameters", make_shared<BfgsMultiDimensions>(f), lik->getBranchLengthParameters().getParameterNames(), 2, MetaOptimizerInfos::IT_TYPE_FULL);
  else
    throw Exception("OptimizationTools::optimizeNumericalParameters. Unknown derivative optimization method: " + optopt.optMethodDeriv);

  // Other parameters

  if (optopt.optMethodModel == OPTIMIZATION_BRENT)
  {
    ParameterList plTmp = lik->getSubstitutionModelParameters();
    plTmp.addParameters(lik->getRootFrequenciesParameters());
    ParameterList plsm = optopt.parameters.getCommonParametersWith(plTmp);

    desc->addOptimizer("Substitution model parameters", make_shared<SimpleMultiDimensions>(f), plsm.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);

    ParameterList plrd = optopt.parameters.getCommonParametersWith(lik->getRateDistributionParameters());
    desc->addOptimizer("Rate distribution parameters", make_shared<SimpleMultiDimensions>(f), plrd.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);

    poptimizer = make_unique<MetaOptimizer>(f, std::move(desc), optopt.nstep);
  }
  else if (optopt.optMethodModel == OPTIMIZATION_BFGS)
  {
    vector<string> vNameDer;
    auto fnum = make_shared<ThreePointsNumericalDerivative>(f);

    ParameterList plsm = optopt.parameters.getCommonParametersWith(lik->getSubstitutionModelParameters());
    vNameDer = plsm.getParameterNames();

    ParameterList plrd = optopt.parameters.getCommonParametersWith(lik->getRateDistributionParameters());

    vector<string> vNameDer2 = plrd.getParameterNames();

    vNameDer.insert(vNameDer.begin(), vNameDer2.begin(), vNameDer2.end());
    fnum->setParametersToDerivate(vNameDer);

    desc->addOptimizer("Rate & model distribution parameters", make_shared<BfgsMultiDimensions>(fnum), vNameDer, 1, MetaOptimizerInfos::IT_TYPE_FULL);
    poptimizer = make_unique<MetaOptimizer>(fnum, std::move(desc), optopt.nstep);
  }
  else
    throw Exception("OptimizationTools::optimizeNumericalParameters. Unknown optimization method: " + optopt.optMethodModel);

  poptimizer->setVerbose(optopt.verbose);
  poptimizer->setProfiler(optopt.profiler);
  poptimizer->setMessageHandler(optopt.messenger);
  poptimizer->setMaximumNumberOfEvaluations(optopt.nbEvalMax);
  poptimizer->getStopCondition()->setTolerance(optopt.tolerance);

  // Optimize TreeLikelihood function:
  poptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  auto nanListener = make_shared<NaNListener>(poptimizer.get(), lik.get());
  poptimizer->addOptimizationListener(nanListener);
  if (optopt.listener)
    poptimizer->addOptimizationListener(optopt.listener);

  poptimizer->init(pl);
  poptimizer->optimize();

  //  optopt.parameters.setAllParametersValues(poptimizer->getParameters());

  if (optopt.verbose > 0)
    ApplicationTools::displayMessage("\n");

  // We're done.
  unsigned int nb = poptimizer->getNumberOfEvaluations();
  return nb;
}


/************************************************************/

unsigned int OptimizationTools::optimizeNumericalParameters2(
    std::shared_ptr<PhyloLikelihoodInterface> lik,
    const OptimizationOptions& optopt)
{
  shared_ptr<SecondOrderDerivable> f = lik;
  ParameterList pl = optopt.parameters;

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
  if (optopt.reparametrization)
  {
    frep.reset(new ReparametrizationDerivableSecondOrderWrapper(f, pl));
    f = frep;

    // Reset parameters to remove constraints:
    pl = f->getParameters().createSubList(pl.getParameterNames());
  }

  // Build optimizer:
  unique_ptr<OptimizerInterface> optimizer;
  shared_ptr<AbstractNumericalDerivative> fnum;

  if (optopt.optMethodDeriv == OPTIMIZATION_GRADIENT)
  {
    fnum = make_shared<TwoPointsNumericalDerivative>(dynamic_pointer_cast<FirstOrderDerivable>(f));
    fnum->setInterval(0.0000001);
    optimizer = make_unique<ConjugateGradientMultiDimensions>(fnum);
  }
  else if (optopt.optMethodDeriv == OPTIMIZATION_NEWTON)
  {
    fnum = make_shared<ThreePointsNumericalDerivative>(f);
    fnum->setInterval(0.0001);
    optimizer = make_unique<PseudoNewtonOptimizer>(fnum);
  }
  else if (optopt.optMethodDeriv == OPTIMIZATION_BFGS)
  {
    fnum = make_shared<TwoPointsNumericalDerivative>(dynamic_pointer_cast<FirstOrderDerivable>(f));
    fnum->setInterval(0.0001);
    optimizer = make_unique<BfgsMultiDimensions>(fnum);
  }
  else
    throw Exception("OptimizationTools::optimizeNumericalParameters2. Unknown optimization method: " + optopt.optMethodDeriv);

  // Numerical derivatives:
  // Variables not derivatived in Likelihood DF but in numerical way
  ParameterList tmp = lik->getParameters();

  // if (useClock)
  //   tmp.addParameters(fclock->getHeightParameters());

  fnum->setParametersToDerivate(tmp.getParameterNames());

  optimizer->setVerbose(optopt.verbose);
  optimizer->setProfiler(optopt.profiler);
  optimizer->setMessageHandler(optopt.messenger);
  optimizer->setMaximumNumberOfEvaluations(optopt.nbEvalMax);
  optimizer->getStopCondition()->setTolerance(optopt.tolerance);

  // Optimize TreeLikelihood function:
  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  auto nanListener = make_shared<NaNListener>(optimizer.get(), lik.get());
  optimizer->addOptimizationListener(nanListener);
  if (optopt.listener)
    optimizer->addOptimizationListener(optopt.listener);

  optimizer->init(pl);
  optimizer->optimize();

  if (optopt.verbose > 0)
    ApplicationTools::displayMessage("\n");

  // We're done.
  return optimizer->getNumberOfEvaluations();
}

/************************************************************/

unsigned int OptimizationTools::optimizeNumericalParameters2(
    std::shared_ptr<SingleProcessPhyloLikelihood> lik,
    const OptimizationOptions& optopt)
{
  shared_ptr<SecondOrderDerivable> f = lik;
  ParameterList pl = optopt.parameters;
  if (optopt.reparametrization)
  {
    // Shall we reparametrize the function to remove constraints?
    f = make_shared<ReparametrizationDerivableSecondOrderWrapper>(f, optopt.parameters);

    // Reset parameters to remove constraints:
    pl = f->getParameters().createSubList(optopt.parameters.getParameterNames());
  }

  // Build optimizer:
  shared_ptr<AbstractNumericalDerivative> fnum;
  unique_ptr<OptimizerInterface> optimizer;

  if (optopt.optMethodDeriv == OPTIMIZATION_GRADIENT)
  {
    lik->likelihoodCalculationSingleProcess().setNumericalDerivateConfiguration(0.00001, NumericalDerivativeType::ThreePoints);

    fnum = make_shared<TwoPointsNumericalDerivative>(dynamic_pointer_cast<FirstOrderDerivable>(f));
    fnum->setInterval(0.0000001);
    optimizer.reset(new ConjugateGradientMultiDimensions(fnum));
  }
  else if (optopt.optMethodDeriv == OPTIMIZATION_NEWTON)
  {
    lik->likelihoodCalculationSingleProcess().setNumericalDerivateConfiguration(0.0001, NumericalDerivativeType::FivePoints);
    fnum = make_shared<ThreePointsNumericalDerivative>(f);
    fnum->setInterval(0.0001);
    optimizer.reset(new PseudoNewtonOptimizer(fnum));
  }
  else if (optopt.optMethodDeriv == OPTIMIZATION_BFGS)
  {
    lik->likelihoodCalculationSingleProcess().setNumericalDerivateConfiguration(0.0001, NumericalDerivativeType::ThreePoints);
    fnum = make_shared<TwoPointsNumericalDerivative>(dynamic_pointer_cast<FirstOrderDerivable>(f));
    fnum->setInterval(0.0001);
    optimizer.reset(new BfgsMultiDimensions(fnum));
  }
  else
    throw Exception("OptimizationTools::optimizeNumericalParameters2. Unknown optimization method: " + optopt.optMethodDeriv);


  // Variables not derivatived in Likelihood DF but in numerical way
  ParameterList tmp = f->getParameters();

  fnum->setParametersToDerivate(tmp.getParameterNames());

  optimizer->setVerbose(optopt.verbose);
  optimizer->setProfiler(optopt.profiler);
  optimizer->setMessageHandler(optopt.messenger);
  optimizer->setMaximumNumberOfEvaluations(optopt.nbEvalMax);
  optimizer->getStopCondition()->setTolerance(optopt.tolerance);

  // Optimize TreeLikelihood function:
  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  auto nanListener = make_shared<NaNListener>(optimizer.get(), lik.get());
  optimizer->addOptimizationListener(nanListener);
  if (optopt.listener)
    optimizer->addOptimizationListener(optopt.listener);

  optimizer->init(pl);
  optimizer->optimize();

  if (optopt.verbose > 0)
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
    ParameterList tmp = estimationMethod.process().getSubstitutionModelParameters(true);
    tmp.addParameters(estimationMethod.process().getRateDistributionParameters(true));
    tmp.addParameters(estimationMethod.process().getRootFrequenciesParameters(true));
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
    const std::string& param,
    OptimizationOptions& optopt)
{
  estimationMethod.resetAdditionalParameters();
  estimationMethod.setVerbose(optopt.verbose);

  if (param == DISTANCEMETHOD_PAIRWISE)
  {
    ParameterList tmp = estimationMethod.process().getSubstitutionModelParameters(true);
    tmp.addParameters(estimationMethod.process().getRateDistributionParameters(true));
    tmp.addParameters(estimationMethod.process().getRootFrequenciesParameters(true));
    tmp = tmp.getCommonParametersWith(optopt.parameters);
    estimationMethod.setAdditionalParameters(tmp);
  }

  unique_ptr<TreeTemplate<Node>> tree = nullptr;

  bool test = true;

  auto process = std::shared_ptr<SubstitutionProcessInterface>(estimationMethod.process().clone());
  auto autoProc = dynamic_pointer_cast<AutonomousSubstitutionProcessInterface>(process);
  auto procMb = dynamic_pointer_cast<SubstitutionProcessCollectionMember>(process);
  if (!autoProc && !procMb)
    throw Exception("OptimizationTools::buildDistanceTree : unknown process type. Ask developpers.");

  // Vector of successive trees, to return the best one
  std::vector<unique_ptr<TreeTemplate<Node>>> vTree;
  std::vector<double> vLik;

  size_t nstep = 0;
  while (test && nstep < 100)
  {
    // Compute matrice:
    if (optopt.verbose > 0)
      ApplicationTools::displayTask("Estimating distance matrix", true);

    estimationMethod.computeMatrix();
    auto matrix = estimationMethod.getMatrix();

    if (estimationMethod.getVerbose() > 0)
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
      return unique_ptr<TreeTemplate<Node>>(new TreeTemplate<Node>(n1));
    }

    if (optopt.verbose > 0)
      ApplicationTools::displayTask("Building tree");

    reconstructionMethod.setDistanceMatrix(*matrix);
    reconstructionMethod.computeTree();

    tree = make_unique<TreeTemplate<Node>>(reconstructionMethod.tree());

    vTree.push_back(std::move(tree));

    if (estimationMethod.getVerbose() > 0)
      ApplicationTools::displayTaskDone();

    size_t nbTree = vTree.size();
    const auto& ltree = vTree[nbTree - 1];

    if (vTree.size() > 1)
    {
      for (size_t iT = 0; iT < nbTree - 1; iT++)
      {
        const auto& pTree = vTree[iT];
        int rf = TreeTools::robinsonFouldsDistance(*pTree, *ltree);
        // if (optopt.verbose > 0)
        ApplicationTools::displayResult("Topo. distance with iteration " + TextTools::toString(iT + 1), TextTools::toString(rf));
        test &= (rf != 0);
        if (!test)
          break;
      }
    }

    if ((param != DISTANCEMETHOD_ITERATIONS) || !test)
      break; // Ends here.

    // Now, re-estimate parameters:
    Context context;

    auto phyloTree  = make_shared<ParametrizablePhyloTree>(*PhyloTreeTools::buildFromTreeTemplate(*ltree));
    if (autoProc)
      autoProc->setPhyloTree(*phyloTree);
    else if (procMb)
    {
      auto& coll = procMb->collection();
      size_t maxTNb = procMb->getTreeNumber();
      coll.replaceTree(phyloTree, maxTNb);
    }

    auto lik     = make_shared<LikelihoodCalculationSingleProcess>(context, estimationMethod.getData(), process);
    auto tl      = make_shared<SingleProcessPhyloLikelihood>(context, lik);

    vLik.push_back(tl->getValue());

    // hide opt verbose
    optopt.verbose = estimationMethod.getVerbose() > 0 ? (unsigned int)(estimationMethod.getVerbose() - 1) : 0;

    optimizeNumericalParameters(tl, optopt);
    process->matchParametersValues(tl->getParameters());

    estimationMethod.matchParametersValues(process->getParameters());

    auto trtemp = std::make_shared<ParametrizablePhyloTree>(*tl->tree());
    const PhyloTree trt2(*trtemp);
    tree.reset(TreeTemplateTools::buildFromPhyloTree(trt2).release());

    if (optopt.verbose > 0)
    {
      auto tmp = process->getSubstitutionModelParameters(true);
      for (unsigned int i = 0; i < tmp.size(); ++i)
      {
        ApplicationTools::displayResult(tmp[i].getName(), TextTools::toString(tmp[i].getValue()));
      }
      tmp = process->getRootFrequenciesParameters(true);
      for (unsigned int i = 0; i < tmp.size(); ++i)
      {
        ApplicationTools::displayResult(tmp[i].getName(), TextTools::toString(tmp[i].getValue()));
      }
      tmp = process->getRateDistributionParameters(true);
      for (unsigned int i = 0; i < tmp.size(); ++i)
      {
        ApplicationTools::displayResult(tmp[i].getName(), TextTools::toString(tmp[i].getValue()));
      }
    }
    nstep++;
  }

  size_t posM = static_cast<size_t>(std::distance(vLik.begin(), std::min_element(vLik.begin(), vLik.end())));

  return std::move(vTree[posM]);
}

/******************************************************************************/
