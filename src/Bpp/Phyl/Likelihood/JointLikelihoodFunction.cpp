#include "JointLikelihoodFunction.h"

// for bpp-core
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Text/KeyvalTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTools.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/MixedSubstitutionModel.h>
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/Codon/RELAX.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Mapping/StochasticMapping.h>

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

// From the STL:
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <vector>
#include <math.h>       /* pow exp */
#include <string>
#include <algorithm>    // std::min
#include <limits>       // std::numeric_limits (for -inf setting)
#include <map>

using namespace bpp;
using namespace std;

/******************************************************************************/

JointLikelihoodFunction::JointLikelihoodFunction(BppApplication* bppml, const Tree* tree, const VectorSiteContainer* characterData, TransitionModel* characterModel, const VectorSiteContainer* sequenceData, MixedSubstitutionModelSet* sequenceModel, DiscreteDistribution* rDist, bool debug):
AbstractParametrizable(""), 
hypothesis_(null),
bppml_(bppml),
characterTreeLikelihood_(),
sequenceTreeLikelihood_(),
jointParameters_(),
previousParametersValues_(),
stocMapping_(),
optimizationScope_(none),
logl_(0),
characterChanged_(false),
sequenceChanged_(false),
previousK_(1),
origTreeLength_(0),
debugDir_(),
debug_(false),
cycleNum_(0)
{
    // get the original total tree length
    const TreeTemplate<Node>* ttree = dynamic_cast<const TreeTemplate<Node>*>(tree);
    
    vector <const Node*> nodes = ttree->getNodes();
    for (size_t b=0; b<nodes.size(); ++b)
    {
        if (nodes[b]->getId() != tree->getRootId())
        
        {
            origTreeLength_ = origTreeLength_ + nodes[b]->getDistanceToFather();
        }
    }
    
    // create the character likelihood function as data member based on the character model, character data and tree
    RHomogeneousTreeLikelihood* characterTreeLikelihood = new RHomogeneousTreeLikelihood(*tree, *characterData, characterModel, rDist, false);
    characterTreeLikelihood->initialize();
    characterTreeLikelihood_ = characterTreeLikelihood;

    // create the sequence tree likelihood of the null model as initial data memeber based on the sequence model, sequence data and tree
    RNonHomogeneousMixedTreeLikelihood* sequenceTreeLikelihood = new RNonHomogeneousMixedTreeLikelihood(*tree, *sequenceData, sequenceModel, rDist, true, true);
    sequenceTreeLikelihood->initialize();
    sequenceTreeLikelihood_ = sequenceTreeLikelihood;

    // create the stochadtic mapping data member based on the character likelihood function data member
    size_t numOfMappings = static_cast<size_t>(ApplicationTools::getIntParameter("character.num_of_mappings", bppml_->getParams(), 1000));
    stocMapping_ = new StochasticMapping(dynamic_cast<TreeLikelihood*>(characterTreeLikelihood_), numOfMappings);

    // maintain a list of character model and sequence model parameters to return from getParameters()
    jointParameters_.addParameters(characterModel->getParameters());
    jointParameters_.addParameters(sequenceModel->getParameters());

    // set the previuos parameters list to the initial ones
    for (size_t i=0; i<jointParameters_.size(); ++i)
    {
        previousParametersValues_[jointParameters_[i].getName()] = jointParameters_[i].getValue();
    }

    // get the expected histories dir from bppml_
    debugDir_ = ApplicationTools::getStringParameter("output.debug.dir", bppml_->getParams(), "", "", true, 1);
    debug_ = debug;
}

/******************************************************************************/

JointLikelihoodFunction::~JointLikelihoodFunction()
{
    // delete the stochastic mapping
    if (stocMapping_) delete stocMapping_;

    // delete the likelihood function
    if (characterTreeLikelihood_) delete characterTreeLikelihood_;

    // delete the sequence likelihood function
    if (sequenceTreeLikelihood_) delete sequenceTreeLikelihood_;

    // deletion of the parameters created via the constrcutor is done by the the deletion of the inheriting class AbstractParametrizable
}

/******************************************************************************/

void JointLikelihoodFunction::updatePreviousParametersValues()
{
    string parName;
    double parValue;
    for (size_t i=0; i < jointParameters_.size(); ++i)
    {
        parName = jointParameters_[i].getName();
        parValue = jointParameters_[i].getValue();
        // check if any re-computationo flag needs to be updated
        if (previousParametersValues_[parName] != parValue)
        {
            if (parName.find("TwoParameterBinary") != std::string::npos)
            {
                characterChanged_ = true;
            }
            else if (parName.find("RELAX") != std::string::npos)
            {
                sequenceChanged_ = true;
            }
        }
        previousParametersValues_[parName] = parValue;
    }
}

/******************************************************************************/

void JointLikelihoodFunction::fireParameterChanged(const ParameterList& pl) 
{  
    /* if a character parameter was changed -> call computeAlternativeJointLikelihood without sequence model optimization
       else, the changed paramet der must be the boolean distating that the sequence model should be optimized */
    switch(hypothesis_)
    {
        case null: 
            computeNullJointLikelihood(); 
            break;
        case alternative: 
            computeAlternativeJointLikelihood(); 
            break;
        default: 
            throw Exception("Error! illegal hypothesis setting");
    }

    // update the log likelihood of the model
    logl_ = characterTreeLikelihood_->getValue() + sequenceTreeLikelihood_->getValue();

    // update the previous parameters list
    updatePreviousParametersValues();
}

/******************************************************************************/

void JointLikelihoodFunction::updateStatesInNodesNames(Tree* mapping)
{
    string label = "state";
    vector<Node*> nodes = (dynamic_cast<TreeTemplate<Node>*>(mapping))->getNodes(); 
    for (int i=0; i < static_cast<int>(nodes.size()); i++) 
    {
        string name = nodes[i]->getName();
        int state = stocMapping_->getNodeState(nodes[i]);
        nodes[i]->setName(name + "{" + TextTools::toString(state) + "}");
    }
}

/******************************************************************************/

void JointLikelihoodFunction::setPartitionByHistory(Tree* history)
{
    sequenceTreeLikelihood_->getSubstitutionModelSet()->resetModelToNodeIds();
    vector<const Node*> nodes = (dynamic_cast<const TreeTemplate<Node>*>(history))->getNodes();
    for (size_t i=0; i<nodes.size(); ++i)
    {
      int nodeId = nodes[i]->getId();
      if (nodes[i]->hasFather())
      {
        int nodeState = stocMapping_->getNodeState(nodes[i]);
        if (nodeState == 0)
        {
            sequenceTreeLikelihood_->getSubstitutionModelSet()->setNodeToModel(0,nodeId);
        }    
        else
        {
            sequenceTreeLikelihood_->getSubstitutionModelSet()->setNodeToModel(1,nodeId);
        }
      }
    }
}

/******************************************************************************/

void JointLikelihoodFunction::updatesequenceTreeLikelihood(const Tree* history)
{
    // extract the input for the next SequenceTreeLikelihood from the previouts one
    const VectorSiteContainer* sequenceData = dynamic_cast<const VectorSiteContainer*>(sequenceTreeLikelihood_->getData()->clone());
    MixedSubstitutionModelSet* sequenceModel = dynamic_cast<MixedSubstitutionModelSet*>(sequenceTreeLikelihood_->getSubstitutionModelSet());
    DiscreteDistribution* rDist = RASTools::getPosteriorRateDistribution(*sequenceTreeLikelihood_);

    // create the new sequenceTreeLikelihood
    RNonHomogeneousMixedTreeLikelihood* sequenceTreeLikelihood = new RNonHomogeneousMixedTreeLikelihood(*history, *sequenceData, sequenceModel, rDist, true, true);
    sequenceTreeLikelihood->initialize();

    // delete the previous SequenceTreeLikelihood instance
    if (sequenceTreeLikelihood_) delete sequenceTreeLikelihood_;

    // assign a new sequenceTreeLikelihood instnace
    sequenceTreeLikelihood_ = sequenceTreeLikelihood;
}

/******************************************************************************/

double JointLikelihoodFunction::getSequenceScalingFactor(bool verbose)
{
    const Tree& tree = sequenceTreeLikelihood_->getTree();
    vector <const Node*> nodes = (dynamic_cast<const TreeTemplate<Node>&>(tree)).getNodes();
    double treeSize = 0;
    for (size_t b=0; b<nodes.size(); ++b)
    {
        if (nodes[b]->getId() != tree.getRootId())

        {
            treeSize = treeSize + nodes[b]->getDistanceToFather();
        }
    }
    double scalingFactor = treeSize / origTreeLength_;
    if (verbose)
    {
        cout << "The tree has been scaled by a sequence scaling factor of: " << scalingFactor << endl;
    }
    return scalingFactor;
}

/******************************************************************************/

map<string,double> JointLikelihoodFunction::getSequenceModelParameters(bool verbose)
{
  MixedSubstitutionModelSet* sequenceModel = dynamic_cast<MixedSubstitutionModelSet*>(sequenceTreeLikelihood_->getSubstitutionModelSet());
  ParameterList parameters;
  map<string,double> parametersValues;
  for (size_t m = 0; m < sequenceModel->getNumberOfModels(); ++m) {
    if (verbose)
      ApplicationTools::displayMessage("\nmodel " + TextTools::toString(m+1) + "\n");
    TransitionModel* model = sequenceModel->getModel(m);
    parameters = model->getParameters();
    for (size_t i = 0; i < parameters.size(); i++)
    {
      if (verbose)
        ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
      if (((parameters[i].getName().compare("RELAX.k") == 0) & (m == 1)) | ((parameters[i].getName().compare("RELAX.k") != 0) & (m == 0)))
        parametersValues[parameters[i].getName() + "_" + TextTools::toString(m+1)] = parameters[i].getValue();
    }
  }
  if (verbose)
  {
    parameters = sequenceTreeLikelihood_->getRateDistributionParameters();
    for (size_t i = 0; i < parameters.size(); i++)
    {
      ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }
  }

  // get the scaling factor
  parametersValues["sequenceScalingFactor"] = getSequenceScalingFactor(false);
  ApplicationTools::displayResult("Sequence scaling factor",  TextTools::toString(parametersValues["sequenceScalingFactor"]));

  double characterLogl = -1.0 * characterTreeLikelihood_->getValue();
  double sequenceLogl = -1.0 * sequenceTreeLikelihood_->getValue();
  double overallLogl = characterLogl + sequenceLogl;
  ApplicationTools::displayResult("\nCharacter Log likelihood", TextTools::toString(characterLogl));
  ApplicationTools::displayResult("Sequence Log likelihood", TextTools::toString(sequenceLogl));
  ApplicationTools::displayResult("Overall Log likelihood", TextTools::toString(overallLogl));
  cout << "\n" << endl;
  parametersValues["logl"] = overallLogl;

  return parametersValues;
}

/******************************************************************************/

void JointLikelihoodFunction::scaleSequenceTree(double factor)
{
    // get a new tree and scale it
    Tree* newTree = sequenceTreeLikelihood_->getTree().clone();
    (dynamic_cast<TreeTemplate<Node>*>(newTree))->scaleTree(factor);
    // switch the new tree with the ols tree in the sequence likelihood function
    updatesequenceTreeLikelihood(newTree);
}

/******************************************************************************/

map<string,double> JointLikelihoodFunction::getModelParameters(bool verbose)
{
    map<string,double> modelParameters = getSequenceModelParameters(false);
    
    // report character model parameters and likelihood
    ParameterList parameters;
    TransitionModel* characterModel = characterTreeLikelihood_->getModel();
    parameters = characterModel->getParameters();
	cout << "\n" << endl;
    for (size_t i = 0; i < parameters.size(); i++)
    {
        if (verbose)
        {
            ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
        }
        modelParameters[parameters[i].getName()] = parameters[i].getValue();
    }
    double charLogL = characterTreeLikelihood_->getValue();
    if (verbose)
    {
        ApplicationTools::displayResult("Character Log likelihood", TextTools::toString(-charLogL, 15));
    }

    // report sequence model parameters
    if (verbose)
    {
        getSequenceScalingFactor();
        MixedSubstitutionModelSet* sequenceModel = dynamic_cast<MixedSubstitutionModelSet*>(sequenceTreeLikelihood_->getSubstitutionModelSet());
        for (size_t m = 0; m < sequenceModel->getNumberOfModels(); ++m) 
        {
            ApplicationTools::displayMessage("\nmodel " + TextTools::toString(m+1) + "\n");
            TransitionModel* model = sequenceModel->getModel(m);
            parameters = model->getParameters();
            for (size_t j = 0; j < parameters.size(); j++)
            {
                ApplicationTools::displayResult(parameters[j].getName(), TextTools::toString(parameters[j].getValue()));
            }
        }
    }
    double seqLogL = sequenceTreeLikelihood_->getValue();
    if (verbose)
    {
        ApplicationTools::displayResult("Sequence Log likelihood", TextTools::toString(-seqLogL, 15));
    }

    // report joint likelihood
    double overallLogL = charLogL + seqLogL; // = log(charLikelihood * seqLikelihood)
    logl_ = overallLogL;
    ApplicationTools::displayResult("Overall Log likelihood", TextTools::toString(-overallLogL, 15));
    cout << "\n" << endl;
    return modelParameters;
}


/******************************************************************************/

void JointLikelihoodFunction::optimizeOmegaParameters()
{
    double kValue = sequenceTreeLikelihood_->getParameterValue("RELAX.k_2");
    double upperFreedom = pow(999, max(kValue, 1.0));
    double lowerFreedom = pow(0.001, max(kValue, 1.0));
    const IntervalConstraint* praramBounds;
    // increase the upper bound of omega2 according to the upper freedom
    sequenceTreeLikelihood_->setParameterBounds("RELAX.omega2_1", 1, upperFreedom);    
    // divide the lower freedom equally between p and omega1 and decrease their lower bounds accordingly
    praramBounds = dynamic_cast<const IntervalConstraint*>(sequenceTreeLikelihood_->getParameter("RELAX.p_1").getConstraint());
    double pLb = praramBounds->getLowerBound();
    double pValue = sequenceTreeLikelihood_->getParameterValue("RELAX.p_1");
    double pNeedOfFreedom = 1 - (pValue - pLb) / (1 - pLb); // (1-pLb) is the range of the parameter, the this value is a measurement of how close the parameter value is to its lower bound (the closer it is, the larger its need for freedom, so the computed measurement is larger)
    praramBounds = dynamic_cast<const IntervalConstraint*>(sequenceTreeLikelihood_->getParameter("RELAX.omega1_1").getConstraint());
    double omega1Lb = praramBounds->getLowerBound();
    double omega1Value = sequenceTreeLikelihood_->getParameterValue("RELAX.omega1_1");
    double omega1NeedOfFreedom = 1 - (omega1Value - omega1Lb) / (1 - omega1Lb);
    pLb = pow(lowerFreedom, (pNeedOfFreedom / (pNeedOfFreedom + omega1NeedOfFreedom)));
    omega1Lb = pow(lowerFreedom, (omega1NeedOfFreedom / (pNeedOfFreedom + omega1NeedOfFreedom)));
    sequenceTreeLikelihood_->setParameterBounds("RELAX.p_1", pLb, 1);
	sequenceTreeLikelihood_->setParameterBounds("RELAX.omega1_1", omega1Lb, 1);
    // create an optimizer for p, omega1 and omega2 and then optimize the three parameters all the others, including k
    string paramsToIgnore = "BrLen,RELAX.kappa_1,RELAX.theta1_1,RELAX.theta2_1,RELAX.k_1,RELAX.k_2,RELAX.1_Full.theta_1,RELAX.1_Full.theta1_1,RELAX.1_Full.theta2_1,RELAX.2_Full.theta_1,RELAX.2_Full.theta1_1,RELAX.2_Full.theta2_1,RELAX.3_Full.theta_1,RELAX.3_Full.theta1_1,RELAX.3_Full.theta2_1";
    bppml_->getParam("optimization.ignore_parameters") = paramsToIgnore;
    try
    {
        PhylogeneticsApplicationTools::optimizeParameters(sequenceTreeLikelihood_, sequenceTreeLikelihood_->getParameters(), bppml_->getParams());
    }
    catch (Exception& e) // if needed, adjust the boundaries of the parameter
    {
        cout << "bug in optimization of omega parameters" << endl;
    }    
}

/******************************************************************************/

void JointLikelihoodFunction::optimizeSequenceOneDimension(const std::string& paramName)
{
    // compute the required lower and upper bound of the examined parameter and update its boundaries accordingly 
    const IntervalConstraint* praramBounds;
    praramBounds = dynamic_cast<const IntervalConstraint*>(sequenceTreeLikelihood_->getParameter(paramName).getConstraint());
    double lb = praramBounds->getLowerBound();
    double ub = praramBounds->getUpperBound();
    if (paramName.compare("RELAX.p_1") == 0)
    {
        lb = pow(0.001, 1/ max(sequenceTreeLikelihood_->getParameterValue("RELAX.k_2"), 1.0)) / sequenceTreeLikelihood_->getParameterValue("RELAX.omega1_1");
        sequenceTreeLikelihood_->setParameterBounds("RELAX.p_1", lb, ub);       
    }
    else if (paramName.compare("RELAX.omega1_1") == 0)
    {
        lb = pow(0.001, 1 / max(sequenceTreeLikelihood_->getParameterValue("RELAX.k_2"), 1.0)) / sequenceTreeLikelihood_->getParameterValue("RELAX.p_1");
        sequenceTreeLikelihood_->setParameterBounds("RELAX.omega1_1", lb, ub);
    }
    else if (paramName.compare("RELAX.omega1_2") == 0)
    {
        ub = pow(999, max(sequenceTreeLikelihood_->getParameterValue("RELAX.k_2"), 1.0));
        sequenceTreeLikelihood_->setParameterBounds("RELAX.omega2_1", lb, ub);
    }
    else if (paramName.compare("RELAX.k_2") == 0)
    {
        // get the values of the smallest_bg_omega (omega0) and largest_bg_omega (omega2) in the background model (i.e., model 1)
        double smallest_bg_omega = sequenceTreeLikelihood_->getParameterValue("RELAX.p_1") * sequenceTreeLikelihood_->getParameterValue("RELAX.omega1_1");
        double largest_bg_omega = sequenceTreeLikelihood_->getParameterValue("RELAX.omega2_1");
        ub = min((log(0.001)/log(smallest_bg_omega)), (log(999)/log(largest_bg_omega)));
        sequenceTreeLikelihood_->setParameterBounds("RELAX.k_2", lb, ub);
    }
    praramBounds = dynamic_cast<const IntervalConstraint*>(sequenceTreeLikelihood_->getParameter(paramName).getConstraint());

    // optimize the sequence likelihood function with respect to the given parameter
    BrentOneDimension* optimizer = new BrentOneDimension(sequenceTreeLikelihood_);
    optimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);
    optimizer->getStopCondition()->setTolerance(0.0001); // set the tolerance to be slighly less strict to account for the instability of the joint likelihood function
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setVerbose(1);
    ParameterList param; 
    param.addParameter(sequenceTreeLikelihood_->getParameter(paramName));
    praramBounds = dynamic_cast<const IntervalConstraint*>(sequenceTreeLikelihood_->getParameter(paramName).getConstraint());
    optimizer->setInitialInterval(praramBounds->getLowerBound(), praramBounds->getUpperBound()); // search within stricter bounds that the actual ones of pi0 to avoid failute of stochasitc mapping
    try
    {
        optimizer->init(param);
        optimizer->optimize();
    }
    catch (Exception& e) // if needed, adjust the boundaries of the parameter
    {
        cout << "bug in optimization of A single parameter" << endl;
    }

    // now constrain the optimized parameter according to its obtained value
    if (paramName.compare("RELAX.k_2") == 0)
    {
        lb = 0;
        ub = sequenceTreeLikelihood_->getParameterValue(paramName);
        sequenceTreeLikelihood_->setParameterBounds(paramName, lb, ub);
        praramBounds = dynamic_cast<const IntervalConstraint*>(sequenceTreeLikelihood_->getParameter("RELAX.k_2").getConstraint());
	    cout << "for parameter k, bounds are: (" << praramBounds->getLowerBound() << ", " << praramBounds->getUpperBound() << ")" << endl;
    }
    else if (paramName.compare("RELAX.omega2_1") == 0)
    {
        lb = 1;
        ub = sequenceTreeLikelihood_->getParameterValue(paramName);
        sequenceTreeLikelihood_->setParameterBounds(paramName, lb, ub);  
    }
    else
    {
        lb = sequenceTreeLikelihood_->getParameterValue(paramName);
        ub = 1;
        sequenceTreeLikelihood_->setParameterBounds(paramName, lb, ub);
    }
}


/******************************************************************************/

bool JointLikelihoodFunction::checkIfParameterRechedBound(const Parameter& parameter)
{
    const IntervalConstraint* praramBounds = dynamic_cast<const IntervalConstraint*>(parameter.getConstraint());
    double lb = praramBounds->getLowerBound();
    double ub = praramBounds->getUpperBound();
    double parValue = parameter.getValue();
    if (parValue == lb || parValue == ub)
    {
        return true;
    }
    return false;
}

/******************************************************************************/

void JointLikelihoodFunction::optimizeSequenceWithDynamicBounds(uint method)
{
    vector<std::string> parametersToOptimize;
    parametersToOptimize.push_back("RELAX.k_2");
    parametersToOptimize.push_back("RELAX.p_1");
    parametersToOptimize.push_back("RELAX.omega1_1");
    parametersToOptimize.push_back("RELAX.omega2_1");

    vector<int> parametersOptimizationOrder;
    parametersOptimizationOrder.push_back(0);
    parametersOptimizationOrder.push_back(1);
    parametersOptimizationOrder.push_back(2);
    parametersOptimizationOrder.push_back(3);

    double tolerance = 0.01;
    double prevLogl = std::numeric_limits<int>::min();
    double currLogl = -1 * sequenceTreeLikelihood_->getValue();
    ParameterList parameters = sequenceTreeLikelihood_->getSubstitutionModelParameters();
    map<string,double> bestParametersValues;
    for (size_t i=0; i<parameters.size(); ++i)
    {
        bestParametersValues[static_cast<string>(parameters[i].getName())] = parameters[i].getValue();
    }
    string logl_key = "logl"; 
    bestParametersValues[logl_key] = currLogl;

    // option 1: repeat iteratively: first optimize k, then fix k and expand the boundaries of the other parameters and optimize all of them at once
    if (method == 0)
    {
        while ((currLogl-prevLogl) > -1*tolerance)
        {
            optimizeSequenceOneDimension("RELAX.k_2");
            optimizeOmegaParameters();

            // update the log likelihood values and the values of the parameters from the current visited point, if the likelihood has improved
            prevLogl = currLogl;
            currLogl = -1* (sequenceTreeLikelihood_->getValue());
            if (currLogl > bestParametersValues["logl"])
            {
                for (size_t i=0; i<parameters.size(); ++i)
                {
                    bestParametersValues[parameters[i].getName()] = parameters[i].getValue();
                }  
            }
        }
    }

    // option 2: repeat iteratively: optimize each parameter seperately while fixing the others and expanding its bounds on expense of the fixed parameters
    else
    {
        while ((currLogl-prevLogl) > -1*tolerance)
        {
            // optimize the sequence likelihood funxtion with repsect to each parameter seperately
            int parameterIndex;
            std::string parameterName;
            for (size_t p=0; p<parametersOptimizationOrder.size(); ++p)
            {
                parameterIndex = parametersOptimizationOrder[p];
                parameterName = parametersToOptimize[parameterIndex];
                optimizeSequenceOneDimension(parameterName);
            }
            // update the log likelihood values and the values of the parameters from the current visited point, if the likelihood has improved
            prevLogl = currLogl;
            currLogl = -1* (sequenceTreeLikelihood_->getValue());
            if (currLogl > bestParametersValues["logl"])
            {
                for (size_t i=0; i<parameters.size(); ++i)
                {
                    bestParametersValues[parameters[i].getName()] = parameters[i].getValue();
                }  
            }

            // change the order of parameters to be optimized
            int firstOptimized = parametersOptimizationOrder[0];
            parametersOptimizationOrder[0] = parametersOptimizationOrder[1];
            parametersOptimizationOrder[1] = parametersOptimizationOrder[2];
            parametersOptimizationOrder[2] = parametersOptimizationOrder[3];
            parametersOptimizationOrder[3] = firstOptimized;
        }
    }

    // if the last visited point is not the best, update the parameters according to the best visited point
    if (bestParametersValues[logl_key] > currLogl)
    {
        for (size_t i=0; i<parameters.size(); ++i)
        {
            try
            {
                parameters[i].setValue(bestParametersValues[parameters[i].getName()]);
            }
            catch (Exception& e) // if needed, adjust the boundaries of the parameter
            {
                double lb, ub;
                // set the new boundaries according to the parameter of interest
                if ((parameters[i].getName()).compare("RELAX.k_2") == 0)
                {
                    lb = 0;
                    ub = bestParametersValues[parameters[i].getName()];
                }
                else if ((parameters[i].getName()).compare("RELAX.omega2_1") == 0)
                {
                    lb = 1;
                    ub = bestParametersValues[parameters[i].getName()];
                }
                else
                {
                    lb = bestParametersValues[parameters[i].getName()];
                    ub = 1;
                }
                sequenceTreeLikelihood_->setParameterBounds(parameters[i].getName(), lb, ub);
                parameters[i].setValue(bestParametersValues[parameters[i].getName()]); 
            }
        }  
    }
}

/******************************************************************************/

void JointLikelihoodFunction::optimizeSequenceModel()
{
    MixedSubstitutionModelSet* sequenceModel = dynamic_cast<MixedSubstitutionModelSet*>(sequenceTreeLikelihood_->getSubstitutionModelSet());
    string paramsToIgnore;
    // extract the user initial value of k for potential later use
    string FGModelInitialValues = ApplicationTools::getStringParameter("model2", bppml_->getParams(), "RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1)", "", true, true);
    string modelName = "RELAX";
    map<string, string> args;
    KeyvalTools::parseProcedure(FGModelInitialValues, modelName, args);
    double userInitialValue = TextTools::toDouble(args["k"]);
    switch(hypothesis_)
    {
        case null: 
            // set the value of RELAX.k_2 to 1 as well and then ignore
            sequenceTreeLikelihood_->setParameterValue("RELAX.k_2", 1);
            // now add RELAX.k_2 to the set of parameters to ignore
            paramsToIgnore = "BrLen,RELAX.k_1,RELAX.k_2,";
            break;      // k = 1 both in model1 (BG) and model2 (FG)
        case alternative: 
            // reset the initial value of RELAX.k_2 to the initial value form the parameters file (will run over the setting as 1 in case of preceding null model optimization)
            sequenceModel->getModel(1)->setParameterValue("k", userInitialValue);
            paramsToIgnore = "BrLen,RELAX.k_1,RELAX.1_Full.theta_1,RELAX.1_Full.theta1_1,RELAX.1_Full.theta2_1,RELAX.2_Full.theta_1,RELAX.2_Full.theta1_1,RELAX.2_Full.theta2_1,RELAX.3_Full.theta_1,RELAX.3_Full.theta1_1,RELAX.3_Full.theta2_1"; // ignore frequency parameters to reduce optimization duration - results in one unit of ll reduction in optimality and 1 minutre reduction in duration
            break;      // k = 1 only in model1 (BG)
        default:
            throw Exception("Error! illegal hypothesis setting");
    }
    bppml_->getParam("optimization.ignore_parameters") = paramsToIgnore;
    // first optimize the scaling parameter: 0 = don't scale, 1 - scale only before optimization, 2 - scale only after optimization, 3 - scale before and after optimization
    int scaleTree = ApplicationTools::getIntParameter("optimization.scale.tree", bppml_->getParams(), 1);
    if (scaleTree == 1)
    {
        OptimizationTools::optimizeTreeScale(sequenceTreeLikelihood_, 0.000001, 1000000, ApplicationTools::message.get(), ApplicationTools::message.get(), 0);
        getSequenceScalingFactor(); // debug
    }

    // now optimize the rest of parameters
    PhylogeneticsApplicationTools::optimizeParameters(sequenceTreeLikelihood_, sequenceTreeLikelihood_->getParameters(), bppml_->getParams());
    

    // if p / omega1 / omega2 / k recieved MLEs at theri boundaries, optimize each parameter seperately while allowing dynamic boundaries, until convergence is obtained
    // bool optimizeWithDynamicBounds = false;
    // ParameterList parameters = sequenceTreeLikelihood_->getSubstitutionModelParameters();
    // for (size_t i = 0; i < parameters.size(); i++)
    // {
    //     if ((parameters[i].getName().compare("RELAX.p_1") == 0) || (parameters[i].getName().compare("RELAX.omega1_1") == 0) || (parameters[i].getName().compare("RELAX.omega2_1") == 0) || (parameters[i].getName().compare("RELAX.k_2") == 0))
    //     {
    //         if (checkIfParameterRechedBound(parameters[i]))
    //         {
    //             optimizeWithDynamicBounds = true;
    //             break;
    //         }
    //     }

    // }
    // 
    // if (optimizeWithDynamicBounds)
    // {
    //     optimizeSequenceWithDynamicBounds();
    // }
    
    // update the log likelihood of the joint model
    double charLogL = characterTreeLikelihood_->getValue();
    double seqLogL = sequenceTreeLikelihood_->getValue();
    double overallLogL = charLogL + seqLogL; // = log(charLikelihood * seqLikelihood)
    logl_ = overallLogL;
}

/******************************************************************************/

void JointLikelihoodFunction::computeAlternativeJointLikelihood()
{  
    if (characterChanged_)
    {
        /* compute likelihood of charcter model */
        characterTreeLikelihood_->computeTreeLikelihood();
        
        /* approximate the expected character history based on numOfMappings sampled stochastic mappings */
        bool useAnalytic =  static_cast<bool>(ApplicationTools::getIntParameter("character.use_analytic_mapping", bppml_->getParams(), 0));
        vector<Tree*> mappings;
        if (debug_ & !useAnalytic)
        {
            cout << "Generating stochastic mappings\n" << endl;
            bppml_->startTimer();
        }
        if (!useAnalytic)
        {
            stocMapping_->generateStochasticMapping(mappings);
        }
        if (debug_ & !useAnalytic)
        {
            cout << "Completed generation of stochastic mappings\n" << endl;
            bppml_->done();
        }
        if (debug_)
        {
            cout << "Generating expected history based on the stochastic mappings\n" << endl;
            bppml_->startTimer();
        }
        Tree* expectedHistory;
        if (useAnalytic)
        {
            expectedHistory = stocMapping_->generateAnalyticExpectedMapping();
        }
        else
        {
            expectedHistory = stocMapping_->generateExpectedMapping(mappings);
        }
        if (debug_)
        {
            cout << "Completed generation of expected history\n" << endl;
            bppml_->done();
        }

        /* debug start */
        // likelihood based analysis
        if (debug_ & !useAnalytic)
        {
            /* compute the log likelihood based on mutiple mappings (exhaustive approximation) */
            cout << "Computing log likelihood based on the exhaustive computation\n" << endl;
            bppml_->startTimer();
            double charLogl = -1*(characterTreeLikelihood_->getValue());
            // for each mapping, set the parittion according ot it and define it as a tree of the cloned sequence likelihood function
            ApplicationTools::displayResult("Character model log likelihood: ", TextTools::toString(charLogl, 15));
            cout << "Computing sequence log likelihoods given the different mappings\n" << endl;
            for (size_t h=0; h<mappings.size(); ++h)
            {
                    setPartitionByHistory(mappings[h]); // induce a partition of the tree based on the epxected character history
                    updatesequenceTreeLikelihood(mappings[h]); // compute the likelihood given the mapping
                    cout << TextTools::toString(charLogl + -1*sequenceTreeLikelihood_->getValue(), 15) << endl;
            
            }
            // the computation in exhaustive approximation will be done via python 
            bppml_->done();
        }
        if (debug_)
        {
            /* compute the likelihood based on the expected history approxmation */
            cout << "Computing log likelihood based on the expected history approximation\n" << endl;
            bppml_->startTimer();
            setPartitionByHistory(expectedHistory);         // induce a partition of the tree based on the epxected character history
            updatesequenceTreeLikelihood(expectedHistory);  // update the tree of the sequence likelihood function to be the expected history and then compute the likelihood
            double charLogl = -1*(characterTreeLikelihood_->getValue());
            double sequenceLogl = -1*(sequenceTreeLikelihood_->getValue());
            double overallLogl = charLogl + sequenceLogl;
            ApplicationTools::displayResult("log likelihood of TraitRELAX model given the expected history", TextTools::toString(overallLogl, 15));
            bppml_->done();
        }
        if (debug_ & !useAnalytic)
        {
            /* distance based analysis will be done via python */
            string treeStr;
            string filepath = debugDir_ + TextTools::toString(mappings.size()) + "_mappings_in_nwk.txt";
            ofstream file (filepath);
            for (size_t h=0; h<mappings.size(); ++h)
            {
                // write the mappings into a file
                updateStatesInNodesNames(mappings[h]);
                // write newick string to file
                treeStr = TreeTools::treeToParenthesis(*mappings[h]);
                file << treeStr << "\n";
            }
        }
        if (debug_)
        {
            // write the expected history to a file
            size_t numOfMappings = static_cast<size_t>(ApplicationTools::getIntParameter("character.num_of_mappings", bppml_->getParams(), 1000));
            string method;
            if (useAnalytic)
            {
                method = "analytic";
            }
            else
            {
                method = "sampling_based";
            }
            Tree* clonedExpectedHistory = expectedHistory->clone();// clone the expected history and add states to its names
            updateStatesInNodesNames(clonedExpectedHistory);
            // set the name of the expected history path according to the current cycle
            string treeStr = TreeTools::treeToParenthesis(*clonedExpectedHistory);
            string filepath = debugDir_ + TextTools::toString(numOfMappings) + "_" + method + "_expected_mapping.nwk";
            ofstream file (filepath);
            
            file << treeStr << "\n";
            if (useAnalytic)
            {
                cout << "tree:\n" << treeStr << endl; 
            }
            if (clonedExpectedHistory) delete clonedExpectedHistory; // delete the clone
            file.close();
        }
        /* debug - end */

        /* compute the likelihood of the sequence model */
        setPartitionByHistory(expectedHistory);         // induce a partition of the tree based on the epxected character history
        updatesequenceTreeLikelihood(expectedHistory);  // update the tree of the sequence likelihood function to be the expected history and then compute the likelihood

        /* if needed, optimize the model */
        switch(optimizationScope_)
        {
            case none: 
                break;
            case onlyCharacter: 
                throw Exception("Error! There is no option to optimize only the character model during alternative model fitting"); 
                break;
            case onlySequence: 
                optimizeSequenceModel(); 
                break;
            case both: 
                throw Exception("Error! Attempt to violate the depedency between the chatacter model and sequnce model during alternative model fitting");
            default:
                throw Exception("Error! illegal optimizationScope setting");
        }

        /* free resources - now some of these parameters were defined localy - need to use friend functions otherwise can't free it */
        size_t numOfMappings = mappings.size();
        for (size_t h=0; h<numOfMappings; ++h)
        {
            if  (mappings[h]) delete mappings[h];
        }
        
        if (expectedHistory) delete expectedHistory; // delete the expectedHistory, that was cloned via updatesequenceTreeLikelihood
    }

    else 
    { // either nothing changed or only the seqluence parameters changed
        if (sequenceChanged_)
        {
            sequenceTreeLikelihood_->computeTreeLikelihood();
        }
    }

    /* for debugging purpose - report character model paramceters */
    TransitionModel* characterModel = characterTreeLikelihood_->getModel();
    ParameterList parameters = characterModel->getParameters();
    cout << "\n" << endl;
    for (size_t i = 0; i < parameters.size(); i++)
    {
        ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }
    ApplicationTools::displayResult("Character log likelihood", TextTools::toString(-1*characterTreeLikelihood_->getValue()));
    ApplicationTools::displayResult("Sequence log likelihood", TextTools::toString(-1*sequenceTreeLikelihood_->getValue()));
}


void JointLikelihoodFunction::setHypothesis(JointLikelihoodFunction::Hypothesis hypothesis)
{
    if (hypothesis == JointLikelihoodFunction::Hypothesis(1))
    {
        characterChanged_ = true; // need to alter this to trigger likelihood computation under the alternative model
        if (hypothesis_ != hypothesis && getParameterValue("RELAX.k_2") == 1)
        {
            setParameterValue("RELAX.k_2", previousK_);
        }
    }
    else
    {
        if (hypothesis_ != hypothesis)
        {
            previousK_ = getParameterValue("RELAX.k_2"); // update the last alternative value of k
            setParameterValue("RELAX.k_2", 1);
            sequenceChanged_ = true; // in case of switch from alternative to null, compute the sequence likelihood again under the null hypothesis
        }
    }
    hypothesis_ = hypothesis;
}

/******************************************************************************/

void JointLikelihoodFunction::optimizeCharacterModel()
{
    // optimize with BrentOneDimension, like in the alternative fitting
    BrentOneDimension* characterParametersOptimizer = new BrentOneDimension(characterTreeLikelihood_);
    characterParametersOptimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);
    characterParametersOptimizer->getStopCondition()->setTolerance(0.01); // set the tolerance to be slighly less strict to account for the instability of the joint likelihood function
    characterParametersOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    characterParametersOptimizer->setProfiler(0);
    characterParametersOptimizer->setMessageHandler(0);
    characterParametersOptimizer->setVerbose(1);
    ParameterList pi0; 
    pi0.addParameter(characterTreeLikelihood_->getParameter("TwoParameterBinary.pi0"));
    const IntervalConstraint* pi0Bounds = dynamic_cast<const IntervalConstraint*>(characterTreeLikelihood_->getParameter("TwoParameterBinary.pi0").getConstraint());
    characterParametersOptimizer->setInitialInterval(pi0Bounds->getLowerBound(), pi0Bounds->getUpperBound()); // search within stricter bounds that the actual ones of pi0 to avoid failute of stochasitc mapping
    characterParametersOptimizer->init(pi0);
    characterParametersOptimizer->optimize();
    ParameterList mu; 
    mu.addParameter(characterTreeLikelihood_->getParameter("TwoParameterBinary.mu"));
    const IntervalConstraint* muBounds = dynamic_cast<const IntervalConstraint*>(characterTreeLikelihood_->getParameter("TwoParameterBinary.mu").getConstraint());
    characterParametersOptimizer->setInitialInterval(muBounds->getLowerBound(), muBounds->getUpperBound()); // search within stricter bounds that the actual ones of pi0 to avoid failute of stochasitc mapping
    characterParametersOptimizer->init(mu);
    characterParametersOptimizer->optimize();
    delete characterParametersOptimizer;
    
    // update the log likelihood of the joint model
    double charLogL = characterTreeLikelihood_->getValue();
    double seqLogL = sequenceTreeLikelihood_->getValue();
    double overallLogL = charLogL + seqLogL; // = log(charLikelihood * seqLikelihood)
    logl_ = overallLogL;
}

/******************************************************************************/

void JointLikelihoodFunction::computeNullJointLikelihood()
{
    TransitionModel* characterModel = characterTreeLikelihood_->getModel();
    ParameterList parameters = characterModel->getParameters();
    switch(optimizationScope_)
    {
        case none:
            // because the sequence model doesn't depend on the character model, the likelihood can be computed separately
            if (characterChanged_)
            {
                characterTreeLikelihood_->computeTreeLikelihood();
            }
            if (sequenceChanged_)
            {
                sequenceTreeLikelihood_->computeTreeLikelihood();
            }
            /* for debugging purpose - report character model paramceters */
            cout << "\n" << endl;
            for (size_t i = 0; i < parameters.size(); i++)
            {
                ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
            }
            ApplicationTools::displayResult("Character log likelihood", TextTools::toString(-1*characterTreeLikelihood_->getValue()));
            ApplicationTools::displayResult("Sequence log likelihood", TextTools::toString(-1*sequenceTreeLikelihood_->getValue()));        
            break;
        case onlyCharacter: 
            optimizeCharacterModel();
            break;
        case onlySequence: 
            optimizeSequenceModel();
            break;
        case both:
            // optimize the character model independetly of the sequence model
            optimizeCharacterModel();
            // optimize the sequence model with hypothesis_ (which is null, since computeNullJointLikelihood was called)
            optimizeSequenceModel();
            break;
        default: 
            throw Exception("Error! illegal optimizationScope setting");
    }
}