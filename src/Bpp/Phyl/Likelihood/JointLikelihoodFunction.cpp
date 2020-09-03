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
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
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
    RHomogeneousTreeLikelihood* characterTreeLikelihood = new RHomogeneousTreeLikelihood(*tree, *characterData, characterModel, rDist, false, true, true);
    characterTreeLikelihood->initialize();
    characterTreeLikelihood_ = characterTreeLikelihood;

    // create the sequence tree likelihood of the null model as initial data memeber based on the sequence model, sequence data and tree
    RNonHomogeneousMixedTreeLikelihood* sequenceTreeLikelihood = new RNonHomogeneousMixedTreeLikelihood(*tree, *sequenceData, sequenceModel, rDist, true, true);
    sequenceTreeLikelihood->initialize();
    sequenceTreeLikelihood_ = sequenceTreeLikelihood;

    // update the logl_ data member
    logl_ = characterTreeLikelihood_->getValue() + sequenceTreeLikelihood_->getValue();

    // create the stochadtic mapping data member based on the character likelihood function data member
    size_t numOfMappings = static_cast<size_t>(ApplicationTools::getIntParameter("character.num_of_mappings", bppml_->getParams(), 1000));
    stocMapping_ = new StochasticMapping(dynamic_cast<TreeLikelihood*>(characterTreeLikelihood_), numOfMappings);

    // set the previuos parameters list to the initial ones
    updatePreviousParametersValues();

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
    // update character model parameters
    ParameterList charParams = characterTreeLikelihood_->getSubstitutionModelParameters();
    for (size_t p=0; p < charParams.size(); ++p)
    {
       parName = charParams[p].getName();
       parValue = charParams[p].getValue();
       previousParametersValues_[parName] = parValue; 
    }
    // update sequence model parameters
    ParameterList seqParams = sequenceTreeLikelihood_->getSubstitutionModelParameters();
    for (size_t p=0; p < seqParams.size(); ++p)
    {
       parName = seqParams[p].getName();
       parValue = seqParams[p].getValue();
       previousParametersValues_[parName] = parValue;
    }
    // update log likelihood 
    previousParametersValues_["logl"] = -logl_;
}

/******************************************************************************/

void JointLikelihoodFunction::fireParameterChanged(const ParameterList& pl) 
{  
    // update the characterChanged_ value
    ParameterList charParams = characterTreeLikelihood_->getModelForSite(0, 0)->getParameters();
    double currValue, prevValue;
    for (size_t p=0; p<charParams.size(); ++p)
    {
        prevValue = previousParametersValues_[charParams[p].getName()];
        currValue = charParams[p].getValue();
        // check if pl[p] is in character model and if the value of the parameter changed
        if (abs(currValue - prevValue) > 0.000001)
        {
            characterChanged_ = true;
        }
    }

    // update the sequenceChanged_ value
    ParameterList seqParams = sequenceTreeLikelihood_->getSubstitutionModelParameters();
    sequenceChanged_ = false;
    for (size_t p=0; p<seqParams.size(); ++p)
    {
        prevValue = previousParametersValues_[seqParams[p].getName()];
        currValue = seqParams[p].getValue();
        // check if pl[p] is in character model and if the value of the parameter changed
        if (currValue != prevValue)
        {
            sequenceChanged_ = true;
        }
    }

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
    getModelParameters(false);

    // update the previous parameters list
    updatePreviousParametersValues();
}

/******************************************************************************/

void JointLikelihoodFunction::updateStatesInNodesNames(Tree* mapping)
{
    string label = "state";
    vector<Node*> nodes = (dynamic_cast<TreeTemplate<Node>*>(mapping))->getNodes(); 
    for (auto node: nodes) 
    {
        string name = node->getName();
        size_t state = StochasticMapping::getNodeState(node);
        node->setName(name + "{" + TextTools::toString(state) + "}");
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
        size_t nodeState = StochasticMapping::getNodeState(nodes[i]);
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
    const VectorSiteContainer* sequenceData = dynamic_cast<const VectorSiteContainer*>(sequenceTreeLikelihood_->getData());
    MixedSubstitutionModelSet* sequenceModel = dynamic_cast<MixedSubstitutionModelSet*>(sequenceTreeLikelihood_->getSubstitutionModelSet());
    DiscreteDistribution* rDist = sequenceTreeLikelihood_->getRateDistribution();

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

void JointLikelihoodFunction::scaleSequenceTree(double factor)
{
    // get a new tree and scale it
	const Tree& origTree = characterTreeLikelihood_->getTree(); // the character tree is taken because it was not affected by any previous scaling
	Tree* newTree = origTree.clone();
	(dynamic_cast<TreeTemplate<Node>*>(newTree))->scaleTree(factor);
    // switch the new tree with the ols tree in the sequence likelihood function
    updatesequenceTreeLikelihood(newTree);;
	getSequenceScalingFactor(true); // for debugging
}

/******************************************************************************/

map<string,double> JointLikelihoodFunction::getModelParameters(bool verbose)
{
    map<string,double> modelParameters;
    ParameterList parameters;
	
	// character model parameters
    TransitionModel* characterModel = characterTreeLikelihood_->getModel();
    parameters = characterModel->getParameters();
    for (size_t i = 0; i < parameters.size(); i++)
    {
        if (verbose)
        {
            ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
        }
        modelParameters[parameters[i].getName()] = parameters[i].getValue();
    }

  // sequence model parameters
  MixedSubstitutionModelSet* sequenceModel = dynamic_cast<MixedSubstitutionModelSet*>(sequenceTreeLikelihood_->getSubstitutionModelSet());
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
        modelParameters[parameters[i].getName() + "_" + TextTools::toString(m+1)] = parameters[i].getValue();
    }
  }

  // get the sequence scaling factor
  modelParameters["sequenceScalingFactor"] = getSequenceScalingFactor(false);
  if (verbose)
	ApplicationTools::displayResult("Sequence scaling factor",  TextTools::toString(modelParameters["sequenceScalingFactor"]));

  // get the likelihood
  double characterLogl = characterTreeLikelihood_->getValue();
  modelParameters["Character Log likelihood"] = characterLogl;
  double sequenceLogl = sequenceTreeLikelihood_->getValue();
  modelParameters["Sequence Log likelihood"] = sequenceLogl;
  logl_ = characterLogl + sequenceLogl;
  modelParameters["Overall Log likelihood"] = logl_;

  // report it regardless of verbose level
  if (verbose)
  {
	ApplicationTools::displayResult("\nCharacter Log likelihood", TextTools::toString(-characterLogl, 15));
	ApplicationTools::displayResult("Sequence Log likelihood", TextTools::toString(-sequenceLogl, 15));
	ApplicationTools::displayResult("Overall Log likelihood", TextTools::toString(-logl_, 15));
  }
  return modelParameters;
}

/*****************************************************************************/

void JointLikelihoodFunction::optimizeSequenceModel()
{
    unsigned int verbose = static_cast<unsigned int>(ApplicationTools::getIntParameter("optimization.verbose", bppml_->getParams(), 0));
	if (verbose)
		cout << "** Optimzing the sequence model **" << endl;

	// reset the parameters to ignore
	bppml_->getParam("optimization.ignore_parameters") = "BrLen";
	// now set the addional parameters you would like to ignore
	string initialParamsToIgnore = ApplicationTools::getStringParameter("optimization.ignore_parameters_additonal", bppml_->getParams(), "", "", true, true);
    string paramsToIgnore;
	bool isAlternative = false;
    // extract the user initial values sd initial starting point
    string BGModelInitialValues = ApplicationTools::getStringParameter("model1", bppml_->getParams(), "RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1)", "", true, true);
    string FGModelInitialValues = ApplicationTools::getStringParameter("model2", bppml_->getParams(), "RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1)", "", true, true);
    string modelName = "RELAX";
    map<string, string> BGArgs, FGArgs;
    KeyvalTools::parseProcedure(BGModelInitialValues, modelName, BGArgs);
    KeyvalTools::parseProcedure(FGModelInitialValues, modelName, FGArgs);
    map<string,double> userInitialValues;
    for (map<string, string>::iterator it = BGArgs.begin(); it != BGArgs.end(); it++)
	{
		if (it->first.compare("frequencies") !=0)
			userInitialValues[modelName + "." + it->first + "_1"] = TextTools::toDouble(it->second);
    }
    userInitialValues[modelName + ".k_2"] = TextTools::toDouble(FGArgs["k"]);
    switch(hypothesis_)
    {
        case null: 
            // set the value of RELAX.k_2 to 1 as well and then ignore
            sequenceTreeLikelihood_->setParameterValue("RELAX.k_2", 1);
            // now add RELAX.k_2 to the set of parameters to ignore
            paramsToIgnore = initialParamsToIgnore + ",BrLen,RELAX.k_1,RELAX.k_2";
            break;      // k = 1 both in model1 (BG) and model2 (FG)
        case alternative: 
            cycleNum_ = cycleNum_ + 1;
            isAlternative = true;
            paramsToIgnore = initialParamsToIgnore + ",BrLen,RELAX.k_1,RELAX.1_Full.theta_1,RELAX.1_Full.theta1_1,RELAX.1_Full.theta2_1,RELAX.2_Full.theta_1,RELAX.2_Full.theta1_1,RELAX.2_Full.theta2_1,RELAX.3_Full.theta_1,RELAX.3_Full.theta1_1,RELAX.3_Full.theta2_1"; // ignore frequency parameters to reduce optimization duration - results in one unit of ll reduction in optimality and 1 minutre reduction in duration
            break;      // k = 1 only in model1 (BG) 
        default:
            throw Exception("Error! illegal hypothesis setting");
    }
    bppml_->getParam("optimization.ignore_parameters") = paramsToIgnore;
	// first optimize the scaling parameter: 0 = don't scale, 1 - scale only before optimization, 2 - scale only after optimization, 3 - scale before and after optimization
    int scaleTree = ApplicationTools::getIntParameter("optimization.scale.tree", bppml_->getParams(), 0);
    if (scaleTree >= 1)
    {
        OptimizationTools::optimizeTreeScale(sequenceTreeLikelihood_, 0.000001, 1000000, ApplicationTools::message.get(), ApplicationTools::message.get(), 0);
        getSequenceScalingFactor(); // debug
    }

    // now optimize the rest of parameters
    map<string,double> inferenceResult;
	int advancedOptimization = ApplicationTools::getIntParameter("optimization.advanced", bppml_->getParams(), 0);
    double prevLogLikelihood, currLogLikelihood;
	size_t index;
    if ((advancedOptimization == 1) & (isAlternative))
    {
        bppml_->getParam("optimization.max_number_f_eval") = "100";
        bppml_->getParam("optimization.tolerance") = "0.01";
        if (cycleNum_ == 0)
        {
            // starting point 1 - results of the null fitting
			prevLogLikelihood = -sequenceTreeLikelihood_->getValue();
			currLogLikelihood = -sequenceTreeLikelihood_->getValue();
			index = 1;
			do
			{
				cout << "Optimization cycle: " << TextTools::toString(index) << endl;
				index = index + 1;
				if (scaleTree >= 1)
				{
					OptimizationTools::optimizeTreeScale(sequenceTreeLikelihood_, 0.000001, 1000000, ApplicationTools::message.get(), ApplicationTools::message.get(), 0);
					getSequenceScalingFactor(); // debug
				}
				PhylogeneticsApplicationTools::optimizeParameters(sequenceTreeLikelihood_, sequenceTreeLikelihood_->getParameters(), bppml_->getParams());
				ApplicationTools::displayResult("Current log likelihood", TextTools::toString(-sequenceTreeLikelihood_->getValue(), 15));
				prevLogLikelihood = currLogLikelihood;
				currLogLikelihood = -sequenceTreeLikelihood_->getValue();
				ApplicationTools::displayResult("Current diff", TextTools::toString((currLogLikelihood-prevLogLikelihood), 15));
			} while (currLogLikelihood - prevLogLikelihood > 0.01);
			if (verbose)
            {
                cout << "iterative optimzation complete" << endl;
                ApplicationTools::displayResult("Log likelihood", TextTools::toString(-sequenceTreeLikelihood_->getValue(), 15));
            }
			double sp1Logl = -sequenceTreeLikelihood_->getValue();
            if (verbose)
            {
                cout << "* Statring point: null fitting result *" << endl;
                ApplicationTools::displayResult("Log likelihood", TextTools::toString(sp1Logl, 15));
            }
            map<string, double> sp1Result = getModelParameters(verbose); // debug - print model parameters
            
            // starting point 2 - user initial values
            for (map<string, double>::iterator it = userInitialValues.begin(); it != userInitialValues.end(); it++)
            {
				if ((it->first.find("RELAX") != std::string::npos))
				{
					sequenceTreeLikelihood_->setParameterValue(it->first, it->second);
				}
            }
						prevLogLikelihood = -sequenceTreeLikelihood_->getValue();
			currLogLikelihood = -sequenceTreeLikelihood_->getValue();
			index = 1;
			do
			{
				if (verbose)
                    cout << "Optimization cycle: " << TextTools::toString(index) << endl;
				index = index + 1;
				if (scaleTree >= 1)
				{
					OptimizationTools::optimizeTreeScale(sequenceTreeLikelihood_, 0.000001, 1000000, ApplicationTools::message.get(), ApplicationTools::message.get(), 0);
					getSequenceScalingFactor(); // debug
				}
				PhylogeneticsApplicationTools::optimizeParameters(sequenceTreeLikelihood_, sequenceTreeLikelihood_->getParameters(), bppml_->getParams());
				ApplicationTools::displayResult("Current log likelihood", TextTools::toString(-sequenceTreeLikelihood_->getValue(), 15));
				prevLogLikelihood = currLogLikelihood;
				currLogLikelihood = -sequenceTreeLikelihood_->getValue();
				ApplicationTools::displayResult("Current diff", TextTools::toString((currLogLikelihood-prevLogLikelihood), 15));
			} while (currLogLikelihood - prevLogLikelihood > 0.01);
			if (verbose)
            {
                cout << "iterative optimzation complete" << endl;
			    ApplicationTools::displayResult("Log likelihood", TextTools::toString(-sequenceTreeLikelihood_->getValue(), 15));
            }
            double sp2Logl = -sequenceTreeLikelihood_->getValue();
            if (verbose)
            {
                cout << "* Statring point: user initial values *" << endl;
                ApplicationTools::displayResult("Log likelihood", TextTools::toString(sp2Logl, 15));
            }
            map<string, double> sp2Result = getModelParameters(verbose); // debug - print model parameters
            
            // determine the winning starting point
            map<string, double> chosenInitialValues = sp1Result;
            if (sp1Logl < sp2Logl)
            {
            if (verbose)
                cout << "Winning starting point: user initial values " << endl;
            chosenInitialValues = sp2Result;
            }
            else
            {
                if (verbose)
                    cout << "Winning starting point: null fitting result " << endl;
                chosenInitialValues = sp1Result;
            }
            // set the values of the starting point in the sequence likelihood function
            for (map<string, double>::iterator it = chosenInitialValues.begin(); it != chosenInitialValues.end(); it++)
            {
                if ((it->first.find("RELAX") != std::string::npos))
                {
                    sequenceTreeLikelihood_->setParameterValue(it->first, it->second);
                }
            }
        }
        
        vector<size_t> startingPointsByCycle = ApplicationTools::getVectorParameter<size_t>("optimization.starting_points_by_cycle", bppml_->getParams(), ',', "10,3,1", "", true);
        size_t numberOfCycles = startingPointsByCycle.size();

		/* step 1: fist cycle of optimizaton: only compute the likelihood over multiple starting points with respect to k */
        // set the starting points with respect to the value of k according to the number given in startingPointsByCycle[0]
        if (verbose)
            cout << "** Step 1: choosing initial starting points from a selection of " << startingPointsByCycle[0] << " points **" << endl;
        size_t numberOfRexalationPoints = static_cast<size_t>(floor(startingPointsByCycle[0]/2));
        size_t numberOfIntensificationPoints = startingPointsByCycle[0] - numberOfRexalationPoints;
        double relaxationInterval = 1 / static_cast<double>(numberOfRexalationPoints) - 0.0001;
        double bgOmega0 = sequenceTreeLikelihood_->getParameterValue("RELAX.p_1") * sequenceTreeLikelihood_->getParameterValue("RELAX.omega1_1");
        double bgOmega2 = sequenceTreeLikelihood_->getParameterValue("RELAX.omega2_1");
        double maxSigK = min(min(max(log(0.0001+0.0001)/log(bgOmega0+0.0001),1.0), max(log(999+0.0001)/log(bgOmega2+0.0001),1.0)), 10.0); // compute the maximal k for which the breakwater in RELAX model implementation is not expressed (any k beyond this value will yield the same likelihood)
        double intensificationInterval = (maxSigK-1) / static_cast<double>(numberOfIntensificationPoints);
        vector<double> startingPoints;
        startingPoints.clear();
        startingPoints.resize(startingPointsByCycle[0]);
        // set the starting points with respect to k s.t 1/2 of them represent relaxation and 1/2 represent intensification
        for (size_t r=0; r<startingPointsByCycle[0]; ++r)
        {
            if (r < numberOfRexalationPoints)
            {
                startingPoints[r] = 0.0001 + static_cast<double>(r)*relaxationInterval;
            }
            else
            {
                startingPoints[r] = 1 + static_cast<double>(r-numberOfRexalationPoints)*intensificationInterval;
            }
        }
        // compute the likelihood of the starting points
        vector<map<string,double>> startingPointsResults;
        startingPointsResults.clear();
        startingPointsResults.resize(startingPointsByCycle[0]); 
        for (size_t l=0; l<startingPointsByCycle[0]; ++l)
        {
            sequenceTreeLikelihood_->setParameterValue("RELAX.k_2", startingPoints[l]);
            sequenceTreeLikelihood_->computeTreeLikelihood();
            startingPointsResults[l] = getModelParameters(verbose);
        }
        
        /* step 2: iteratively optimize superficially the best starting points */
        if (verbose)
            cout << "** Step 2: iteratively optimizing superficially the best starting points for " << numberOfCycles-2 << " cycles **" << endl;
        bppml_->getParam("optimization.max_number_f_eval") = "100";
        bppml_->getParam("optimization.tolerance") = "0.01";
        vector<map<string,double>> bestStartingPoints;
        for (size_t c=1; c<numberOfCycles-1; ++c) // exclude the first cycle (choosing of starting pionts) and last cycle(deep optimizaiton on the best starting point)
        {
            if(verbose)
                cout << "* Step 2 Cycle " << c << " *" << endl;
            size_t numberOfBestStartingPoints = startingPointsByCycle[c];
            bestStartingPoints.clear();
            sort(startingPointsResults.begin(), startingPointsResults.end(), JointLikelihoodFunction::sortStartingPointsFunction);
            // select the best starting points with respect to k
            for (size_t p=0; p<numberOfBestStartingPoints; ++p)
            {
                bestStartingPoints.push_back(startingPointsResults[p]);
            }
            // optimize each point and save the output of its optimization
            for (size_t p=0; p<bestStartingPoints.size(); ++p)
            {
                for (map<string, double>::iterator it = bestStartingPoints[p].begin(); it != bestStartingPoints[p].end(); it++)
                {
                    if ((it->first.find("RELAX") != std::string::npos))
                    {
                        sequenceTreeLikelihood_->setParameterValue(it->first, it->second);
                    }
                }
                if (verbose)
                    cout << "Optimizing starting point " << (p+1) << "..." << endl;
				prevLogLikelihood = -sequenceTreeLikelihood_->getValue();
				currLogLikelihood = -sequenceTreeLikelihood_->getValue();
				index = 1;
				do
				{
					if (verbose)
                        cout << "Optimization cycle: " << TextTools::toString(index) << endl;
					index = index + 1;
					if (scaleTree >= 1)
					{
						OptimizationTools::optimizeTreeScale(sequenceTreeLikelihood_, 0.000001, 1000000, ApplicationTools::message.get(), ApplicationTools::message.get(), 0);
						getSequenceScalingFactor(); // debug
					}
					PhylogeneticsApplicationTools::optimizeParameters(sequenceTreeLikelihood_, sequenceTreeLikelihood_->getParameters(), bppml_->getParams());
					if (verbose)
                        ApplicationTools::displayResult("Current log likelihood", TextTools::toString(-sequenceTreeLikelihood_->getValue(), 15));
					prevLogLikelihood = currLogLikelihood;
					currLogLikelihood = -sequenceTreeLikelihood_->getValue();
					if (verbose)
                        ApplicationTools::displayResult("Current diff", TextTools::toString((currLogLikelihood-prevLogLikelihood), 15));
				} while (currLogLikelihood - prevLogLikelihood > 0.01);
				if (verbose)
                    cout << "iterative optimzation complete" << endl;
				ApplicationTools::displayResult("Log likelihood", TextTools::toString(-sequenceTreeLikelihood_->getValue(), 15));
				bestStartingPoints[p] = getModelParameters(verbose);
            }
            startingPointsResults.clear();
            startingPointsResults = bestStartingPoints; 
        }

        /* step 3: deep optimization of the best startingPointsByCycle[startingPointsByCycle.size()-1] starting points */
        if (verbose)
            cout << "** Step 3: optimizating deeply the optimal starting points **" << endl;
        bppml_->getParam("optimization.max_number_f_eval") = "10000";
        bppml_->getParam("optimization.tolerance") = "0.000001";
        size_t numberOfBestStartingPoints = startingPointsByCycle[startingPointsByCycle.size()-1];
        bestStartingPoints.clear();
        sort(startingPointsResults.begin(), startingPointsResults.end(), JointLikelihoodFunction::sortStartingPointsFunction);
        // select the best starting points with respect to k
        for (size_t p=0; p<numberOfBestStartingPoints; ++p)
        {
            bestStartingPoints.push_back(startingPointsResults[p]);
            // optimize deeply each of the chosen starting points
            for (map<string, double>::iterator it = startingPointsResults[p].begin(); it != startingPointsResults[p].end(); it++)
            {
                if ((it->first.find("RELAX") != std::string::npos))
                {
                    sequenceTreeLikelihood_->setParameterValue(it->first, it->second);
                }
            }
            if (verbose)
                cout << "Optimizing starting point " << (p+1) << endl;
			prevLogLikelihood = -sequenceTreeLikelihood_->getValue();
			currLogLikelihood = -sequenceTreeLikelihood_->getValue();
			index = 1;
			do
			{
				if (verbose)
                    cout << "Optimization cycle: " << TextTools::toString(index) << endl;
				index = index + 1;
				if (scaleTree >= 1)
				{
					OptimizationTools::optimizeTreeScale(sequenceTreeLikelihood_, 0.000001, 1000000, ApplicationTools::message.get(), ApplicationTools::message.get(), 0);
					getSequenceScalingFactor(); // debug
				}
				PhylogeneticsApplicationTools::optimizeParameters(sequenceTreeLikelihood_, sequenceTreeLikelihood_->getParameters(), bppml_->getParams());
				if (verbose)
                    ApplicationTools::displayResult("Current log likelihood", TextTools::toString(-sequenceTreeLikelihood_->getValue(), 15));
				prevLogLikelihood = currLogLikelihood;
				currLogLikelihood = -sequenceTreeLikelihood_->getValue();
				if (verbose)
                    ApplicationTools::displayResult("Current diff", TextTools::toString((currLogLikelihood-prevLogLikelihood), 15));
			} while (currLogLikelihood - prevLogLikelihood > 0.01);
			if (verbose)
            {
                cout << "iterative optimzation complete" << endl;
			    ApplicationTools::displayResult("Log likelihood", TextTools::toString(-sequenceTreeLikelihood_->getValue(), 15));bestStartingPoints[p] = getModelParameters(verbose);
            }
        }

        /* step 4: select the best starting point and report its values */ 
        sort(bestStartingPoints.begin(), bestStartingPoints.end(), sortStartingPointsFunction);
        inferenceResult = bestStartingPoints[0];
    }
	else
    {
		prevLogLikelihood = -sequenceTreeLikelihood_->getValue();
		currLogLikelihood = -sequenceTreeLikelihood_->getValue();
		index = 1;
		do
		{
			if (verbose)
                cout << "Optimization cycle: " << TextTools::toString(index) << endl;
			index = index + 1;
			if (scaleTree >= 1)
			{
				OptimizationTools::optimizeTreeScale(sequenceTreeLikelihood_, 0.000001, 1000000, ApplicationTools::message.get(), ApplicationTools::message.get(), 0);
				getSequenceScalingFactor(); // debug
			}
			PhylogeneticsApplicationTools::optimizeParameters(sequenceTreeLikelihood_, sequenceTreeLikelihood_->getParameters(), bppml_->getParams());
			if (verbose)
                ApplicationTools::displayResult("Current log likelihood", TextTools::toString(-sequenceTreeLikelihood_->getValue(), 15));
			prevLogLikelihood = currLogLikelihood;
			currLogLikelihood = -sequenceTreeLikelihood_->getValue();
			if (verbose)
                ApplicationTools::displayResult("Current diff", TextTools::toString((currLogLikelihood-prevLogLikelihood), 15));
		} while (currLogLikelihood - prevLogLikelihood > 0.01);
		if (verbose)
            cout << "iterative optimzation complete" << endl;
		ApplicationTools::displayResult("Log likelihood", TextTools::toString(-sequenceTreeLikelihood_->getValue(), 15));
        inferenceResult = getModelParameters(verbose);
    }

    // set the optimal parameters to the function and then report tge results of the joint model
    for (map<string, double>::iterator it = inferenceResult.begin(); it != inferenceResult.end(); it++)
    {
        if ((it->first.find("RELAX") != std::string::npos))
        {
            sequenceTreeLikelihood_->setParameterValue(it->first, it->second);
        }  
    }
    sequenceTreeLikelihood_->computeTreeLikelihood();
	if (verbose)
	{
		cout << "\n** Model parameters after sequence model optimizaiton ** " << endl;
		getModelParameters(verbose);
	}
	
    // update the log likelihood of the joint model
    double charLogL = characterTreeLikelihood_->getValue();
    double seqLogL = sequenceTreeLikelihood_->getValue();
    double overallLogL = charLogL + seqLogL; // = log(charLikelihood * seqLikelihood)
    logl_ = overallLogL;
}

/******************************************************************************/

void JointLikelihoodFunction::computeAlternativeJointLikelihood()
{  

	unsigned int verbose = static_cast<unsigned int>(ApplicationTools::getIntParameter("optimization.verbose", bppml_->getParams(), 0));
	
    if (characterChanged_)
    {
        /* compute likelihood of charcter model */
        characterTreeLikelihood_->computeTreeLikelihood();
        
        /* approximate the expected character history based on numOfMappings sampled stochastic mappings */
        bool useAnalytic =  static_cast<bool>(ApplicationTools::getIntParameter("character.use_analytic_mapping", bppml_->getParams(), 1));
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
            // for each mapping, set the parittion according ot it and define it as a tree of the cloned sequence likelihood function
            ApplicationTools::displayResult("Character model log likelihood: ", TextTools::toString(-characterTreeLikelihood_->getValue(), 15));
            cout << "Computing sequence log likelihoods given the different mappings\n" << endl;
            for (size_t h=0; h<mappings.size(); ++h)
            {
                    setPartitionByHistory(mappings[h]); // induce a partition of the tree based on the epxected character history
                    updatesequenceTreeLikelihood(mappings[h]); // compute the likelihood given the mapping
                    cout << TextTools::toString(-characterTreeLikelihood_->getValue() - sequenceTreeLikelihood_->getValue(), 15) << endl;
            
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
            double charLogl = -characterTreeLikelihood_->getValue();
            double sequenceLogl = -sequenceTreeLikelihood_->getValue();
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
		sequenceChanged_ = true; 						// since the partition changed, the sequence likelihood has also changed
	
		/* free resources - now some of these parameters were defined localy - need to use friend functions otherwise can't free it */
		size_t numOfMappings = mappings.size();
		for (size_t h=0; h<numOfMappings; ++h)
		{
			if  (mappings[h]) delete mappings[h];
		}
		
		if (expectedHistory) delete expectedHistory; // delete the expectedHistory, that was cloned via updatesequenceTreeLikelihood

		// either nothing changed or only the sequence parameters changed
		if (sequenceChanged_)
		{
			sequenceTreeLikelihood_->computeTreeLikelihood();
		}
	}
	
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
	
	// report for debugging
	if ((characterChanged_) & (optimizationScope_ == JointLikelihoodFunction::OptimizationScope(0)))
	{
		/* for debugging purpose - report character model parameters */
		TransitionModel* characterModel = characterTreeLikelihood_->getModel();
		ParameterList parameters = characterModel->getParameters();
		for (size_t i = 0; i < parameters.size(); i++)
		{
			ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue(), 15));
		}
		// get the likelihood
		double characterLogl = characterTreeLikelihood_->getValue();
		double sequenceLogl = sequenceTreeLikelihood_->getValue();
		logl_ = characterLogl + sequenceLogl;
        if (verbose)
        {
            ApplicationTools::displayResult("Character Log likelihood", TextTools::toString(-characterLogl, 15));
		    ApplicationTools::displayResult("Sequence Log likelihood", TextTools::toString(-sequenceLogl, 15));
            ApplicationTools::displayResult("Overall Log likelihood", TextTools::toString(-logl_, 15));
        }
	}
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
	unsigned int verbose = static_cast<unsigned int>(ApplicationTools::getIntParameter("optimization.verbose", bppml_->getParams(), 0));
	if (verbose)
		cout << "** Optimzing the character model **" << endl; 

	double prevLogLikelihood = -characterTreeLikelihood_->getValue();
    double currLogLikelihood = -characterTreeLikelihood_->getValue();
	// set optimization method to FullD(Newton) which showed convergence abilities
	string initialParamsToIgnore = ApplicationTools::getStringParameter("optimization.ignore_parameters", bppml_->getParams(), "", "", true, true);
	bppml_->getParam("optimization.ignore_parameters") = initialParamsToIgnore + ",BrLen";
	string prevOptimization = ApplicationTools::getStringParameter("optimization", bppml_->getParams(), "FullD(derivatives=Newton)", "", true, true);
	bppml_->getParam("optimization") = "FullD(derivatives=Newton)";
	do
	{
		prevLogLikelihood = -characterTreeLikelihood_->getValue();
		PhylogeneticsApplicationTools::optimizeParameters(characterTreeLikelihood_, characterTreeLikelihood_->getParameters(), bppml_->getParams());
		currLogLikelihood = -characterTreeLikelihood_->getValue();
	} while (currLogLikelihood - prevLogLikelihood > 0.01);
	// switch back to input optimization method which will be used for the sequence model
	bppml_->getParam("optimization") = prevOptimization;
	
	    // update the log likelihood of the joint model
    double charLogL = characterTreeLikelihood_->getValue();
    double seqLogL = sequenceTreeLikelihood_->getValue();
    double overallLogL = charLogL + seqLogL; // = log(charLikelihood * seqLikelihood)
    logl_ = overallLogL;
	
	if (verbose)
	{
		cout << "\n" << endl;
		ParameterList parameters = characterTreeLikelihood_->getSubstitutionModelParameters();
		for (size_t p=0; p<parameters.size(); ++p)
			ApplicationTools::displayResult(parameters[p].getName(), TextTools::toString(parameters[p].getValue()));
		ApplicationTools::displayResult("Character log likelihood", TextTools::toString(-characterTreeLikelihood_->getValue()));
	}
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
            ApplicationTools::displayResult("Character log likelihood", TextTools::toString(-characterTreeLikelihood_->getValue(), 15));
            ApplicationTools::displayResult("Sequence log likelihood", TextTools::toString(-sequenceTreeLikelihood_->getValue(), 15));        
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

double JointLikelihoodFunction::getLikelihood()
{
    double seqLikelihood = sequenceTreeLikelihood_->getLikelihood();
    double charLikelihood = characterTreeLikelihood_->getLikelihood();
    return seqLikelihood + charLikelihood;
}

vector<double> JointLikelihoodFunction::getLikelihoodForEachSite()
{
    vector<double> seqLikelihoodBySite = sequenceTreeLikelihood_->getLikelihoodForEachSite();
    double charLikelihood = characterTreeLikelihood_->getLikelihood();
    vector<double> jointLikelihoodBySite;
    for (size_t s=0; s<seqLikelihoodBySite.size(); ++s)
    {
        jointLikelihoodBySite.push_back(charLikelihood*seqLikelihoodBySite[s]);
    }
    return jointLikelihoodBySite;
}