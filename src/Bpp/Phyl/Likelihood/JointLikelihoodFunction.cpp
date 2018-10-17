#include "JointLikelihoodFunction.h"

// for bpp-core
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/MixedSubstitutionModel.h>
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
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
#include <math.h>       /* pow */
#include <string>

using namespace bpp;
using namespace std;

/******************************************************************************/

JointLikelihoodFunction::JointLikelihoodFunction(BppApplication* bppml, const Tree* tree, const VectorSiteContainer* characterData, TransitionModel* characterModel, const VectorSiteContainer* sequenceData, MixedSubstitutionModelSet* sequenceModel, DiscreteDistribution* rDist):
AbstractParametrizable(""), 
hypothesis_(null),
bppml_(bppml),
characterTreeLikelihood_(),
sequenceTreeLikelihood_(),
stocMapping_(),
optimizationScope_(none),
logl_(0)
{
    // create the character likelihood function as data member based on the character model, character data and tree
    RHomogeneousTreeLikelihood* characterTreeLikelihood = new RHomogeneousTreeLikelihood(*tree, *characterData, characterModel, rDist, false);
    characterTreeLikelihood->initialize();
    characterTreeLikelihood_ = characterTreeLikelihood;

    // create the sequence tree likelihood of the null model as initial data memeber based on the sequence model, sequence data and tree
    RNonHomogeneousMixedTreeLikelihood* sequenceTreeLikelihood = new RNonHomogeneousMixedTreeLikelihood(*tree, *sequenceData, sequenceModel, rDist, true, true);
    sequenceTreeLikelihood->initialize();
    sequenceTreeLikelihood_ = sequenceTreeLikelihood;

    // create the stochadtic mapping data member based on the character likelihood function data member
    size_t numOfMappings = static_cast<size_t>(ApplicationTools::getIntParameter("character.num_of_mappings", bppml_->getParams(), 100));
    stocMapping_ = new StochasticMapping(dynamic_cast<TreeLikelihood*>(characterTreeLikelihood_), numOfMappings);

    // set the parameters of the joint likelihood function instance as clones of the character model and sequence model parameters
    ParameterList characterParameters = characterModel->getParameters();
    for (size_t i=0; i<characterParameters.size(); ++i)
    {
        addParameter_(characterParameters[i].clone());
    }
    ParameterList sequenceParameters = sequenceModel->getParameters();
    for (size_t j=0; j<sequenceParameters.size(); ++j)
    {
        addParameter_(sequenceParameters[j].clone());
    }
    ParameterList paramsToUpdate = getParameters();
    fireParameterChanged(paramsToUpdate);
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

    // deletion of the parameters created via the constrcutor us done by the the deletion of the inheriting class AbstractParametrizable
}

/******************************************************************************/

void JointLikelihoodFunction::fireParameterChanged(const ParameterList& pl) 
{
    // update the values of the character model and sequence model, if required, before calling the function that computes the joint likelihood
    for (size_t i=0; i<pl.size(); ++i)
    {
        if (hasParameter((pl[i]).getName()))
        {
            if ((pl[i]).getName().compare("Binary.rate") == 0) // enters here all the time
            {
                characterTreeLikelihood_->setParameterValue(pl[i].getName(),(pl[i]).getValue());
            }
            else if ((pl[i]).getName().compare("Binary.kappa") == 0)
            {
                characterTreeLikelihood_->setParameterValue(pl[i].getName(),(pl[i]).getValue());
            } else { // the parameter belongs to the seuqence model
                // if the parameter name ends with _1 (i.e, belongs to model1), or it is the selection intensity parameter, set is value (all other parameters are aliased and thus their value needen't be set)
                string paramName = pl[i].getName();
                if (paramName.substr(paramName.length()-2, 2).compare("_1") == 0 || paramName.substr(0,8).compare("RELAX.k_") == 0)
                {
                    sequenceTreeLikelihood_->setParameterValue(pl[i].getName(),(pl[i]).getValue());
                }
            }
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
}

/******************************************************************************/

void JointLikelihoodFunction::setPartitionByHistory(const Tree* history)
{
    sequenceTreeLikelihood_->getSubstitutionModelSet()->resetModelToNodels();
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

void JointLikelihoodFunction::reportResults()
{
    // report character model parameters and likelihood
    ParameterList parameters;

    TransitionModel* characterModel = characterTreeLikelihood_->getModel();
    parameters = characterModel->getParameters();
    bool printToStd = false;
    switch(optimizationScope_)
    {
        case none: 
            break; 
        case onlyCharacter:
            break;
        case onlySequence: 
            printToStd = true;
            break;
        case both: 
            printToStd = true;
            break;
        default:
            throw Exception("Error! illegal hypothesis setting");
    }
    if (printToStd)
    {
        for (size_t i = 0; i < parameters.size(); i++)
        {
            ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
        }
    }
    double charLogL = characterTreeLikelihood_->getValue();
    ApplicationTools::displayResult("Character Log likelihood", TextTools::toString(-charLogL, 15));

    // report sequence model paraneters
    MixedSubstitutionModelSet* sequenceModel = dynamic_cast<MixedSubstitutionModelSet*>(sequenceTreeLikelihood_->getSubstitutionModelSet());
    if (printToStd)
    {
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
    ApplicationTools::displayResult("Sequence Log likelihood", TextTools::toString(-seqLogL, 15));

    // report joint likelihood
    double overallLogL = charLogL + seqLogL; // = log(charLikelihood * seqLikelihood)
    logl_ = overallLogL;
    ApplicationTools::displayResult("Overall Log likelihood", TextTools::toString(-overallLogL, 15));
    cout << "\n" << endl;
}

/******************************************************************************/

void JointLikelihoodFunction::optimizeSequenceModel()
{
    string paramValue;
    try
    {
        paramValue = bppml_->getParam("sequence.optimization.method");
    }
    catch (exception& e)
    {
        paramValue = "D-BFGS(derivatives=Gradient,nstep=10)";
    }
    string paramsToIgnore;
    switch(hypothesis_)
    {
        case null: 
            paramsToIgnore = "BrLen,RELAX.k_1,RELAX.k_2";
            break;      // k = 1 both in model1 (BG) and model2 (FG)
        case alternative: 
            paramsToIgnore = "BrLen,RELAX.k_1,RELAX.1_Full.theta_1,RELAX.1_Full.theta1_1,RELAX.1_Full.theta2_1,RELAX.2_Full.theta_1,RELAX.2_Full.theta1_1,RELAX.2_Full.theta2_1,RELAX.3_Full.theta_1,RELAX.3_Full.theta1_1,RELAX.3_Full.theta2_1"; // ignore frequency parameters to reduce optimization duration - results in one unit of ll reduction in optimality and 1 minutre reduction in duration
            break;      // k = 1 only in model1 (BG)
        default:
            throw Exception("Error! illegal hypothesis setting");
    }
    bppml_->getParam("optimization.ignore_parameters") = paramsToIgnore;
    PhylogeneticsApplicationTools::optimizeParameters(sequenceTreeLikelihood_, sequenceTreeLikelihood_->getParameters(), bppml_->getParams());
}

/******************************************************************************/

void JointLikelihoodFunction::computeAlternativeJointLikelihood()
{  
    /* approximate the expected character history based on numOfMappings sampled stochastic mappings */
    vector<Tree*> mappings;
    stocMapping_->generateStochasticMapping(mappings);
    Tree* expectedHistory = stocMapping_->generateExpectedMapping(mappings);

    setPartitionByHistory(expectedHistory); // induce a partition of the tree based on the epxected character history

    /* compute the likelihood of the sequence model */
    updatesequenceTreeLikelihood(expectedHistory);

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

    reportResults();

    /* free resources - now some of these parameters were defined localy - need to use friend functions otherwise can't free it */
    size_t numOfMappings = mappings.size();
    for (size_t h=0; h<numOfMappings; ++h)
    {
        delete mappings[h];
    }
    delete expectedHistory; // delete the expectedHistory, that was cloned via updatesequenceTreeLikelihood
}

/******************************************************************************/

void JointLikelihoodFunction::optimizeCharacterModel()
{
        // fill in missing parameter for optimization in bppml_
        string paramValue;
        try
        {
            paramValue = bppml_->getParam("sequence.optimization.method");
        }
        catch (exception& e)
        {
            paramValue = "D-Brent(derivatives=Newton, nstep=10)";
        }
        bppml_->getParam("optimization.ignore_parameters") = "BrLen"; // branch lengths should always be ignored
        PhylogeneticsApplicationTools::optimizeParameters(characterTreeLikelihood_, characterTreeLikelihood_->getParameters(), bppml_->getParams());
}

/******************************************************************************/

void JointLikelihoodFunction::computeNullJointLikelihood()
{
    switch(optimizationScope_)
    {
        case none:
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

    reportResults();
}