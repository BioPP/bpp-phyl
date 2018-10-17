// Tester for StochasticMapping implementation

// From the STL:
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// From bpp-core:
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h> 
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SequenceTools.h>

// From bpp-phyl
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Mapping/StochasticMapping.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Simulation/DetailedSiteSimulator.h>
#include <Bpp/Phyl/Simulation/NonHomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Simulation/SequenceSimulationTools.h>

using namespace bpp;
using namespace std;

#define STATE "state"
#define WRITE_PATH "/media/sf_sf_Biopp/data/mappings/"


/******************************************************************************/

void checkIfMappingLegal(const StochasticMapping* stocMapping, const Tree* mapping, const TreeTemplate<Node>* baseTree, const TreeLikelihood* tl)
{
    // make sure that the states are the leafs are set properly
    vector<const Node*> origNodes = baseTree->getNodes();
    vector<const Node*> mappingNodes = dynamic_cast<const TreeTemplate<Node>*>(mapping)->getNodes();
    const SiteContainer* leafsStates = tl->getData();
    for (size_t i=0; i<origNodes.size(); ++i)
    {
        if (origNodes[i]->isLeaf())
        {
            string nodeName = origNodes[i]->getName();
            int origNodeState = static_cast<int>(tl->getAlphabetStateAsInt(leafsStates->getSequence(nodeName).getValue(0)));
            const Node* nodeInMapping = (dynamic_cast<const TreeTemplate<Node>*>(mapping))->getNode(nodeName);
            int nodeStateInMapping = stocMapping->getNodeState(nodeInMapping);
            if (origNodeState != nodeStateInMapping)
            {
                 throw Exception("Leafs states not maintained in mapping");
            }
        }
    }
            
    // make sure that branch lengths are maintained in the mapping
    for (size_t i=0; i<origNodes.size(); ++i)
    {
        if (origNodes[i]->hasFather())
        {
            double origBranchLength = origNodes[i]->getDistanceToFather();
            const Node* nodeInMapping = (dynamic_cast<const TreeTemplate<Node>*>(mapping))->getNode(origNodes[i]->getName());
            const Node* origNodeFatherInMapping = (dynamic_cast<const TreeTemplate<Node>*>(mapping))->getNode(origNodes[i]->getFather()->getName());
            double branchLengthInMapping = 0;
            const Node* curNode = nodeInMapping;
            while (curNode != origNodeFatherInMapping)
            {
                if (curNode->hasFather())
                {
                    branchLengthInMapping += curNode->getDistanceToFather(); // failes on node with id 9, but returns exception on node with id 6 (which is the root)
                    curNode = curNode->getFather();
                }
                else
                {
                    branchLengthInMapping = origBranchLength;
                    curNode = origNodeFatherInMapping;
                }
            }
            if (branchLengthInMapping > origBranchLength + 0.01 || branchLengthInMapping < origBranchLength - 0.01) // expected history fails in the root - but the branch length of the branch coming from the root is insignificant as it isn't included in the tree. what happened here?
            {
                throw Exception("branch lengths not maintained in mapping");
            }
        }
    }

    // make sure there are no two nodes in a raw such that both don't exist in the base tree and both recieve the same state (indicator of illegal transition)
    for (size_t i=0; i<origNodes.size(); ++i)
    {
        if (origNodes[i]->hasFather())
        {
            const Node* nodeInMapping = (dynamic_cast<const TreeTemplate<Node>*>(mapping))->getNode(origNodes[i]->getName());
            string origFatherName = origNodes[i]->getFather()->getName();
            const Node* origNodeFatherInMapping = (dynamic_cast<const TreeTemplate<Node>*>(mapping))->getNode(origFatherName);
            const Node* curNode = nodeInMapping;
            while (curNode->getFather() != origNodeFatherInMapping)
            {
                string curNodeFatherName = curNode->getFather()->getName();
                string origNodeFatherInMappingName = origNodeFatherInMapping->getName();
                int curNodeState = stocMapping->getNodeState(curNode);
                int nextNodeState = stocMapping->getNodeState(curNode->getFather());
                if (curNodeState == nextNodeState && origNodes[i]->getFather()->hasFather()) // such transition is permitted in case the father is the root
                {
                    throw Exception("illegal transitions in the mapping");
                }
                curNode = curNode->getFather();
            }
        }
    }
}


void printMapping(const StochasticMapping* stocMapping, const Tree* mapping) 
{
    ApplicationTools::displayMessage("\n\nPrinting mapping\n");
    string label = "state";
    vector<const Node*> nodes = (dynamic_cast<const TreeTemplate<Node>*>(mapping))->getNodes(); 
    for (int i=0; i < static_cast<int>(nodes.size()); i++) {
        if (nodes[i]->hasFather()) 
        {
        string name = nodes[i]->getName();
        double bl = nodes[i]->getDistanceToFather();
        int state = stocMapping->getNodeState(nodes[i]);
        ApplicationTools::displayMessage("Node (name: " + TextTools::toString(name) + ", branch_length: " + TextTools::toString(bl) + ", state: " + TextTools::toString(state) + ")");
        }
    }
}


void printModelParameters(const DiscreteRatesAcrossSitesTreeLikelihood* tl)
{
    ParameterList parameters = tl->getSubstitutionModelParameters();
    for (size_t i = 0; i < parameters.size(); i++)
    {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }
    parameters = tl->getRateDistributionParameters();
    for (size_t i = 0; i < parameters.size(); i++)
    {
      ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }
    ParameterList pl = tl->getBranchLengthsParameters();
    for(unsigned int i = 0; i < pl.size(); i++)
    {
      ApplicationTools::displayResult(pl[i].getName(), TextTools::toString(pl[i].getValue())); 
    }
}


void giveNamesToInternalNodes(Tree* tree)
{
    TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(tree);
    vector<Node*> nodes = ttree->getNodes();
    for (size_t i=0; i<nodes.size(); ++i) {
        if (!nodes[i]->hasName())
            nodes[i]->setName("_baseInternal_" + TextTools::toString(nodes[i]->getId()));
    }  
}


int getNodeState(const Node* node)
{
    return (dynamic_cast<const BppInteger*>(node->getNodeProperty(STATE)))->getValue(); // exception on root on the true history - why didn't the root recieve a state?
}


void setNodeState(Node* node, size_t state)
{
    node->setNodeProperty(STATE,BppInteger(static_cast<int>(state)));
}


void updateStatesInNodesNames(Tree* mapping)
{
    string label = "state";
    vector<Node*> nodes = (dynamic_cast<TreeTemplate<Node>*>(mapping))->getNodes(); 
    for (int i=0; i < static_cast<int>(nodes.size()); i++) 
    {
        string name = nodes[i]->getName();
        int state = getNodeState(nodes[i]);
        nodes[i]->setName(name + "{" + TextTools::toString(state) + "}");
    }
}


// since the function is fitted to true history, the mutation path may exceed the original branch length, in which case the dutration until the last event should be reduced to fit the original branch length
void updateBranchMapping(Node* son, const MutationPath& branchMapping, size_t initial_count)
{ 
    static size_t nodesCounter_ = initial_count;
    const vector<size_t> states = branchMapping.getStates();
    const VDouble times = branchMapping.getTimes();
    Node* curNode = son;
    double origBranchLength = son->getDistanceToFather();
    Node* nextNode;
    int eventsNum = static_cast<int>(branchMapping.getNumberOfEvents());
    if (eventsNum == 0)  // if there are no events -> return nothing
    {
        return;
    }
    else
    {
        // update the branch length of the node to be as per the time of the last event
        double duration = 0;
        for (int i=eventsNum-2; i>-1; --i) 
        { // add a new node to represent the transition
            nodesCounter_ += 1;
            const string name = "_mappingInternal" + TextTools::toString(nodesCounter_) + "_";
            nextNode = new Node(static_cast<int>(nodesCounter_), name); // nodes counter isn't incremented properly
            setNodeState(nextNode,states[i]);
            if (i == 0)
            {
                duration = times[i];    
            }
            else
            {
                duration = times[i]-times[i-1]; // the duration is the gap between the two adjacent times
            }
            nextNode->setDistanceToFather(duration);    // according to the simulator, the documented time is the accumulated time until the transition and NOT the time since the last transition
            
            // set the father to no longer be the father of curNode
            if (curNode->hasFather()) // as long as we haven't reached the root -> upate the father of the current node
            {
                Node* originalFather = curNode->getFather();
                originalFather->removeSon(curNode); // also removes originalFather as the father of curNode
                // set nextNode to be the new father of curNode
                curNode->setFather(nextNode); // also adds curNode to the sons of nextNode
                // set curNode's original father ot be the father of nextNode
                nextNode->setFather(originalFather); // also adds nextNode to the sons of originalFather - or not? make sure this is the father at all times
                curNode = nextNode;
            }
        }
        // if the sum of distances doesn't add up to the length of the original branch, add an event or change the duration of the last event, depending on the last assigned state
        double lastEventDuration = origBranchLength - times[eventsNum-2]; 
        if (lastEventDuration > 0)
        {
            son->setDistanceToFather(lastEventDuration);
        }
    }
    return;
}



// translate the simulation data into a tree format, similar to the one given by the stochastic mappings
// in the results, there is a mutations paths field that can be returned for each node getMutationPath(int nodeId)
// then, it is enough to use the function I built in StochasticMapping
Tree* extractTrueHistory(RASiteSimulationResult* simulationData, Tree* baseTree)
{
    Tree* history = baseTree->clone();
    giveNamesToInternalNodes(history);
    vector<Node*> nodes = (dynamic_cast<TreeTemplate<Node>*>(history))->getNodes();
    vector<Node*> debugNodes;
    // set the ancestral states for all the nodes
    for (size_t j=0; j<nodes.size(); ++j)
    {
        size_t state = simulationData->getAncestralState(nodes[j]->getId());
        setNodeState(nodes[j], state);
    }
    for (size_t i=0; i<nodes.size(); ++i)
    {
        string name = nodes[i]->getName();
        if (nodes[i]->hasFather())
        {
            const MutationPath branchHistory = simulationData->getMutationPath(nodes[i]->getId());
            updateBranchMapping(nodes[i], branchHistory, nodes.size()-1);
        } 
    }
    // now, add the label of each internal node in the updated history to the node's name
    updateStatesInNodesNames(history);
    debugNodes = (dynamic_cast<TreeTemplate<Node>*>(history))->getNodes(); // for debugging
    return history;
}


void evaluateApproxmation(RHomogeneousTreeLikelihood* tl, const size_t NUM_OF_MAPPINGS,  TreeTemplate<Node>* tree)
{
        std::map<std::string, std::string> writeParams;
        string file_path;
        writeParams["output.tree.format"] = "Newick";
        cout << "**** generating stochasting mappings ****" << endl;
        StochasticMapping* stocMapping = new StochasticMapping(dynamic_cast<TreeLikelihood*>(tl), NUM_OF_MAPPINGS);
        vector<Tree*> mappings;
        stocMapping->generateStochasticMapping(mappings);

        // make sure the mapping is legal. and if it is, print it
        for (size_t i=0; i<NUM_OF_MAPPINGS; ++i)
        {
            checkIfMappingLegal(stocMapping, mappings[i], tree, tl);
            //printMapping(stocMapping, mappings[i]);
        }
        cout << "" << endl;
        
        // average the stochastic mappings
        cout << "**** generating expected mapping from the sampled stochasting mappings by method ((Pf/Pf+Ps),(Ps/Pf+Ps)) on division of dwelling time of state shared by son and father ****" << endl;
        Tree* expectedMapping_1 = stocMapping->generateExpectedMapping(mappings, 0); // bug in the expected history -branch lengths are not maintained
        checkIfMappingLegal(stocMapping, expectedMapping_1, tree, tl);
        printMapping(stocMapping, expectedMapping_1);
        cout << "" << endl;

        cout << "**** generating expected mapping from the sampled stochasting mappings by method (0%,100%) on division of dwelling time of state shared by son and father ****" << endl;
        Tree* expectedMapping_2 = stocMapping->generateExpectedMapping(mappings, 1); // bug in the expected history -branch lengths are not maintained
        checkIfMappingLegal(stocMapping, expectedMapping_2, tree, tl);
        //printMapping(stocMapping, expectedMapping_2);
        cout << "" << endl;
 
 
        // write and then delete all the mappings
        for (size_t i=0; i<NUM_OF_MAPPINGS; ++i)
        {
            updateStatesInNodesNames(mappings[i]);
            file_path = WRITE_PATH + TextTools::toString(NUM_OF_MAPPINGS) + string("_histories/sm/history_") + TextTools::toString(i) + string(".nwk");
            writeParams["output.tree.file"] = file_path;
            PhylogeneticsApplicationTools::writeTree(*mappings[i], writeParams);
            delete mappings[i];
        }
        mappings.clear();
        // write and then delete the expected mapping
        updateStatesInNodesNames(expectedMapping_1);
        file_path =  WRITE_PATH + TextTools::toString(NUM_OF_MAPPINGS) + string("_histories/expected_history_div_Ps_Pf.nwk");
        writeParams["output.tree.file"] = file_path;
        PhylogeneticsApplicationTools::writeTree(*expectedMapping_1, writeParams); // node properties are not written - need to create a function that adds the properties to the nodes names (hopefully internal nodes names are written as well)
        delete expectedMapping_1;

        updateStatesInNodesNames(expectedMapping_2);
        file_path =  WRITE_PATH + TextTools::toString(NUM_OF_MAPPINGS) + string("_histories/expected_history_div_0_100.nwk");
        writeParams["output.tree.file"] = file_path;
        PhylogeneticsApplicationTools::writeTree(*expectedMapping_2, writeParams); // node properties are not written - need to create a function that adds the properties to the nodes names (hopefully internal nodes names are written as well)
        delete expectedMapping_2;
        delete stocMapping;
}


int main() 
{
    try
    {
        // process the base tree
        //TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree("(((((taxon_76:0.496929,((taxon_107:0.318441,(taxon_97:0.0159542,taxon_112:0.0159542):0.302486):0.0913119,(taxon_115:0.0372323,taxon_64:0.0372323):0.37252):0.0871764):0.512717,((taxon_100:0.348635,taxon_16:0.348635):0.586186,taxon_30:0.934821):0.074825):0.922849,((((taxon_84:0.165705,taxon_90:0.165705):0.878595,((taxon_63:0.831182,taxon_34:0.831182):0.0861682,taxon_49:0.91735):0.12695):0.0704198,((taxon_54:0.308919,taxon_103:0.308919):0.00428333,taxon_35:0.313203):0.801518):0.535358,(((((((taxon_62:0.0114359,taxon_113:0.0114359):0.183922,taxon_5:0.195357):0.277249,taxon_27:0.472607):0.149748,((taxon_89:0.379836,(taxon_105:0.116668,taxon_79:0.116668):0.263168):0.186337,taxon_96:0.566173):0.056181):0.26954,(taxon_21:0.653842,taxon_68:0.653842):0.238051):0.440175,taxon_47:1.33207):0.108835,(((((taxon_102:0.310821,(((((taxon_117:0.078882,taxon_109:0.078882):0.0435369,taxon_24:0.122419):0.0257292,taxon_40:0.148148):0.1129,(taxon_93:0.114194,(taxon_31:0.0114579,taxon_110:0.0114579):0.102736):0.146854):0.00819116,taxon_120:0.269239):0.0415812):0.159443,((taxon_116:0.333268,(taxon_8:0.274037,(taxon_38:0.10085,(taxon_127:0.0628677,(taxon_2:0.0302643,taxon_86:0.0302643):0.0326033):0.0379826):0.173187):0.0592309):0.0527145,(taxon_19:0.164283,taxon_125:0.164283):0.221699):0.0842811):0.538943,((((taxon_83:0.125342,taxon_60:0.125342):0.198194,(taxon_122:0.248157,(taxon_33:0.0833234,taxon_57:0.0833234):0.164834):0.0753784):0.0645332,taxon_106:0.388069):0.604844,((taxon_41:0.436068,taxon_69:0.436068):0.524348,taxon_51:0.960416):0.0324962):0.0162942):0.144357,taxon_92:1.15356):0.0883542,((taxon_108:0.0902716,taxon_11:0.0902716):0.330309,taxon_87:0.420581):0.821337):0.198986):0.209175):0.282416):0.0675047,((((((((taxon_99:0.0682924,taxon_50:0.0682924):0.643115,(((taxon_104:0.287754,taxon_65:0.287754):0.362259,((taxon_14:0.384113,(taxon_91:0.32686,taxon_101:0.32686):0.0572534):0.00244688,(taxon_25:0.11595,taxon_81:0.11595):0.27061):0.263453):0.0290233,((taxon_124:0.334751,(taxon_121:0.126325,taxon_85:0.126325):0.208425):0.187795,(taxon_48:0.120423,(taxon_123:0.037176,taxon_58:0.037176):0.0832468):0.402123):0.156491):0.032371):0.122432,(taxon_3:0.190266,taxon_67:0.190266):0.643573):0.343819,(taxon_45:0.0614213,taxon_28:0.0614213):1.11624):0.0281366,((taxon_32:0.344014,(taxon_7:0.258017,taxon_74:0.258017):0.0859975):0.501839,(taxon_111:0.66751,(((taxon_94:0.0880479,taxon_15:0.0880479):0.515293,((taxon_73:0.121536,(taxon_52:0.00849519,taxon_70:0.00849519):0.113041):0.00709992,taxon_23:0.128636):0.474705):0.0516189,(taxon_80:0.650887,(taxon_59:0.519428,(taxon_78:0.461878,(taxon_88:0.0870246,taxon_61:0.0870246):0.374853):0.0575505):0.131458):0.00407329):0.0125497):0.178344):0.359941):0.0329343,(((taxon_46:0.0215311,taxon_6:0.0215311):0.103302,taxon_20:0.124833):0.719305,(taxon_13:0.422967,(taxon_37:0.293719,taxon_66:0.293719):0.129248):0.421171):0.39459):0.0578141,((((taxon_82:0.451795,((((taxon_56:0.106496,((taxon_72:1.08453e-05,taxon_126:1.08453e-05):0.00491296,taxon_55:0.0049238):0.101572):0.0144704,taxon_75:0.120966):0.0983091,(taxon_26:0.0578083,(taxon_1:0.0332433,taxon_119:0.0332433):0.0245649):0.161467):0.0318726,((taxon_9:0.123689,taxon_39:0.123689):0.076446,taxon_95:0.200135):0.0510127):0.200647):0.172515,(taxon_18:0.315929,(taxon_53:0.226351,taxon_42:0.226351):0.0895782):0.308381):0.0930432,taxon_17:0.717353):0.499559,(((taxon_22:0.0629585,taxon_98:0.0629585):0.0807454,taxon_71:0.143704):0.869988,(taxon_43:0.312529,taxon_77:0.312529):0.701163):0.20322):0.0796314):0.649234,(((taxon_128:0.104818,taxon_10:0.104818):0.509612,(taxon_12:0.110661,taxon_36:0.110661):0.503769):0.858453,((((taxon_4:0.141369,taxon_118:0.141369):0.267746,taxon_44:0.409116):0.195081,taxon_29:0.604197):0.0821644,taxon_114:0.686361):0.786522):0.472894):0.0542229):0.407829);");
        TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree("((S1:1,S2:1):1,(S3:1,S4:1):1):1;");
        Tree& ttree = dynamic_cast<Tree&>(*tree);
        giveNamesToInternalNodes(&ttree);
        // write the base tree to the main analysis path 
        std::map<std::string, std::string> writeParams;
        string file_path = WRITE_PATH + string("base_tree.nwk");
        writeParams["output.tree.file"] = file_path;
        writeParams["output.tree.format"] = "Newick";
        PhylogeneticsApplicationTools::writeTree(ttree, writeParams);

        vector<Node*> nodes = tree->getNodes();

        // create a binary model
        const BinaryAlphabet* alphabet = new BinaryAlphabet();
        double mu = 2.;
        double pi0 = 0.5;
        SubstitutionModel* model = new TwoParameterBinarySubstitutionModel(alphabet,mu,pi0); // second arguent stands for mu
        DiscreteDistribution* rdist = new ConstantRateDistribution();

        // 15.8.18 - simulate character history using a simulator over a simple binary model
        NonHomogeneousSequenceSimulator* simulator = new NonHomogeneousSequenceSimulator(model, rdist, &ttree);
        //simulator->outputInternalSequences(true); // get the internal sequences in the sites container as well for later comparison
        vector<string> seqNames = tree->getLeavesNames();
        VectorSiteContainer sites(seqNames, alphabet);
        
        RASiteSimulationResult* result = simulator->dSimulateSite();
        unique_ptr<Site> site(result->getSite(*simulator->getSubstitutionModelSet()->getModel(0)));
        site->setPosition(0);
        sites.addSite(*site, false);
        Tree* trueHistory = extractTrueHistory(result, &ttree);
        // write the true history to a file
        file_path = WRITE_PATH + string("true_history.nwk");
        writeParams["output.tree.file"] = file_path;
        PhylogeneticsApplicationTools::writeTree(*trueHistory, writeParams);
        delete trueHistory;
        // write the character data into a file
        file_path = WRITE_PATH + string("char_data.fas");
        writeParams["output.sequence.file"] = file_path;
        writeParams["output.sequence.format"] = "Fasta";
        SequenceApplicationTools::writeSequenceFile(*(dynamic_cast<SequenceContainer*>(&sites)), writeParams, "", false, 1);

        // create a likelihood function
        RHomogeneousTreeLikelihood* tl = new RHomogeneousTreeLikelihood(ttree, dynamic_cast<const SiteContainer&>(sites), dynamic_cast<TransitionModel*>(model), rdist,false);
        tl->initialize();
        ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl->getValue(), 15));
        printModelParameters(tl);
        
        // evaluate the sotchastic mappinfs based approxmations
        cout<< "***** approximation based on 100 stochastic mappings *****" << endl;
        evaluateApproxmation(tl, 100, tree);

        cout<< "***** approximation based on 1000 stochastic mappings *****" << endl;
        evaluateApproxmation(tl, 1000, tree);

        cout<< "***** approximation based on 10000 stochastic mappings *****" << endl;
        evaluateApproxmation(tl, 10000, tree);

        // delete all the created instances
        delete tl;
        delete rdist;
        delete tree;
        delete model;
        delete alphabet;
        delete result;
        delete simulator;
    }
    catch (exception & e)
    {
        cout << e.what() << endl;
        return 1;
    }
    return 0;
}