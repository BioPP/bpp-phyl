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
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h> 
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>
#include <Bpp/Seq/Alphabet/NumericAlphabet.h>

// From bpp-phyl
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/G2001.h>
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Mapping/StochasticMapping.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Simulation/DetailedSiteSimulator.h>
#include <Bpp/Phyl/Simulation/NonHomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Simulation/SequenceSimulationTools.h>
#include <Bpp/Phyl/Mapping/RewardMappingTools.h>
#include <Bpp/Phyl/Mapping/Reward.h>
#include <Bpp/Phyl/Mapping/DecompositionReward.h>



using namespace bpp;
using namespace std;

#define STATE "state"
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
            size_t origNodeState = static_cast<size_t>(tl->getAlphabetStateAsInt(leafsStates->getSequence(nodeName).getValue(0)));
            const Node* nodeInMapping = (dynamic_cast<const TreeTemplate<Node>*>(mapping))->getNode(nodeName);
            size_t nodeStateInMapping = StochasticMapping::getNodeState(nodeInMapping);
            if ((nodeStateInMapping < 2) & (origNodeState != nodeStateInMapping)) // if the node's state corresaponds to a concrete character state (not unknown character or a combination of state and rate in the case of a markov modulated model)
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
            if (abs(branchLengthInMapping - origBranchLength) > 0.0001) // expected history fails in the root - but the branch length of the branch coming from the root is insignificant as it isn't included in the tree. what happened here?
            {
                throw Exception("branch lengths not maintained in mapping");
            }
        }
    }
    // make sure there are no two nodes in a row such that both don't exist in the base tree and both recieve the same state (indicator of illegal transition)
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
                size_t curNodeState = StochasticMapping::getNodeState(curNode);
                size_t nextNodeState = StochasticMapping::getNodeState(curNode->getFather());
                if (curNodeState == nextNodeState && ((curNode->getName()).find("mapping") != std::string::npos) && ((curNode->getFather()->getName()).find("mapping") != std::string::npos)) // such transition is permitted in case the father is the root

                {
                    throw Exception("illegal transitions in the mapping");
                }
                curNode = curNode->getFather();
            }
        }
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

void computePosteriors(VVDouble& posteriorProbabilities, Tree* baseTree, RHomogeneousTreeLikelihood* tl, map<int, size_t> nodeIdToIndex)
{
    // some auxiliiary variables

    size_t statesNum = tl->getNumberOfStates();
    const TransitionModel* model = tl->getModelForSite(0,0); // this calls assumes that all the sites and all the branches are assoiacted with the same node

    const SiteContainer* leafsStates = tl->getData();
    TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(baseTree); 
    vector<Node*> nodes = ttree->getNodes();
    // compute the posterior probabilities according to Felsenstein prunnig algorithm: for each node nodes[i] and state s compute: P(Data[leafs under node[i]]|node[i] has state s] 
    for (size_t i=0; i<nodes.size(); ++i) // traverse the tree in post-order

    {
        int nodeId = nodes[i]->getId();
        string nodeName = nodes[i]->getName();
        if (nodes[i]->isLeaf()) // if the node is a leaf, set the posterior probability of its state to 1, and the rest ot 0

        {
            size_t leafState = static_cast<int>(tl->getAlphabetStateAsInt(leafsStates->getSequence(nodeName).getValue(0)));
            for (size_t nodeState=0; nodeState<statesNum; ++nodeState)
		    {
                if (nodeState != leafState)
                {
                    posteriorProbabilities[nodeIdToIndex[nodeId]][nodeState] = 0;
                }
                else {
                    posteriorProbabilities[nodeIdToIndex[nodeId]][nodeState] = 1; 
                }
            }   
        }
        else                   // if the node is internal, follow the Felesenstein computation rule to compute the posterior probability

        {   
            double dataProb = 0;
            for (size_t nodeState=0; nodeState<statesNum; ++nodeState)
		    {
                double fullProb = 1;
                for (size_t j=0; j<(nodes[i]->getNumberOfSons()); ++j) // for each son of the node, sum over the probabilities of all its assignments given its father's state (i.e, nodeState)
                {
                    double sonProb = 0;
                    double bl = nodes[i]->getSon(j)->getDistanceToFather();
                    for(size_t sonState=0; sonState<statesNum; ++sonState)
                    {
                        sonProb += model->Pij_t(nodeState, sonState, bl) * posteriorProbabilities[nodeIdToIndex[nodes[i]->getSon(j)->getId()]][sonState];
                    }
                    fullProb *= sonProb;
                }
                posteriorProbabilities[nodeIdToIndex[nodeId]][nodeState] = fullProb;
                dataProb += posteriorProbabilities[nodeIdToIndex[nodeId]][nodeState];
            }
            // now, compute from the so far compued partial likelihoods the posterior probabilities by dividing by the probability of the data (prior(data)=1 in ML world)
            // because the sum of partial likelihoods is in fact the probablity of the data, it is sufficient to standardize the vector

            for (size_t nodeState=0; nodeState<statesNum; ++nodeState)
		    {
                posteriorProbabilities[nodeId][nodeState] = posteriorProbabilities[nodeIdToIndex[nodeId]][nodeState] / dataProb;
            }    
		}
    }
}


int main() 
{
    try

    {
        //fix seed for debugging purposes
        double seedUb = 10000000;
        double mySeed = RandomTools::giveRandomNumberBetweenZeroAndEntry(seedUb);
        RandomTools::setSeed(static_cast<long int>(mySeed));
        cout << "seed: " << mySeed << endl; // for debugging purposes in case the tester fails
        
        // create a binary model
        const BinaryAlphabet* alphabet = new BinaryAlphabet();
        double mu = 1.;
        double pi0 = 0.5;
        ReversibleSubstitutionModel* nestedModel = dynamic_cast<ReversibleSubstitutionModel*>(new TwoParameterBinarySubstitutionModel(alphabet,mu,pi0));
        DiscreteDistribution* rDist = new GammaDiscreteDistribution(2, 1, 1, true);
        SubstitutionModel* model = new G2001 (nestedModel, rDist);
        
        // process tree
        TreeTemplate<Node>* ttree = TreeTemplateTools::parenthesisToTree("(S15:0.85385,((S19:0.0854569,S16:0.139158):0.248594,(((((S12:0.0215813,S14:0.0122578):0.00733911,((S20:0.0133406,S18:0.02058):0.00622244,S8:0.0616991):0.00855007):0.0194517,S21:0.0361841):0.0260926,(S10:6.01257,S9:0.0572114):0.00432963):0.0582364,((((S11:0.00192042,((S13:0.00546429,S22:0.00413541):1e-06,S7:0.00544892):0.00223313):0.0224013,S6:0.0147796):0.0012621,(S24:1e-06,S23:1e-06):0.020303):0.0480321,((S2:0.0212492,((S1:0.029627,S3:0.322449):1e-06,S17:0.0303775):1e-06):0.0311297,(S5:0.00337913,S4:1e-06):0.0451854):0.00880453):0.0445887):0.133367):0.85385);");
        Tree* tree = dynamic_cast<Tree*>(ttree); 
        giveNamesToInternalNodes(tree); // give internal names to nodes in post-order
        vector<Node*> nodes = ttree->getNodes();

        map<int,size_t> nodeIdToIndex;
        for (size_t i=0; i<nodes.size(); ++i)
        {
            nodeIdToIndex[nodes[i]->getId()] = i; 
        }

        // process character data
        VectorSiteContainer sites(alphabet);
		sites.addSequence(BasicSequence("S1", "1", alphabet));
		sites.addSequence(BasicSequence("S2", "0", alphabet));
		sites.addSequence(BasicSequence("S3", "1", alphabet));
		sites.addSequence(BasicSequence("S4", "0", alphabet));
		sites.addSequence(BasicSequence("S5", "0", alphabet));
		sites.addSequence(BasicSequence("S6", "0", alphabet));
		sites.addSequence(BasicSequence("S7", "0", alphabet));
		sites.addSequence(BasicSequence("S8", "1", alphabet));
		sites.addSequence(BasicSequence("S9", "0", alphabet));
		sites.addSequence(BasicSequence("S10", "0", alphabet));
		sites.addSequence(BasicSequence("S11", "0", alphabet));
		sites.addSequence(BasicSequence("S12", "0", alphabet));
		sites.addSequence(BasicSequence("S13", "0", alphabet));
		sites.addSequence(BasicSequence("S14", "1", alphabet));
		sites.addSequence(BasicSequence("S15", "1", alphabet));
		sites.addSequence(BasicSequence("S16", "1", alphabet));
		sites.addSequence(BasicSequence("S17", "0", alphabet));
		sites.addSequence(BasicSequence("S18", "1", alphabet));
		sites.addSequence(BasicSequence("S19", "0", alphabet));
		sites.addSequence(BasicSequence("S20", "1", alphabet));
		sites.addSequence(BasicSequence("S21", "1", alphabet));
		sites.addSequence(BasicSequence("S22", "0", alphabet));
		sites.addSequence(BasicSequence("S23", "0", alphabet));
		sites.addSequence(BasicSequence("S24", "-", alphabet));
        SiteContainerTools::changeGapsToUnknownCharacters(sites);
        
        // create tree likelihood function
        DiscreteDistribution* LfrDist = new ConstantRateDistribution();
        RHomogeneousTreeLikelihood* characterTreeLikelihood = new RHomogeneousTreeLikelihood(*tree, dynamic_cast<const SiteContainer&>(sites), dynamic_cast<TransitionModel*>(model), LfrDist, false);
        characterTreeLikelihood->initialize();
        
        // generate 1000 sotchastic mappings
        unsigned int mappingsNum = 1000;
        StochasticMapping* stocMapping = new StochasticMapping(dynamic_cast<TreeLikelihood*>(characterTreeLikelihood), mappingsNum);
        vector<Tree*> mappings;
        stocMapping->generateStochasticMapping(mappings);
        
        // make sure all the mappings are legal
        for (size_t i=0; i<mappingsNum; ++i)
        {
            checkIfMappingLegal(stocMapping, mappings[i], ttree, characterTreeLikelihood);
        }

		// compute ancestral frequencies over the stochastic mappings
        VVDouble ancestralFrequencies;
        ancestralFrequencies.clear();
        size_t statesNum = characterTreeLikelihood->getNumberOfStates();
        ancestralFrequencies.resize(nodes.size(), VDouble(statesNum));
        stocMapping->computeStatesFrequencies(ancestralFrequencies, mappings);
        
		// generate an expected history
        Tree* expectedHistory = stocMapping->generateExpectedMapping(mappings);
        string treeStr = TreeTools::treeToParenthesis(*expectedHistory); // for debugging

        checkIfMappingLegal(stocMapping, expectedHistory, ttree, characterTreeLikelihood);
        
		// compute the average dwelling times at node S5
        VDouble AverageDwellingTimes;
        AverageDwellingTimes.clear();
        AverageDwellingTimes.resize(statesNum,0);
        
		// compute the average dwelling times of all the states
        for (size_t i=0; i<mappings.size(); ++i)
        {
            TreeTemplate<Node>* mapping =  dynamic_cast<TreeTemplate<Node>*>(mappings[i]);
            // get the pointers to the node and its father in the i'th mapping

            Node* curNode = mapping->getNode("S19");
            
            Node* father = mapping->getNode((ttree->getNode("S19"))->getFather()->getName()); // the original father of the node (according to the base tree) in the mapping

            while (curNode != father)
            {
                AverageDwellingTimes[StochasticMapping::getNodeState(curNode)] += curNode->getDistanceToFather();
                curNode = curNode->getFather();
            }
        }
        for (size_t s=0; s<statesNum; ++s)
        {
            AverageDwellingTimes[s] = 1.0* AverageDwellingTimes[s] / mappingsNum;
        }

		// check state of father of S5: if the state of the father is 1 -> also make sure the division of dwelling time under state 1 corresponds to the frequency of 1 in the father
        Node* son = dynamic_cast<TreeTemplate<Node>*>(expectedHistory)->getNode("S19");
		string fatherName = (ttree->getNode("S19"))->getFather()->getName();
		Node* father = dynamic_cast<TreeTemplate<Node>*>(expectedHistory)->getNode(fatherName);
        size_t fatherState = StochasticMapping::getNodeState(father);
        size_t sonState = StochasticMapping::getNodeState(son);
		double splitFromSon, splitToFather;
        map<size_t, double> stateToDwelling;
        splitFromSon = son->getDistanceToFather();
        stateToDwelling.clear();
        stateToDwelling[StochasticMapping::getNodeState(son->getFather())] = (son->getFather())->getDistanceToFather();
        stateToDwelling[StochasticMapping::getNodeState(son->getFather()->getFather())] = ((son->getFather())->getFather())->getDistanceToFather();  
        stateToDwelling[StochasticMapping::getNodeState(son->getFather()->getFather()->getFather())] = ((son->getFather())->getFather())->getFather()->getDistanceToFather(); 
        splitToFather = ((son->getFather())->getFather())->getFather()->getFather()->getDistanceToFather();  
        if (fatherState == sonState)
        {
            stateToDwelling[sonState] = splitFromSon + splitToFather;
            // compute division of dwelling time according to the states freuqncies at the father
            double fatherFrequency = ancestralFrequencies[nodeIdToIndex[(ttree->getNode("S19"))->getFather()->getId()]][sonState];
			double sonFrequency = ancestralFrequencies[nodeIdToIndex[(ttree->getNode("S19"))->getId()]][sonState];
            double fatherShare = fatherFrequency / (sonFrequency+fatherFrequency) * AverageDwellingTimes[sonState];
            double sonShare = AverageDwellingTimes[sonState] - fatherShare;

            // expected: ...S5{1}:AverageDwellingTimes[1]-fatherShare)mappingInternal_1{0}:AverageDwellingTimes[0])mappingInternal_2{1}:fatherShare)father{1}
            if (abs(splitFromSon - sonShare) > 0.0001)
            {
                cout << "Error in dwelling time division between father and son. Branch of son is of length " << splitFromSon << " instead of " << sonShare << endl;
                return 1;
            }
            if (abs(splitToFather - fatherShare) > 0.0001)
            {
                cout << "Error in dwelling time division between father and son. Branch beneath father is of length " << splitToFather << " instead of " << fatherShare << endl;
                return 1;               
            }
        }
        else
        {
            stateToDwelling[sonState] = splitFromSon;
            stateToDwelling[fatherState] = splitToFather;
            if (abs(stateToDwelling[sonState] - AverageDwellingTimes[sonState]) > 0.0001)
            {
                cout << "Error in dwelling time from son. Branch of son is of length " << stateToDwelling[sonState] << " instead of " << AverageDwellingTimes[sonState] << endl;
                return 1;
            }
            if (abs(stateToDwelling[fatherState] - AverageDwellingTimes[fatherState]) > 0.0001)
            {
                cout << "Error in dwelling time to father. Branch of son is of length " << stateToDwelling[fatherState] << " instead of " << AverageDwellingTimes[fatherState] << endl;
                return 1;
            }
        }
        // check the rest of the splits
        for (size_t s=0; s<statesNum; ++s)
        {
            if ((s != sonState) & (s != fatherState))
            {
                if (abs(stateToDwelling[s] - AverageDwellingTimes[s]) > 0.0001)
                {
                    cout << "Error in dwelling time in transition from state " << s << "Branch is of length " << stateToDwelling[s] << " instead of " << AverageDwellingTimes[s] << endl;
                    return 1;   
                }
            }
        }

        // delete all the created stochastic mappings
        for (size_t i=0; i<mappings.size(); ++i)
        {
            delete mappings[i];
        }

		// repeat the same tests for the analytic expected history
        Tree* analyticExpectedHistory = stocMapping->generateAnalyticExpectedMapping();
        checkIfMappingLegal(stocMapping, analyticExpectedHistory, ttree, characterTreeLikelihood);

        // create another analytic expected mapping and make sure its equal to the former one (reconstruction is deterministic)
        Tree* analyticExpectedHistory2 = stocMapping->generateAnalyticExpectedMapping();
        
		// compare for each node in post-order traversal: name, state and distance to father
        vector<Node*> hist1Nodes = dynamic_cast<TreeTemplate<Node>*>(analyticExpectedHistory)->getNodes();
        vector<Node*> hist2Nodes = dynamic_cast<TreeTemplate<Node>*>(analyticExpectedHistory2)->getNodes();
        if (hist1Nodes.size() != hist2Nodes.size())
        {
            cerr << "Error! in repeated reconstruction of analytic expected history the number of nodes varies" << endl;
            return 1;
        }
        for (size_t n=0; n<hist1Nodes.size(); ++n)
        {
            if (hist1Nodes[n]->getId() != hist2Nodes[n]->getId())
            {
                cerr << "Error! the nodes IDs in the two reconstrcued analyitc histories don't match for index " << n << endl;
                return 1;
            }
            if (StochasticMapping::getNodeState(hist1Nodes[n]) != StochasticMapping::getNodeState(hist2Nodes[n]))
            {
                cerr << "Error! the nodes states in the two reconstrcued analyitc histories don't match for node id " << hist1Nodes[n]->getId() << endl;
                return 1;
            }
            if (hist1Nodes[n]->getId() != analyticExpectedHistory->getRootId())
            {
                if (abs(hist1Nodes[n]->getDistanceToFather() - hist2Nodes[n]->getDistanceToFather()) > 0.0001)
                {
                    cerr << "Error! the branch lengths in the two reconstrcued analyitc histories don't match for node id " << hist1Nodes[n]->getId() << endl;
                    return 1;
                }
            }
        }
         
        delete characterTreeLikelihood;
        delete expectedHistory;
        delete analyticExpectedHistory;
        delete analyticExpectedHistory2;
        delete ttree;
        delete model; // rDist will be deleted via model as the destructor of G2001 deletes its rDist_ data member
        delete alphabet;
        delete stocMapping;
        delete LfrDist;

    }
    catch (exception & e)
    {
        cout << e.what() << endl;
        return 1;
    }
    return 0;
}