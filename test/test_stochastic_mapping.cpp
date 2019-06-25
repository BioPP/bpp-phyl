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
#include <Bpp/Numeric/Random/RandomTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h> 
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>

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
                int curNodeState = stocMapping->getNodeState(curNode);
                int nextNodeState = stocMapping->getNodeState(curNode->getFather());
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





void setNodeState(Node* node, size_t state)
{
    BppInteger* stateProperty = new BppInteger(static_cast<int>(state));
    node->setNodeProperty(STATE, *stateProperty);
    delete stateProperty;
}



void computePosteriors(VVDouble& posteriorProbabilities, Tree* baseTree, RHomogeneousTreeLikelihood* tl)
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
                    posteriorProbabilities[nodeId][nodeState] = 0;
                }
                else {
                    posteriorProbabilities[nodeId][nodeState] = 1; 
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
                        sonProb += model->Pij_t(nodeState, sonState, bl) * posteriorProbabilities[nodes[i]->getSon(j)->getId()][sonState];
                    }
                    fullProb *= sonProb;
                }
                posteriorProbabilities[nodeId][nodeState] = fullProb;
                dataProb += posteriorProbabilities[nodeId][nodeState];
            }
            // now, compute from the so far compued partial likelihoods the posterior probabilities by dividing by the probability of the data (prior(data)=1 in ML world)
            // because the sum of partial likelihoods is in fact the probablity of the data, it is sufficient to standardize the vector

            for (size_t nodeState=0; nodeState<statesNum; ++nodeState)
		    {
                posteriorProbabilities[nodeId][nodeState] = posteriorProbabilities[nodeId][nodeState] / dataProb;
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
        // create a binary model

        const BinaryAlphabet* alphabet = new BinaryAlphabet();
        double mu = 1.;
        double pi0 = 0.5;
        SubstitutionModel* model = new TwoParameterBinarySubstitutionModel(alphabet,mu,pi0);
        // process tree

        TreeTemplate<Node>* ttree = TreeTemplateTools::parenthesisToTree("(((S1:0,S2:0):1,S3:2):1,(S4:1,S5:1):2);");
        Tree* tree = dynamic_cast<Tree*>(ttree); 
        giveNamesToInternalNodes(tree); // give internal names to nodes in post-order

        vector<Node*> nodes = ttree->getNodes();
        vector<int> nodeIds = ttree->getNodesId();
  
        // process character data

        VectorSiteContainer sites(alphabet);
        sites.addSequence(BasicSequence("S1", "0", alphabet));
        sites.addSequence(BasicSequence("S2", "0", alphabet));
        sites.addSequence(BasicSequence("S3", "0", alphabet));
        sites.addSequence(BasicSequence("S4", "0", alphabet));
        sites.addSequence(BasicSequence("S5", "1", alphabet));
        // create tree likelihood function

        DiscreteDistribution* rDist = new ConstantRateDistribution();
        RHomogeneousTreeLikelihood* characterTreeLikelihood = new RHomogeneousTreeLikelihood(*tree, dynamic_cast<const SiteContainer&>(sites), dynamic_cast<TransitionModel*>(model), rDist, false);
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
        // compute posterior probabilies

        VVDouble posteriorProbabilities;
        posteriorProbabilities.clear();
        posteriorProbabilities.resize(nodes.size(), VDouble(2));
        computePosteriors(posteriorProbabilities, tree, characterTreeLikelihood);
        // make sure that the posterior probability of state 0 in the parent of S1,S2 is 1

        Node* parent = (ttree->getNode("S1"))->getFather();
        double probability = posteriorProbabilities[parent->getId()][0];
        if (probability != 1)
        {
            cout << "Error in computation of posterior probability at node (S1,S2):X. Computed probability is " << probability << " instead of 1" << endl;
        }
        // compute ancestral frequencies over the stochastic mappings

        VVDouble ancestralFrequencies;
        ancestralFrequencies.clear();
        ancestralFrequencies.resize(nodes.size(), VDouble(2));
        double singletonFraction = 1.0/mappingsNum;
        for (size_t i=0; i<mappingsNum; ++i)
        {
            for (size_t j=0; j < nodes.size(); ++j)
            {
                Node* mappingNode = dynamic_cast<TreeTemplate<Node>*>(mappings[i])->getNode(nodes[j]->getName());
                int state = stocMapping->getNodeState(mappingNode);
                if (state == 0)
                {
                    ancestralFrequencies[j][0] = ancestralFrequencies[j][0] + singletonFraction; 
                }
                else

                {
                    ancestralFrequencies[j][1] = ancestralFrequencies[j][1] + singletonFraction; 
                }
            }
        }
        // generate an expected history

        Tree* expectedHistory = stocMapping->generateExpectedMapping(mappings);
        string treeStr = TreeTools::treeToParenthesis(*expectedHistory); // for debugging

        checkIfMappingLegal(stocMapping, expectedHistory, ttree, characterTreeLikelihood);
        Node* expectedMappingNode;
        int state;
        double stateFrequency, stateProbability;
        // make sure ancestral assignments correspond to frequencies in the mappings and the the posterior probabilities

        for (size_t j=0; j < nodes.size(); ++j)
        {
           if (!nodes[j]->isLeaf())
           {
            expectedMappingNode = dynamic_cast<TreeTemplate<Node>*>(expectedHistory)->getNode(nodes[j]->getName()); 
            state = stocMapping->getNodeState(expectedMappingNode);
            stateFrequency = ancestralFrequencies[j][state];
            if (stateFrequency < 0.5)
            {
                    cout << "Failed to assign ancestral state to node " << expectedMappingNode->getName() << " according to the frequency: Assigned state is " << state << " while its frequency is " << stateFrequency << endl;
                    return 1;
            }
            stateProbability = posteriorProbabilities[j][state];
            if (stateProbability < 0.5)
            {
                cout << "Frequency of state " << state << " at node " << expectedMappingNode->getName() << " doesn't agree with the posterior probability: frequency is " << stateFrequency << " while probability is " << stateProbability << endl;
                return 1;
            }
           }
        }
        // compute the average dwelling times at node S5

        VDouble AverageDwellingTimes;
        AverageDwellingTimes.clear();
        AverageDwellingTimes.resize(2,0);
        // compute the average dwelling times of all the states

        for (size_t i=0; i<mappings.size(); ++i)
        {
            TreeTemplate<Node>* mapping =  dynamic_cast<TreeTemplate<Node>*>(mappings[i]);
            // get the pointers to the node and its father in the i'th mapping

            Node* curNode = mapping->getNode("S5");
            
            Node* father = mapping->getNode((ttree->getNode("S5"))->getFather()->getName()); // the original father of the node (according to the base tree) in the mapping

            while (curNode != father)
            {
                AverageDwellingTimes[stocMapping->getNodeState(curNode)] += curNode->getDistanceToFather();
                curNode = curNode->getFather();
            }
        }
        for (size_t s=0; s<2; ++s)
        {
            AverageDwellingTimes[s] = 1.0* AverageDwellingTimes[s] / mappingsNum;
        }
        // check state of father of S5: if the state of the father is 1 -> also make sure the division of dwekking time under state 1 cirresponds to the frequency of 1 in the father

        Node* son = dynamic_cast<TreeTemplate<Node>*>(expectedHistory)->getNode("S5");
        string fatherName = (ttree->getNode("S5"))->getFather()->getName();
        Node* father = dynamic_cast<TreeTemplate<Node>*>(expectedHistory)->getNode(fatherName);
        int fatherState = stocMapping->getNodeState(father);
        int sonState = stocMapping->getNodeState(son);
        double split1, split2, split3;
        if (fatherState == sonState)
        {
            split1 = son->getDistanceToFather();
            split2 = (son->getFather())->getDistanceToFather();
            split3 = ((son->getFather())->getFather())->getDistanceToFather();      
            // compute division of dwelling time according to the states freuqncies at the father

            double fatherFrequency = ancestralFrequencies[(ttree->getNode("S5"))->getFather()->getId()][1];
            double fatherShare = fatherFrequency / (1+fatherFrequency) * AverageDwellingTimes[1];
            double sonShare = AverageDwellingTimes[1] - fatherShare;
            // expected: ...S5{1}:AverageDwellingTimes[1]-fatherShare)mappingInternal_1{0}:AverageDwellingTimes[0])mappingInternal_2{1}:fatherShare)father{1}
            if (abs(split1 - sonShare) > 0.0001)
            {
                cout << "Error in dwelling time division between father and son. Branch of son is of length " << split1 << " instead of " << sonShare << endl;
                return 1;
            }
            if (split2 != AverageDwellingTimes[0])
            {
                cout << "Error in dwelling time assignment under state 0 to splitting point of branch S5: branch length is " << split2 << " instead of " << AverageDwellingTimes[0] << endl;
                return 1;
            }
            if (abs(split3 - fatherShare) > 0.0001)
            {
                cout << "Error in dwelling time division between father and son. Branch beneath father is of length " << split3 << " instead of " << fatherShare << endl;
                return 1;               
            }
        }
        else

        {
            split1 = son->getDistanceToFather();
            split2 = (son->getFather())->getDistanceToFather();
            // expected: ...S5{1}:AverageDwellingTimes[1])mappingInternal{0}:AverageDwellingTimes[0])father{0}
            if (split1 != AverageDwellingTimes[1])
            {
                cout << "Error in dwelling time assignment under state 1 to splitting point of branch S5: branch length is " << split1 << " instead of " << AverageDwellingTimes[1] << endl;
                return 1;
            }
            if (split2 != AverageDwellingTimes[0])
            {
                cout << "Error in dwelling time assignment under state 0 to splitting point of branch S5: branch length is " << split2 << " instead of " << AverageDwellingTimes[0] << endl;
                return 1;
            }
        }
        // repeat the same tests for the analytic expected history

        Tree* analyticExpectedHistory = stocMapping->generateAnalyticExpectedMapping();
        checkIfMappingLegal(stocMapping, analyticExpectedHistory, ttree, characterTreeLikelihood);
        // make sure ancestral assignments correspond to the posterior probabilities

        for (size_t j=0; j < nodes.size(); ++j)
        {
           if (!nodes[j]->isLeaf())
           {
            expectedMappingNode = dynamic_cast<TreeTemplate<Node>*>(analyticExpectedHistory)->getNode(nodes[j]->getName()); 
            state = stocMapping->getNodeState(expectedMappingNode);
            stateProbability = posteriorProbabilities[j][state];
            if (stateProbability < 0.5)
            {
                cout << "Starte assignment at node " << expectedMappingNode->getName() << " doesn't agree with the posterior probability: probability is " << stateProbability << endl;
                return 1;
            }
           }
        }
        // compute the expected dwelling times at node S5

        VDouble ExpectedDwellingTimes;
        ExpectedDwellingTimes.clear();
        ExpectedDwellingTimes.resize(2,0);
        UserAlphabetIndex1* alpha = new UserAlphabetIndex1(characterTreeLikelihood->getAlphabet());
        Node* node = dynamic_cast<TreeTemplate<Node>*>(tree)->getNode("S5");
        DRTreeLikelihood* drtl = new DRHomogeneousTreeLikelihood(*tree, dynamic_cast<const SiteContainer&>(sites), dynamic_cast<TransitionModel*>(model), rDist, false);
        drtl->initialize();
        vector <int> ids;
        ids.push_back(node->getId());
        const vector<int> states =  drtl->getAlphabetStates();
        for (uint s=0; s<states.size(); ++s)
        {
            alpha->setIndex(s,1);
            for (uint m=0; m<states.size(); ++m)
            {
                if (m != s)
                {
                    alpha->setIndex(m,0);
                }
            } 
            unique_ptr<Reward> reward(new DecompositionReward(dynamic_cast<const SubstitutionModel*>(model), alpha));
            unique_ptr<ProbabilisticRewardMapping> mapping(RewardMappingTools::computeRewardVectors(*drtl, ids, *reward, false));  
            ExpectedDwellingTimes[s] = mapping->getReward(node->getId(), 0);
        }
        // make sure the sum of the expected dwelling times sums up to the branch length

        double sumOfExpectedDwellingTimes = ExpectedDwellingTimes[0] + ExpectedDwellingTimes[1];
        double branchLength = node->getDistanceToFather();
        if (sumOfExpectedDwellingTimes != branchLength)
        {
            cerr << "Error! sum of expected dwelling times under all states of the model doesn't add up to the original length of the branch. Sum is: " << sumOfExpectedDwellingTimes << "instead of " << branchLength << endl;
            return 1; 
        }
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
            if (stocMapping->getNodeState(hist1Nodes[n]) != stocMapping->getNodeState(hist2Nodes[n]))
            {
                cerr << "Error! the nodes states in the two reconstrcued analyitc histories don't match for node id " << hist1Nodes[n]->getId() << endl;
                return 1;
            }
            if (hist1Nodes[n]->getId() != analyticExpectedHistory->getRootId())
            {
                if (hist1Nodes[n]->getDistanceToFather() != hist2Nodes[n]->getDistanceToFather())
                {
                    cerr << "Error! the branch lengths in the two reconstrcued analyitc histories don't match for node id " << hist1Nodes[n]->getId() << endl;
                    return 1;
                }
            }
        }
         
    
        // delete all the created instances

        for (size_t i=0; i<mappings.size(); ++i)
        {
            delete mappings[i];
        }
        delete rDist;
        delete characterTreeLikelihood;
        delete expectedHistory;
        delete analyticExpectedHistory;
        delete analyticExpectedHistory2;
        delete ttree;
        delete model;
        delete alphabet;
        delete stocMapping;
    }
    catch (exception & e)
    {
        cout << e.what() << endl;
        return 1;
    }
    return 0;
}