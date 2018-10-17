#include "StochasticMapping.h"
#include "../Likelihood/TreeLikelihoodTools.h"
#include "../Likelihood/TreeLikelihood.h"
#include "../Simulation/MutationProcess.h"
#include "../Node.h"
#include "../App/PhylogeneticsApplicationTools.h"
#include "../Io/Newick.h"
#include "../TreeIterator.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/Numeric/Random/RandomTools.h>

#include <iostream>
#include <fstream>
#include <algorithm>

using namespace bpp;

using namespace std;

#define STATE "state"

/******************************************************************************/

StochasticMapping::StochasticMapping(const TreeLikelihood* tl, size_t numOfMappings):
    mappingParameters_(),
    baseTree_(),
    tl_(),
    ConditionalProbabilities_(),
    nodesCounter_(0),
    numOfMappings_(numOfMappings)
{
    tl_ = tl;
    baseTree_ = tl_->getTree().clone();                      // this calls clone - but for some reason upson deletion a segnetation fault occurs
    giveNamesToInternalNodes(baseTree_);                     // set names for the internal nodes of the tree, in case of absence                       
    const SubstitutionModel* model = dynamic_cast<const SubstitutionModel*>(tl_->getModelForSite(0,0));
    mappingParameters_ = new SimpleMutationProcess(model);   // the procedure assumes that the same model applies to all the branches of the tree
    ComputeConditionals();                                               
}

/******************************************************************************/

StochasticMapping::~StochasticMapping()
{
    if (mappingParameters_) delete mappingParameters_;
    if (baseTree_) delete baseTree_; // delete the cloned base tree
}

/******************************************************************************/

void StochasticMapping::generateStochasticMapping(vector<Tree*>& mappings) 
{
 
    for (size_t i=0; i<numOfMappings_; ++i)
    {
        // clone the base tree to acheive the skeleton in which the mapping will be represented
        Tree* mapping = baseTree_->clone();
        setLeafsStates(mapping);
        
        /* step 2: simulate a set of ancestral states, based on the fractional likelihoods from step 1 */
        sampleAncestrals(mapping);

        /* step 3: simulate mutational history of each lineage of the phylogeny, conditional on the ancestral states */
        sampleMutationsGivenAncestrals(mapping);

        // add the mapping to the vector of mapping
        mappings.push_back(mapping);
    }
}

/******************************************************************************/

void StochasticMapping::setExpectedAncestrals(Tree* expectedMapping)
{
    TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(expectedMapping); 
    PreOrderTreeIterator* treeIt = new PreOrderTreeIterator(*ttree);
    for (Node* node = treeIt->begin(); node != treeIt->end(); node = treeIt->next()) {
        if (!node->isLeaf()) {
            int nodeId = node->getId();
            if (!node->hasFather())
            {
                size_t rootState = distance(ConditionalProbabilities_[nodeId][0].begin(), max_element(ConditionalProbabilities_[nodeId][0].begin(), ConditionalProbabilities_[nodeId][0].end()));
                setNodeState(node,rootState);
            }
            else
            {
                int fatherState = getNodeState(node->getFather());
                size_t sonState = distance(ConditionalProbabilities_[nodeId][fatherState].begin(), max_element(ConditionalProbabilities_[nodeId][fatherState].begin(), ConditionalProbabilities_[nodeId][fatherState].end()));
                setNodeState(node,sonState); // in the case of a leaf, the assigned state must be sampled
            }
        }
	}
    delete(treeIt);
}

/******************************************************************************/

Tree* StochasticMapping::generateExpectedMapping(vector<Tree*>& mappings, size_t divMethod)
{
    // initialize the expected history
    nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size()-1;
    Tree* expectedMapping = baseTree_->clone();
    setLeafsStates(expectedMapping);
    
    // set the ancestral states accrdonig to the maximal posterior (i.e, conditional) probability
    setExpectedAncestrals(expectedMapping);
    
    // update the expected history with the dwelling times
    size_t statesNum = tl_->getNumberOfStates();
    vector<Node*> nodes = (dynamic_cast<TreeTemplate<Node>*>(expectedMapping))->getNodes(); // get all the base nodes in post order ahead (as a post order iterator would also encounter nodes that are added during the procedure)
    for (size_t n=0; n< nodes.size(); ++n) {
        Node* node = nodes[n];
        if (node->hasFather()) // for any node except to the root
            {
            // initialize vector of average dwelling times for the branch stemming from node
            VDouble AverageDwellingTimes;
            AverageDwellingTimes.clear();
            AverageDwellingTimes.resize(statesNum,0);
            // compute the average dwelling times of all the states
            for (size_t i=0; i<mappings.size(); ++i)
            {
                TreeTemplate<Node>* mapping =  dynamic_cast<TreeTemplate<Node>*>(mappings[i]);
                // get the pointers to the node and its father in the i'th mapping
                Node* curNode = mapping->getNode(node->getName());
                Node* father = mapping->getNode(node->getFather()->getName()); // the original father of the node (according to the base tree) in the mapping
                while (curNode != father)
                {
                    AverageDwellingTimes[getNodeState(curNode)] += curNode->getDistanceToFather();
                    curNode = curNode->getFather();
                }
            }
            double branchLength = node->getDistanceToFather();   // this is the length of the original branch in the base tree
            bool updateBranch = true;
            for (size_t state=0; state<statesNum; ++state)
            {
                AverageDwellingTimes[state] /= static_cast<double>(mappings.size());
                if (AverageDwellingTimes[state] == branchLength) // if one of the dwelling times equals the branch length, then there is only one state along te branch and there is no need to edit it
                {
                    updateBranch = false;
                }  
            }
            // break the branch according to average dwelling times
            if (updateBranch)
            {
                updateBranchByDwellingTimes(node, AverageDwellingTimes, divMethod);
            }
        }
    }
    return expectedMapping;
}

/******************************************************************************/

int StochasticMapping::getNodeState(const Node* node) const
{
    return (dynamic_cast<const BppInteger*>(node->getNodeProperty(STATE)))->getValue();
}

/******************************************************************************/

void StochasticMapping::setNodeState(Node* node, size_t state)
{
    node->setNodeProperty(STATE,BppInteger(static_cast<int>(state)));
}

/******************************************************************************/

void StochasticMapping::giveNamesToInternalNodes(Tree* tree)
{
    TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(tree);
    vector<Node*> nodes = ttree->getNodes();
    for (size_t i=0; i<nodes.size(); ++i) {
        if (!nodes[i]->hasName())
            nodes[i]->setName("_baseInternal_" + TextTools::toString(nodes[i]->getId()));
    }  
}

/******************************************************************************/

void StochasticMapping::setLeafsStates(Tree* mapping)
{
    const SiteContainer* leafsStates = tl_->getData();
    TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(mapping); 
    vector<Node*> nodes = ttree->getNodes();

    for (size_t i=0; i<nodes.size(); ++i)
    {
        if (nodes[i]->isLeaf()) {
            string nodeName = nodes[i]->getName();
            size_t leafState = static_cast<int>(tl_->getAlphabetStateAsInt(leafsStates->getSequence(nodeName).getValue(0)));
            setNodeState(nodes[i],leafState);
        }
    }
}

/******************************************************************************/

void StochasticMapping::computeFractionals(VVDouble& fractionalProbabilities)
{
    // some auxiliiary variables
    size_t statesNum = tl_->getNumberOfStates();
    const TransitionModel* model = tl_->getModelForSite(0,0); // this calls assumes that all the sites and all the branches are assoiacted with the same node
    const SiteContainer* leafsStates = tl_->getData();
    TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(baseTree_); 
    vector<Node*> nodes = ttree->getNodes();

    // compute the fractional probabilities according to Felsenstein prunnig algorithm: for each node nodes[i] and state s compute: P(Data[leafs under node[i]]|node[i] has state s] 
    for (size_t i=0; i<nodes.size(); ++i) // traverse the tree in post-order
    {
        int nodeId = nodes[i]->getId();
        string nodeName = nodes[i]->getName();
        if (nodes[i]->isLeaf()) // if the node is a leaf, set the fractional probability of its state to 1, and the rest ot 0
        {
            size_t leafState = static_cast<int>(tl_->getAlphabetStateAsInt(leafsStates->getSequence(nodeName).getValue(0)));
            for (size_t nodeState=0; nodeState<statesNum; ++nodeState)
		    {
                if (nodeState != leafState)
                {
                    fractionalProbabilities[nodeId][nodeState] = 0;
                }
                else {
                    fractionalProbabilities[nodeId][nodeState] = 1; 
                }
            }   
        }
        else                // if the node is internal, follow the Felesenstein computation rule to compute the fractional probability
        {   
            for (size_t nodeState=0; nodeState<statesNum; ++nodeState)
		    {
                double fullProb = 1;
                for (size_t j=0; j<(nodes[i]->getNumberOfSons()); ++j) // for each son of the node, sum over the probabilities of all its assignments given its father's state (i.e, nodeState)
                {
                    double sonProb = 0;
                    double bl = nodes[i]->getSon(j)->getDistanceToFather();
                    for(size_t sonState=0; sonState<statesNum; ++sonState)
                    {
                        sonProb += model->Pij_t(nodeState, sonState, bl) * fractionalProbabilities[nodes[i]->getSon(j)->getId()][sonState];
                    }
                    fullProb *= sonProb;
                }
                fractionalProbabilities[nodeId][nodeState] = fullProb;
            }
		}
    }
}

/******************************************************************************/

void StochasticMapping::ComputeConditionals()
{
    // some auxiliiary variables
    size_t statesNum = tl_->getNumberOfStates();
    VDouble rootProbabilities = tl_->getRootFrequencies(0); 
    const TransitionModel* model = tl_->getModelForSite(0,0);
    TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(baseTree_); 
    vector<Node*> nodes = ttree->getNodes(); // here - how many nodes? only 4! what happened?

    /* compute the fractional probabilities for each node and state */
    VVDouble fractionalProbabilities;
    fractionalProbabilities.clear();
    fractionalProbabilities.resize(nodes.size(), VDouble(statesNum));
    computeFractionals(fractionalProbabilities);
    
    /*  compute the conditional probabilities: for each combination of nodes son, father, compute Pr(son recieves sonState | father has fatherState) */
    ConditionalProbabilities_.clear();
    ConditionalProbabilities_.resize((nodes.size()));
    
    for (size_t i=0; i<nodes.size(); ++i)
    {
        if (!nodes[i]->isLeaf() || !nodes[i]->hasFather()) { // the second condition will catch the root even if it has a single child (in which case, isLeaf() returns true)
            int nodeId = nodes[i]->getId();
            ConditionalProbabilities_[nodes[i]->getId()].resize(statesNum, VDouble(statesNum)); // use aux variable numOfStates
            if (!(nodes[i]->hasFather())) { // if the node is the root -> set the conditional probability to be same for all "fatherStates"
                double sum = 0.0;
                for (size_t sonState=0; sonState<statesNum; ++sonState)
                {
                    double stateConditionalNominator = fractionalProbabilities[nodeId][sonState] * rootProbabilities[sonState];
                    for (size_t fatherState=0; fatherState<statesNum; ++fatherState)
                        ConditionalProbabilities_[nodeId][fatherState][sonState] = stateConditionalNominator;
                    sum += stateConditionalNominator;
                }
                for (size_t fatherState=0; fatherState<statesNum; ++fatherState)
                {
                    for (size_t sonState=0; sonState<statesNum; ++sonState)
                        ConditionalProbabilities_[nodeId][fatherState][sonState] /= sum;
                }
            }
            else {                          // else -> follfow equation (10) from the paper to compute the consitional assingment probabilities given the ones of his father
                for (size_t fatherState=0; fatherState<statesNum; ++fatherState)
                {
                    double sum = 0.0; 
                    for (size_t sonState=0; sonState<statesNum; ++sonState)
                    {
                        double stateConditionalNominator = fractionalProbabilities[nodes[i]->getFatherId()][fatherState] * model->Pij_t(fatherState,sonState,nodes[i]->getDistanceToFather());
                        ConditionalProbabilities_[nodeId][fatherState][sonState] = stateConditionalNominator;
                        sum += stateConditionalNominator;
                    }
                    for (size_t sonState=0; sonState<statesNum; ++sonState)
                        ConditionalProbabilities_[nodeId][fatherState][sonState] /= sum;
                }
            }
        }
    }
}

/******************************************************************************/

size_t StochasticMapping::sampleState(const VDouble& distibution)
{
    size_t state=0;        // the default state is 0
    double prob = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);

    for (size_t i=0; i<distibution.size(); ++i) {
        prob-=distibution[i];
        if (prob < 0) { // if the the sampled probability is smaller than the probability to choose state i -> set state to be i
            state = i; 
            break;
        }
    }
    return state;
}

/******************************************************************************/

void StochasticMapping::sampleAncestrals(Tree* mapping)
{
    TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(mapping); 
    PreOrderTreeIterator* treeIt = new PreOrderTreeIterator(*ttree);
    for (Node* node = treeIt->begin(); node != treeIt->end(); node = treeIt->next()) {
        if (!node->isLeaf()) {
            int nodeId = node->getId();
            if (!node->hasFather())
            {
                size_t rootState = sampleState(ConditionalProbabilities_[nodeId][0]); // set father state to 0 (all the entries in the fatherState level are the same anyway)
                setNodeState(node,rootState);
            }
            else
            {
                int fatherState = getNodeState(node->getFather());
                size_t sonState = sampleState(ConditionalProbabilities_[nodeId][fatherState]);
                setNodeState(node,sonState); // in the case of a leaf, the assigned state must be sampled
            }
        }
	}
    delete(treeIt);
}

/******************************************************************************/

void StochasticMapping::sampleMutationsGivenAncestrals(Tree* mapping) 
{
    TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(mapping); 
    vector<Node*> nodes = ttree->getNodes();
    nodesCounter_ = nodes.size()-1; // intialize nodesCounter_ according ot the number of nodes in the base tree
    for (size_t i=0; i <nodes.size(); ++i)
    {
        
        Node* son = nodes[i];
        if (son->hasFather()) {
            sampleMutationsGivenAncestralsPerBranch(son);
        }
    }
}

/******************************************************************************/

void StochasticMapping::updateBranchMapping(Node* son, const MutationPath& branchMapping)
{  
    const vector<size_t> states = branchMapping.getStates();
    const VDouble times = branchMapping.getTimes();
    Node* curNode = son;
    Node* nextNode;
    int eventsNum = static_cast<int>(branchMapping.getNumberOfEvents());
    
    if (eventsNum == 0) // if there are no events >-return nothing
        return;
    else
    {
        for (int i=eventsNum-1; i>-1; --i) { // add a new node to represent the transition
            nodesCounter_ += 1;
            const string name = "_mappingInternal" + TextTools::toString(nodesCounter_) + "_";
            nextNode = new Node(static_cast<int>(nodesCounter_), name);
            setNodeState(nextNode,states[i]);
            nextNode->setDistanceToFather(times[i]);
            
            // set the father to no longer be the father of curNode
            Node* originalFather = curNode->getFather();
            originalFather->removeSon(curNode); // also removes originalFather as the father of curNode
            // set nextNode to be the new father of curNode
            curNode->setFather(nextNode); // also adds curNode to the sons of nextNode
            // set curNode's original father ot be the father of nextNode
            nextNode->setFather(originalFather); // also adds nextNode to the sons of originalFather - or not? make sure this is the father at all times
            curNode = nextNode;
        }
        return;
    }
}

/******************************************************************************/

void StochasticMapping::sampleMutationsGivenAncestralsPerBranch(Node* son, size_t maxIterNum)
{
    Node* father = son->getFather();
    double branchLength = son->getDistanceToFather();
    size_t fatherState = getNodeState(father); 
    size_t sonState = getNodeState(son);

    /* simulate mapping on a branch until you manage to finish at the son's state */
	for (size_t i = 0; i < maxIterNum; ++i) 
	{
		double disFromNode = 0.0;
		size_t curState = fatherState;
        MutationPath tryMapping(mappingParameters_->getSubstitutionModel()->getAlphabet(), fatherState, branchLength);

		double timeTillChange;
        // if the father's state is not the same as the son's state -> use the correction corresponding to equation (11) in the paper
        if (fatherState != sonState)
		{   // sample timeTillChange conditional on it being smaller than branchLength
			double u = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
            double waitingTimeParam = -1 * mappingParameters_->getSubstitutionModel()->Qij(fatherState, fatherState); // get the parameter for the exoponential distribution to draw the waiting time from
            double tmp = u * (1.0 - exp(branchLength * -waitingTimeParam));
            timeTillChange =  -log(1.0 - tmp) / waitingTimeParam;
            assert (timeTillChange < branchLength);
		} else {
            timeTillChange = mappingParameters_->getTimeBeforeNextMutationEvent(fatherState); // draw the time until a transition from exponential distribution with the rate of leaving fatherState   
        }
        
		while (disFromNode + timeTillChange < branchLength)  // a jump occured but not passed the whole branch -> 
		{
            tryMapping.addEvent(curState, timeTillChange);                                        // add the current state and time to branch history
			
            disFromNode += timeTillChange;			
            timeTillChange = mappingParameters_->getTimeBeforeNextMutationEvent(curState);        // draw the time until a transition from exponential distribution with the rate of leaving curState
       	    curState = mappingParameters_->mutate(curState);                                      // draw the state to transition to after from initial state curState based on the relative tranistion rates distribution (see MutationProcess.cpp line 50)	
        }
		// the last jump passed the length of the branch -> finish the simulation and check if it's sucessfull (i.e, mapping is finished at the son's state)
		if (curState != sonState) // if the simulation failed, try again
		{      
			continue;
		}
		else                      // if the simulation was sucessfully, add it to the build mapping
		{
			double timeOfJump = branchLength - disFromNode;
            son->setDistanceToFather(timeOfJump);
            updateBranchMapping(son, tryMapping);     // add the successfull simulation to the build mapping
            return;
		}
	}
	// if all simulations failed -> throw an exception
    throw Exception("could not produce simulations with father = " + TextTools::toString(fatherState) + " son " + TextTools::toString(sonState) + " branch length = " + TextTools::toString(branchLength));
}

/******************************************************************************/

void StochasticMapping::updateBranchByDwellingTimes(Node* node, VDouble& dwellingTimes, size_t divMethod)
{
    /* first, convert the dwelling times vector to a mutation path of the branch */
    size_t statesNum = tl_->getNumberOfStates();
    size_t sonState = getNodeState(node);
    size_t fatherState = getNodeState(node->getFather());
    double branchLength = node->getDistanceToFather();
    double Pf = 1;
    double Ps = 1;
    MutationPath branchMapping(mappingParameters_->getSubstitutionModel()->getAlphabet(), fatherState, branchLength);
    // set the first event with the dwelling time that matches the state of the father
    if (fatherState == sonState)
    {
        if (divMethod == 0) // here - do what Itay suggested
        {
            int ancestorState = 0;
            if (node->getFather()->hasFather())
            {
                ancestorState = getNodeState(node->getFather()->getFather());
            }
            if (node->hasFather())
            {
                Pf = ConditionalProbabilities_[node->getFather()->getId()][ancestorState][fatherState]; // segmentation fault on node id 30
            }
            Pf = Ps = 1;
            if (!node->isLeaf())
            {
                Ps = ConditionalProbabilities_[node->getId()][fatherState][sonState];
            }
            branchMapping.addEvent(fatherState, dwellingTimes[fatherState]*(Pf/(Pf+Ps)));
        }
        else
        {
            branchMapping.addEvent(fatherState, 0);
        }
    }
    else
    {
         branchMapping.addEvent(fatherState, dwellingTimes[fatherState]);
    }
    // set all events except for the one entering the son
    for (size_t state=0; state<statesNum; ++state)
    {
        if (state != fatherState && state != sonState && dwellingTimes[state] > 0) // if the state matches an event which is not the first or the last -> add it
        {
            branchMapping.addEvent(state, dwellingTimes[state]);
        }  
    }
    // change the length of the branch whose bottom node is the son according to the dwelling time of the relevant state
    if (fatherState == sonState)
    {
        if (divMethod == 0)
        {
            node->setDistanceToFather(dwellingTimes[sonState]*(Ps/(Pf+Ps)));
        }
        else
        {
            node->setDistanceToFather(dwellingTimes[sonState]);
        }
    }
    else
    {
        node->setDistanceToFather(dwellingTimes[sonState]);
    }
    
    /* secondly, update the expected history with the dwelling times-based mutation path */
    updateBranchMapping(node, branchMapping); 
}