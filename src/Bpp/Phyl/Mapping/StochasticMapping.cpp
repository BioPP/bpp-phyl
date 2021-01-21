#include "StochasticMapping.h"
#include "../Likelihood/TreeLikelihoodTools.h"
#include "../Likelihood/DRHomogeneousTreeLikelihood.h"
#include "../Likelihood/TreeLikelihood.h"
#include "../Simulation/MutationProcess.h"
#include "../Node.h"
#include "../App/PhylogeneticsApplicationTools.h"
#include "../Io/Newick.h"
#include "../TreeIterator.h"
#include "RewardMappingTools.h"
#include "Reward.h"
#include "DecompositionReward.h"
#include "ProbabilisticRewardMapping.h"
#include "../Model/RateDistribution/ConstantRateDistribution.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>
#include <Bpp/Seq/Alphabet/NumericAlphabet.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric> // to sum over items in a vector

using namespace bpp;
using namespace std;

#define STATE "state"

/******************************************************************************/

StochasticMapping::StochasticMapping(const TreeLikelihood* tl, size_t numOfMappings) :
  mappingParameters_(),
  baseTree_(),
  tl_(),
  fractionalProbabilities_(),
  ConditionalProbabilities_(),
  nodesCounter_(0),
  numOfMappings_(numOfMappings),
  nodeIdToIndex_()
{

  tl_ = tl;
  baseTree_ = tl_->getTree().clone();                      // this calls clone - but for some reason upson deletion a segnetation fault occurs
  vector<Node*> nodes = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes();
  giveNamesToInternalNodes(baseTree_);                     // set names for the internal nodes of the tree, in case of absence
  for (size_t i=0; i<nodes.size(); ++i)
  {
    nodeIdToIndex_[nodes[i]->getId()] = i; 
  }
  const SubstitutionModel* model = dynamic_cast<const SubstitutionModel*>(tl_->getModelForSite(0, 0));
  mappingParameters_ = new SimpleMutationProcess(model);   // the procedure assumes that the same model applies to all the branches of the tree
  ComputeConditionals();
}

/******************************************************************************/

StochasticMapping::~StochasticMapping()
{
  if (mappingParameters_)
    delete mappingParameters_;
  if (baseTree_)
    delete baseTree_; // delete the cloned base tree
}

/******************************************************************************/

void StochasticMapping::generateStochasticMapping(vector<Tree*>& mappings)
{
  for (size_t i = 0; i < numOfMappings_; ++i)
  {
    // clone the base tree to acheive the skeleton in which the mapping will be represented
    Tree* mapping = baseTree_->clone();
    map<int,vector<size_t>> leafIdToStates = setLeafsStates(mapping);

    /* step 2: simulate a set of ancestral states, based on the fractional likelihoods from step 1 */
    sampleAncestrals(mapping, leafIdToStates);

    /* step 3: simulate mutational history of each lineage of the phylogeny, conditional on the ancestral states */
    sampleMutationsGivenAncestrals(mapping);

    // add the mapping to the vector of mapping
    mappings.push_back(mapping);

    // reset the nodes counter
    nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;
  }
}

/******************************************************************************/

void StochasticMapping::setExpectedAncestrals(Tree* expectedMapping, VVDouble& ancestralStatesFrequencies)
{
  TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(expectedMapping);
  vector<Node*> nodes = ttree->getNodes();
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    Node* node = nodes[i];
    auto d = distance(ancestralStatesFrequencies[nodeIdToIndex_[node->getId()]].begin(), max_element(ancestralStatesFrequencies[nodeIdToIndex_[node->getId()]].begin(), ancestralStatesFrequencies[nodeIdToIndex_[node->getId()]].end()));
    size_t state = static_cast<size_t>(d);
    setNodeState(node, state);
  }
}

/******************************************************************************/

Tree* StochasticMapping::generateExpectedMapping(const vector<Tree*>& mappings, size_t divMethod)
{
  // initialize the expected history
  nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;
  Tree* expectedMapping = baseTree_->clone();
  map<int,vector<size_t>> leafIdToStates = setLeafsStates(expectedMapping);

  // compute a vector of the posterior asssignment probabilities for each inner node
  VVDouble ancestralStatesFrequencies;
  // initialize the vector
  ancestralStatesFrequencies.clear();
  vector<Node*> nodes = dynamic_cast<TreeTemplate<Node>*>(expectedMapping)->getNodes();
  size_t statesNum = tl_->getNumberOfStates();
  ancestralStatesFrequencies.resize(nodes.size(), VDouble(statesNum));
  
  computeStatesFrequencies(ancestralStatesFrequencies, mappings);

  // set the ancestral states accrdonig to the maximal posterior (i.e, conditional) probability
  setExpectedAncestrals(expectedMapping, ancestralStatesFrequencies);

  // update the expected history with the dwelling times
  for (size_t n = 0; n < nodes.size(); ++n)
  {
    Node* node = nodes[n];
    if (node->hasFather()) // for any node except to the root
    {
      // initialize vector of average dwelling times for the branch stemming from node
      VDouble AverageDwellingTimes;
      AverageDwellingTimes.clear();
      AverageDwellingTimes.resize(statesNum, 0);
      // compute the average dwelling times of all the states
      for (size_t i = 0; i < mappings.size(); ++i)
      {
        const TreeTemplate<Node>* mapping =  dynamic_cast<const TreeTemplate<Node>*>(mappings[i]);
        // get the pointers to the node and its father in the i'th mapping
        const Node* curNode = mapping->getNode(node->getName());
        const Node* father = mapping->getNode(node->getFather()->getName()); // the original father of the node (according to the base tree) in the mapping
        while (curNode != father)
        {
          AverageDwellingTimes[getNodeState(curNode)] += curNode->getDistanceToFather(); // only model states, which are non-negative, are considered
          curNode = curNode->getFather();
        }
      }
      double branchLength = node->getDistanceToFather();   // this is the length of the original branch in the base tree
      bool updateBranch = true;
      for (size_t state = 0; state < statesNum; ++state)
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
        updateBranchByDwellingTimes(node, AverageDwellingTimes, ancestralStatesFrequencies, divMethod);
      }
    }
  }
  nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;
  return expectedMapping;
}

/******************************************************************************/

Tree* StochasticMapping::generateAnalyticExpectedMapping(size_t divMethod)
{
  /* Compute the posterior assignment probabilities to internal nodes, based on the fractional probablities computed earlier */
  size_t statesNum = tl_->getNumberOfStates();
  vector<int> nodeIds = baseTree_->getNodesId();
  int nodeId;
  VVDouble posteriorProbabilities;
  posteriorProbabilities.clear();
  posteriorProbabilities.resize(baseTree_->getNumberOfNodes(), VDouble(statesNum));
  double nodeDataProb;
  // because the sum of partial likelihoods (i.e, the fractional probabilities) is in fact the probablity of the data, it is sufficient to standardize the vector of fractional probabilires for each node to obtain the posterior probabilities
  for (size_t n = 0; n < baseTree_->getNumberOfNodes(); ++n)
  {
    nodeId = nodeIds[n];
    nodeDataProb = 0;
    for (size_t s = 0; s < statesNum; ++s)
    {
      nodeDataProb = nodeDataProb + fractionalProbabilities_[nodeIdToIndex_[nodeId]][s];
    }
    for (size_t nodeState = 0; nodeState < statesNum; ++nodeState)
    {
      posteriorProbabilities[nodeIdToIndex_[nodeId]][nodeState] = fractionalProbabilities_[nodeIdToIndex_[nodeId]][nodeState] / nodeDataProb;
    }
  }

  /* Assign states to internal nodes based on the majority rule over the posterior probabilities */
  Tree* expectedMapping = baseTree_->clone();
  setLeafsStates(expectedMapping);
  setExpectedAncestrals(expectedMapping, posteriorProbabilities);

  /* Compute the reward per state per site - expect two entries per site (that is, two entries in total).
     Let r0 be the reward of state 0 nd r1 the reward of state 1. */
  // create a numeric alphabet whose states correspond to the model states
  UserAlphabetIndex1* alpha = new UserAlphabetIndex1(tl_->getAlphabet());
  DiscreteDistribution* rDist = new ConstantRateDistribution();
  TransitionModel* model = tl_->getModelForSite(0, 0)->clone();
  if (dynamic_cast<MarkovModulatedSubstitutionModel*>(model) != nullptr)
  {
    ReversibleSubstitutionModel* nestedModel = (dynamic_cast<MarkovModulatedSubstitutionModel*>(model))->getNestedModel()->clone(); // here, use the nested model with which the rewards methods is equiped to deal with
    delete model;
    model = dynamic_cast<TransitionModel*>(nestedModel);
  }

  DRTreeLikelihood* drtl = new DRHomogeneousTreeLikelihood(*baseTree_, *(tl_->getData()), model, rDist, false);
  drtl->initialize();
  vector<int> ids = baseTree_->getNodesId();

  /* Compute the expected dwelling times per branch and state as follows:
     For branch b of length t, the average welling time in state 0 is r0*t (based on Minin and Suchard paper).
     The average dwelling time in state 1 should complement to t (make sure of it!) */
  vector<Node*> nodes = dynamic_cast<TreeTemplate<Node>*>(expectedMapping)->getNodes();
  double branchLength;
  Node* node;
  statesNum = tl_->getNumberOfStates();
  VVDouble expectedDwellingTimes;
  expectedDwellingTimes.clear();
  expectedDwellingTimes.resize(nodes.size(), VDouble(statesNum));
  map <size_t, int> alphabetStatesToModelStates; 
  size_t alphabetStatesNum = tl_->getAlphabet()->getNumberOfStates();
  vector<string> resolvedStates = tl_->getAlphabet()->getResolvedChars(); // here, only treat resolved states of the alphabet
  //for (size_t s = 0; s<alphabetStatesNum; ++s)
  for (size_t s = 0; s<resolvedStates.size(); ++s)
  {
    //int character_state_a = tl_->getAlphabet()->getStateAt(s).getNum(); // this state could represent a gap or an unknown character directly and not correspond to the set of mappable states of the alphabet
    int character_state_a = tl_->getAlphabet()->getState(resolvedStates[s]).getNum();
	// special case for gaps - treat as unknown (i.e., the last state in the alphabet)    
    if (character_state_a < 0) // this section should not be visited if only resolved states are regarded
      character_state_a = tl_->getAlphabet()->getStateAt(alphabetStatesNum-1).getNum();
    vector<int> alphabetCorrespondingStates_a = tl_->getAlphabet()->getAlias(character_state_a);
    for (size_t cs=0; cs<alphabetCorrespondingStates_a.size(); ++cs)
    {
      int state_a = alphabetCorrespondingStates_a[cs];
      alpha->setIndex(state_a, 1); // set the reward of the state as 1 and the reward for the rest of the states as 0
      //for (size_t m = 0; m < alphabetStatesNum; ++m)
      for (size_t m=0; m<resolvedStates.size(); ++m)
	  {
        //int character_state_b = tl_->getAlphabet()->getStateAt(m).getNum(); // this state could represent a gap or an unknown character directly and not correspond to the set of mappable states of the alphabet
        int character_state_b = tl_->getAlphabet()->getState(resolvedStates[m]).getNum();
		// special case for gaps - treat as unknown (i.e., the last state in the alphabet)
        if (character_state_b < 0)
          character_state_b = tl_->getAlphabet()->getStateAt(alphabetStatesNum-1).getNum();
        vector<int> alphabetCorrespondingStates_b = tl_->getAlphabet()->getAlias(character_state_b);
        for (size_t cm=0; cm<alphabetCorrespondingStates_b.size(); ++cm)
        {
            int state_b = alphabetCorrespondingStates_b[cm];
            bool in_a = false;
            for (size_t i=0; i<alphabetCorrespondingStates_a.size(); ++i)
            {
              if (state_b == alphabetCorrespondingStates_a[i])
                in_a = true;
            }
            if (!in_a)
            {
              alpha->setIndex(state_b, 0);
            }
        }
      }
	  unique_ptr<Reward> reward(new DecompositionReward(dynamic_cast<const SubstitutionModel*>(model), alpha));
      unique_ptr<ProbabilisticRewardMapping> mapping(RewardMappingTools::computeRewardVectors(*drtl, ids, *reward, false));
      
	  for (size_t n = 0; n < nodes.size(); ++n)
      {
		node = nodes[n];
		if (node->hasFather()) // for any node except to the root
        {
          double dwellingTime =  mapping->getReward(node->getId(), 0);
          vector <size_t> correspondingModelStates = model->getModelStates(state_a);
          for (size_t ms=0; ms<correspondingModelStates.size(); ++ms)
            expectedDwellingTimes[nodeIdToIndex_[nodes[n]->getId()]][correspondingModelStates[ms]] = dwellingTime / static_cast<double>(correspondingModelStates.size());
		}
      }
    }
  }

  /* free the resources - all od these are deleted via reward and mapping variables*/
  delete alpha;
  delete rDist;
  delete model;
  delete drtl;


  // standardize expected dwelling itmes, if needed, and update the mapping accorgingly
  double sumOfDwellingTimes;
  bool updateBranch;
  nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;
  for (size_t n = 0; n < nodes.size(); ++n)
  {
    node = nodes[n];
    if (node->hasFather()) // for any node except to the root
    {
      branchLength = node->getDistanceToFather();
      sumOfDwellingTimes = 0;
      updateBranch = true;
      for (size_t s = 0; s < statesNum; ++s)
      {
        if (expectedDwellingTimes[nodeIdToIndex_[node->getId()]][s] == 0)

        {
          updateBranch = false;
        }
        sumOfDwellingTimes = sumOfDwellingTimes + expectedDwellingTimes[nodeIdToIndex_[node->getId()]][s];
      }

      if (branchLength < 0.00001) // branch length is 0 -> no need to update mapping on the branch
      {
        node->setDistanceToFather(branchLength);
        updateBranch = false;
      }
      else
      {
        if (sumOfDwellingTimes != branchLength)
        {
          for (size_t s = 0; s < statesNum; ++s)
          {
            expectedDwellingTimes[nodeIdToIndex_[node->getId()]][s] =  branchLength * (expectedDwellingTimes[nodeIdToIndex_[node->getId()]][s]) / sumOfDwellingTimes;
          }
        }
      }

      if (updateBranch)
      {
        updateBranchByDwellingTimes(node, expectedDwellingTimes[nodeIdToIndex_[node->getId()]], posteriorProbabilities, divMethod);
      }
    }
  }
  nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;

  return expectedMapping;
}
/******************************************************************************/

size_t StochasticMapping::getNodeState(const Node* node)
{
  return static_cast<size_t>((dynamic_cast<const BppInteger*>(node->getNodeProperty(STATE)))->getValue());
}

/******************************************************************************/

void StochasticMapping::setNodeState(Node* node, size_t state)
{
  BppInteger* stateProperty = new BppInteger(static_cast<int>(state));
  node->setNodeProperty(STATE, *stateProperty);
  delete stateProperty;
}

/******************************************************************************/

void StochasticMapping::giveNamesToInternalNodes(Tree* tree)
{
  TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(tree);
  vector<Node*> nodes = ttree->getNodes();
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    if (!nodes[i]->hasName())
      nodes[i]->setName("_baseInternal_" + TextTools::toString(nodes[i]->getId()));
  }
}

/******************************************************************************/

vector<size_t> StochasticMapping::getLeafModelStates(Node* node)
{
  const SiteContainer* leafsStates = tl_->getData();
  string nodeName = node->getName();
  int leafAlphabetState = leafsStates->getSequence(nodeName).getValue(0);
  vector<int> leafCharacterStates = leafsStates->getAlphabet()->getAlias(leafAlphabetState); 
  const SubstitutionModel* model = dynamic_cast<const SubstitutionModel*>(tl_->getModelForSite(0, 0)); // this call assumes that all the sites and all the branches are assoiacted with the same node
  vector<size_t> leafModelStates;
  leafModelStates.clear();
  for (size_t s=0; s<leafCharacterStates.size(); ++s)
  {
    vector <size_t> correspondingModelStates = model->getModelStates(leafCharacterStates[s]);
    for (size_t m=0; m<correspondingModelStates.size(); ++m)
    {
      leafModelStates.push_back(correspondingModelStates[m]);
    }
  }
  return leafModelStates;
}
/******************************************************************************/

map<int,vector<size_t>> StochasticMapping::setLeafsStates(Tree* mapping)
{
  TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(mapping);
  vector<Node*> nodes = ttree->getNodes();
  map<int,vector<size_t>> leafIdToStates;
  for (auto node: nodes)
  {
    if (node->isLeaf())
    {
      leafIdToStates[node->getId()] = getLeafModelStates(node);
    }
  }
  return leafIdToStates;
}

/******************************************************************************/


void StochasticMapping::computeFractionals()
{
  // some auxiliiary variables
  size_t statesNum = tl_->getNumberOfStates();
  const TransitionModel* model = tl_->getModelForSite(0, 0); // this calls assumes that all the sites and all the branches are assoiacted with the same node
  TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(baseTree_);
  vector<Node*> nodes = ttree->getNodes();

  // compute the fractional probabilities according to Felsenstein prunnig algorithm: for each node nodes[i] and state s compute: P(Data[leafs under node[i]]|node[i] has state s]
  for (size_t i = 0; i < nodes.size(); ++i) // traverse the tree in post-order
  {
    size_t nodeIndex = nodeIdToIndex_[nodes[i]->getId()];
    string nodeName = nodes[i]->getName();
    if (nodes[i]->isLeaf()) // if the node is a leaf, set the fractional probability of its state to 1, and the rest ot 0
    {
      vector<size_t> leafModelStates = getLeafModelStates(nodes[i]);
      for (size_t s = 0; s < leafModelStates.size(); ++s)
      {
        fractionalProbabilities_[nodeIndex][leafModelStates[s]] = 1;
      }
    }
    else                // if the node is internal, follow the Felesenstein computation rule to compute the fractional probability
    {
      for (size_t nodeState = 0; nodeState < statesNum; ++nodeState)
      {
        double fullProb = 1;
        for (size_t j = 0; j < (nodes[i]->getNumberOfSons()); ++j) // for each son of the node, sum over the probabilities of all its assignments given its father's state (i.e, nodeState)
        {
          double sonProb = 0;
          double bl = nodes[i]->getSon(j)->getDistanceToFather();
          for (size_t sonState = 0; sonState < statesNum; ++sonState)
          {
            sonProb += model->Pij_t(nodeState, sonState, bl) * fractionalProbabilities_[nodeIdToIndex_[nodes[i]->getSon(j)->getId()]][sonState];
          }
          fullProb *= sonProb;
        }
        fractionalProbabilities_[nodeIndex][nodeState] = fullProb;
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
  const TransitionModel* model = tl_->getModelForSite(0, 0);
  TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(baseTree_);
  vector<Node*> nodes = ttree->getNodes(); // here - how many nodes? only 4! what happened?

  /* compute the fractional probabilities for each node and state */
  fractionalProbabilities_.clear();
  fractionalProbabilities_.resize(nodes.size(), VDouble(statesNum));
  computeFractionals();

  /*  compute the conditional probabilities: for each combination of nodes son, father, compute Pr(son recieves sonState | father has fatherState) */
  ConditionalProbabilities_.clear();
  ConditionalProbabilities_.resize(nodes.size(), VVDouble(statesNum, VDouble(statesNum)));
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    size_t nodeIndex = nodeIdToIndex_[nodes[i]->getId()];
    if (!nodes[i]->isLeaf() || !nodes[i]->hasFather())  // the second condition will catch the root even if it has a single child (in which case, isLeaf() returns true)
    {
      if (!(nodes[i]->hasFather()))  // if the node is the root -> set the conditional probability to be same for all "fatherStates"
      {
        double sum = 0.0;
        for (size_t sonState = 0; sonState < statesNum; ++sonState)
        {
          double stateConditionalNominator = fractionalProbabilities_[nodeIndex][sonState] * rootProbabilities[sonState];
          for (size_t fatherState = 0; fatherState < statesNum; ++fatherState)
          {
            ConditionalProbabilities_[nodeIndex][fatherState][sonState] = stateConditionalNominator;
          }
          sum += stateConditionalNominator;
        }
        for (size_t fatherState = 0; fatherState < statesNum; ++fatherState)
        {
          for (size_t sonState = 0; sonState < statesNum; ++sonState)
          {
            ConditionalProbabilities_[nodeIndex][fatherState][sonState] /= sum;
          }
        }
      }
      else                            // else -> follfow equation (10) from the paper to compute the consitional assingment probabilities given the ones of his father
      {
        for (size_t fatherState = 0; fatherState < statesNum; ++fatherState)
        {
          double sum = 0.0;
          for (size_t sonState = 0; sonState < statesNum; ++sonState)
          {
            double stateConditionalNominator = fractionalProbabilities_[nodeIdToIndex_[nodes[i]->getId()]][sonState] * model->Pij_t(fatherState, sonState, nodes[i]->getDistanceToFather());
            ConditionalProbabilities_[nodeIndex][fatherState][sonState] = stateConditionalNominator;
            sum += stateConditionalNominator;
          }
          for (size_t sonState = 0; sonState < statesNum; ++sonState)
          {
            ConditionalProbabilities_[nodeIndex][fatherState][sonState] /= sum;
          }
        }
      }
    }
    else // the node is a leaf, so the conditional probabilities should correspond to the possible leaf states assignments
    {
      vector<size_t> leafModelStates = getLeafModelStates(nodes[i]);  // gap integer is negative and cannot be cast to size_t instance 
      for (size_t fatherState = 0; fatherState < statesNum; ++fatherState)
      {
        double sum = 0.0;
        for (size_t s = 0; s < leafModelStates.size(); ++s)
        {
          size_t sonState = leafModelStates[s];
          double stateConditionalNominator = fractionalProbabilities_[nodeIdToIndex_[nodes[i]->getId()]][sonState] * model->Pij_t(fatherState, sonState, nodes[i]->getDistanceToFather());
          ConditionalProbabilities_[nodeIndex][fatherState][sonState] = stateConditionalNominator;
          sum += stateConditionalNominator;
        }
        for (size_t s = 0; s < leafModelStates.size(); ++s)
        {
          ConditionalProbabilities_[nodeIndex][fatherState][leafModelStates[s]] /= sum;
        }
      }
    }
  }
  
}

/******************************************************************************/

void StochasticMapping::computeStatesFrequencies(VVDouble& ancestralStatesFrequencies, const vector<Tree*>& mappings)
{
  // initialize the vector
  vector<int> nodeIds = baseTree_->getNodesId();
  size_t statesNum = tl_->getNumberOfStates();

  // compute the node assignment probabilities based on their frequency in the mappings
  for (size_t i = 0; i < nodeIds.size(); ++i)
  {
    size_t nodeIndex = nodeIdToIndex_[nodeIds[i]];
    string nodeName = baseTree_->getNodeName(nodeIds[i]);
    // go over all the mappings and collect the number of states assignment per node (inclusing leafs which can have ambiguous states like N, in which case they will not be consistent across mappings)
    fill(ancestralStatesFrequencies[nodeIndex].begin(), ancestralStatesFrequencies[nodeIndex].end(), 0); // reset all the values to 0
    for (size_t h = 0; h < mappings.size(); ++h)
    {
      Node* nodeInMapping = dynamic_cast<TreeTemplate<Node>*>(mappings[h])->getNode(nodeName); // the ID of this node is either the original node ID or it is a new ID which is necessarily positive
      ancestralStatesFrequencies[nodeIndex][getNodeState(nodeInMapping)]++; // the node index is non-negative value mapped to the node id, which could be negative
    }
    // now divide the vector entries by the number of mappings
    for (size_t nodeState = 0; nodeState < statesNum; ++nodeState)
    {
      ancestralStatesFrequencies[nodeIndex][nodeState] = ancestralStatesFrequencies[nodeIndex][nodeState] / static_cast<int>(mappings.size());
    }
  }
}


/******************************************************************************/

size_t StochasticMapping::sampleState(const VDouble& distibution)
{
  size_t state = 0;        // the default state is 0
  double prob = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);

  for (size_t i = 0; i < distibution.size(); ++i)
  {
    prob -= distibution[i];
    if (prob < 0)  // if the the sampled probability is smaller than the probability to choose state i -> set state to be i
    {
      state = i;
      break;
    }
  }
  return state;
}

/******************************************************************************/

void StochasticMapping::sampleAncestrals(Tree* mapping, map<int,vector<size_t>> leafIdToStates)
{
  TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(mapping);
  PreOrderTreeIterator* treeIt = new PreOrderTreeIterator(*ttree);
  for (const Node* node = treeIt->begin(); node != treeIt->end(); node = treeIt->next())
  {
    size_t nodeIndex = nodeIdToIndex_[node->getId()];
    if (!node->isLeaf())
    {
      if (!node->hasFather())
      {
        size_t rootState = sampleState(ConditionalProbabilities_[nodeIndex][0]); // set father state to 0 (all the entries in the fatherState level are the same anyway)
        setNodeState(ttree->getNode(node->getId()), rootState);
      }
      else
      {
        size_t fatherState = getNodeState(node->getFather());
        size_t sonState = sampleState(ConditionalProbabilities_[nodeIndex][fatherState]);
        setNodeState(ttree->getNode(node->getId()), sonState); 
      }
    }
    else
    {
      vector<size_t> leafStates = leafIdToStates[node->getId()];
      if (leafStates.size() == 1)
      {
        setNodeState(ttree->getNode(node->getId()), leafStates[0]);
      }
      else // in the case of a leaf wih multiple possible assignments, the assigned state must be sampled
      {
        size_t fatherState = getNodeState(node->getFather());
        size_t sonState = sampleState(ConditionalProbabilities_[nodeIndex][fatherState]);
        setNodeState(ttree->getNode(node->getId()), sonState);
      }
    }
  }
  delete treeIt;
}

/******************************************************************************/

void StochasticMapping::sampleMutationsGivenAncestrals(Tree* mapping)
{
  TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(mapping);
  vector<Node*> nodes = ttree->getNodes();
  nodesCounter_ = nodes.size() - 1; // intialize nodesCounter_ according ot the number of nodes in the base tree
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    Node* son = nodes[i];
    if (son->hasFather())
    {
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
  size_t eventsNum = branchMapping.getNumberOfEvents();

  if (eventsNum == 0) // if there are no events >-return nothing
    return;
  else
  {
    for (size_t i = eventsNum; i > 0; --i) // add a new node to represent the transition
    {
      nodesCounter_ = nodesCounter_ + 1;
      const string name = "_mappingInternal" + TextTools::toString(nodesCounter_) + "_";
      nextNode = new Node(static_cast<int>(nodesCounter_), name);
      setNodeState(nextNode, states[i - 1]);
      nextNode->setDistanceToFather(times[i - 1]);

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
    }
    else
    {
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

void StochasticMapping::updateBranchByDwellingTimes(Node* node, VDouble& dwellingTimes, VVDouble& ancestralStatesFrequencies, size_t divMethod)
{
  
  // first, convert the dwelling times vector to a mutation path of the branch
  size_t statesNum = tl_->getNumberOfStates();
  size_t sonState = getNodeState(node);
  size_t fatherState = getNodeState(node->getFather());
  double branchLength = node->getDistanceToFather();
  double Pf = 1;
  double Ps = 1;
  double shareOfFather = 0;
  double shareOfSon = 0;
  MutationPath* branchMapping = new MutationPath(mappingParameters_->getSubstitutionModel()->getAlphabet(), fatherState, branchLength);
  // set the first event with the dwelling time that matches the state of the father
  if (fatherState == sonState)
  {
    if (divMethod == 0)
    {
      if (node->hasFather())
      {
        Pf = ancestralStatesFrequencies[nodeIdToIndex_[node->getFather()->getId()]][fatherState];
      }
      Ps = ancestralStatesFrequencies[nodeIdToIndex_[node->getId()]][sonState];
      shareOfFather = Pf / (Pf + Ps);
      branchMapping->addEvent(fatherState, dwellingTimes[fatherState] * shareOfFather);
    }
    else
    {
      branchMapping->addEvent(fatherState, 0);
    }
  }
  else
  {
    branchMapping->addEvent(fatherState, dwellingTimes[fatherState]);
  }
  // set all events except for the one entering the son
  for (size_t state = 0; state < statesNum; ++state)
  {
    if (state != fatherState && state != sonState && dwellingTimes[state] > 0) // if the state matches an event which is not the first or the last -> add it
    {
      branchMapping->addEvent(state, dwellingTimes[state]);
    }
  }
  // change the length of the branch whose bottom node is the son according to the dwelling time of the relevant state
  if (fatherState == sonState)
  {
    if (divMethod == 0)
    {
      shareOfSon = 1 - shareOfFather;
      node->setDistanceToFather(dwellingTimes[sonState] * shareOfSon);
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
  updateBranchMapping(node, *branchMapping);
  delete branchMapping;
}
