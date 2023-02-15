//
// File: StochasticMapping.cpp
// Authors:
//

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Seq/Alphabet/NumericAlphabet.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>
#include <Bpp/Text/TextTools.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric> // to sum over items in a vector

#include "../Simulation/MutationProcess.h"
#include "DecompositionReward.h"
#include "ProbabilisticRewardMapping.h"
#include "Reward.h"
#include "RewardMappingTools.h"
#include "StochasticMapping.h"

using namespace bpp;
using namespace std;

#define STATE "state"

/******************************************************************************/

StochasticMapping::StochasticMapping(std::shared_ptr<LikelihoodCalculationSingleProcess> drl, size_t numOfMappings) :
  likelihood_(drl),
  tree_(make_shared<PhyloTree>(drl->substitutionProcess().parametrizablePhyloTree())),
//  mappingParameters_(drl->getSubstitutionProcess()),
  fractionalProbabilities_(),
  ConditionalProbabilities_(),
  nodesCounter_(0),
  numOfMappings_(numOfMappings)// ,
  // nodeIdToIndex_()
{
  giveNamesToInternalNodes(*tree_);                     // set names for the internal nodes of the tree, in case of absence
  ComputeConditionals();
}

/******************************************************************************/

StochasticMapping::~StochasticMapping()
{}

/******************************************************************************/

void StochasticMapping::generateStochasticMapping(vector<shared_ptr<PhyloTree> >& mappings)
{
  for (size_t i = 0; i < numOfMappings_; ++i)
  {
//    auto mapping = std::shared_ptr<PhyloTree>(*baseTree);
//   setLeafsStates(mapping);

    //   /* step 2: simulate a set of ancestral states, based on the fractional likelihoods from step 1 */
    //   sampleAncestrals(mapping);

    //   /* step 3: simulate mutational history of each lineage of the phylogeny, conditional on the ancestral states */
    //   sampleMutationsGivenAncestrals(mapping);

    //   // add the mapping to the vector of mapping
    //   mappings.push_back(mapping);

    //   // reset the nodes counter
    //   nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;
  }
}

/******************************************************************************/

void StochasticMapping::setExpectedAncestrals(shared_ptr<PhyloTree> expectedMapping, VVDouble& ancestralStatesFrequencies)
{
  // TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(expectedMapping);
  // vector<Node*> nodes = ttree->getNodes();
  // for (size_t i = 0; i < nodes.size(); ++i)
  // {
  //   Node* node = nodes[i];
  //   int nodeId = node->getId();
  //   size_t j = static_cast<size_t>(nodeId); //Note@Laurent (Julien 17/06/20): is this really intended, as nodeIds can be discontinuous? Should there be some can of index instead?
  //   auto d = distance(ancestralStatesFrequencies[j].begin(), max_element(ancestralStatesFrequencies[j].begin(), ancestralStatesFrequencies[j].end()));
  //   size_t state = static_cast<size_t>(d); //Note@Laurent (Julien 17/06/20): assimuming this is always positive, is that so?
  //   setNodeState(node, state); // in the case of a leaf, the assigned state must be sampled
  // }
}

/******************************************************************************/

shared_ptr<PhyloTree> StochasticMapping::generateExpectedMapping(vector<shared_ptr<PhyloTree> >& mappings, size_t divMethod)
{
  // // initialize the expected history
  // nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;
  shared_ptr<PhyloTree> expectedMapping(make_shared<PhyloTree>(*tree_));
  // setLeafsStates(expectedMapping);

  // // compute a vector of the posterior asssignment probabilities for each inner node
  // VVDouble ancestralStatesFrequencies;
  // ancestralStatesFrequencies.clear();
  // vector<Node*> nodes = dynamic_cast<TreeTemplate<Node>*>(expectedMapping)->getNodes();
  // size_t statesNum = tl_->getNumberOfStates();
  // ancestralStatesFrequencies.resize(nodes.size(), VDouble(statesNum));
  // computeStatesFrequencies(ancestralStatesFrequencies, mappings);

  // // set the ancestral states accrdonig to the maximal posterior (i.e, conditional) probability
  // setExpectedAncestrals(expectedMapping, ancestralStatesFrequencies);

  // // update the expected history with the dwelling times
  // for (size_t n = 0; n < nodes.size(); ++n)
  // {
  //   Node* node = nodes[n];
  //   if (node->hasFather()) // for any node except to the root
  //   {
  //     // initialize vector of average dwelling times for the branch stemming from node
  //     VDouble AverageDwellingTimes;
  //     AverageDwellingTimes.clear();
  //     AverageDwellingTimes.resize(statesNum, 0);
  //     // compute the average dwelling times of all the states
  //     for (size_t i = 0; i < mappings.size(); ++i)
  //     {
  //       TreeTemplate<Node>* mapping =  dynamic_cast<TreeTemplate<Node>*>(mappings[i]);
  //       // get the pointers to the node and its father in the i'th mapping
  //       Node* curNode = mapping->getNode(node->getName());
  //       Node* father = mapping->getNode(node->getFather()->getName()); // the original father of the node (according to the base tree) in the mapping
  //       while (curNode != father)
  //       {
  //         AverageDwellingTimes[static_cast<size_t>(getNodeState(curNode))] += curNode->getDistanceToFather(); //Note@Laurent (Julien 17/06/20): assuming state is positive, is that so?
  //         curNode = curNode->getFather();
  //       }
  //     }
  //     double branchLength = node->getDistanceToFather();   // this is the length of the original branch in the base tree
  //     bool updateBranch = true;
  //     for (size_t state = 0; state < statesNum; ++state)
  //     {
  //       AverageDwellingTimes[state] /= static_cast<double>(mappings.size());
  //       if (AverageDwellingTimes[state] == branchLength) // if one of the dwelling times equals the branch length, then there is only one state along te branch and there is no need to edit it
  //       {
  //         updateBranch = false;
  //       }
  //     }
  //     // break the branch according to average dwelling times
  //     if (updateBranch)
  //     {
  //       updateBranchByDwellingTimes(node, AverageDwellingTimes, ancestralStatesFrequencies, divMethod);
  //     }
  //   }
  // }
  // nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;
  return expectedMapping;
}

/******************************************************************************/

shared_ptr<PhyloTree> StochasticMapping::generateAnalyticExpectedMapping(size_t divMethod)
{
  // /* Compute the posterior assignment probabilities to internal nodes, based on the fractional probablities computed earlier */
  // const vector<int> states =  tl_->getAlphabetStates();
  // vector<int> nodeIds = baseTree_->getNodesId();
  // size_t nodeId;
  // VVDouble posteriorProbabilities;
  // posteriorProbabilities.clear();
  // posteriorProbabilities.resize(baseTree_->getNumberOfNodes(), VDouble(states.size()));
  // double nodeDataProb;
  // // because the sum of partial likelihoods (i.e, the fractional probabilities) is in fact the probablity of the data, it is sufficient to standardize the vector of fractional probabilires for each node to obtain the posterior probabilities
  // for (size_t n = 0; n < baseTree_->getNumberOfNodes(); ++n)
  // {
  //   nodeId = static_cast<size_t>(nodeIds[n]); //Note@Laurent (Julien 17/06/20): what is nodeId is negative?
  //   nodeDataProb = 0;
  //   for (size_t s = 0; s < states.size(); ++s)
  //   {
  //     nodeDataProb = nodeDataProb + fractionalProbabilities_[nodeId][s];
  //   }
  //   for (size_t nodeState = 0; nodeState < states.size(); ++nodeState)
  //   {
  //     posteriorProbabilities[nodeId][nodeState] = fractionalProbabilities_[nodeId][nodeState] / nodeDataProb;
  //   }
  // }

  // /* Assign states to internal nodes based on the majority rule over the posterior probabilities */
  shared_ptr<PhyloTree> expectedMapping(make_shared<PhyloTree>(*tree_));

  // setLeafsStates(expectedMapping);
  // setExpectedAncestrals(expectedMapping, posteriorProbabilities);

  // /* Compute the reward per state per site - expect two entries per site (that is, two entries in total).
  //    Let r0 be the reward of state 0 nd r1 the reward of state 1. */
  // UserAlphabetIndex1* alpha = new UserAlphabetIndex1(tl_->getAlphabet());
  // DiscreteDistribution* rDist = new ConstantRateDistribution();
  // TransitionModel* tlModel = tl_->getModelForSite(0, 0)->clone();
  // DRTreeLikelihood* drtl = new DRHomogeneousTreeLikelihood(*baseTree_, *(tl_->getData()), tlModel, rDist, false);
  // drtl->initialize();
  // vector<int> ids = baseTree_->getNodesId();
  // const SubstitutionModel* model = dynamic_cast<const SubstitutionModel*>(tl_->getModelForSite(0, 0));

  // /* Compute the expected dwelling times per branch and state as follows:
  //    For branch b of length t, the average welling time in state 0 is r0*t (based on Minin and Suchard paper).
  //    The average dwelling time in state 1 should complement to t (make sure of it!) */
  // vector<Node*> nodes = dynamic_cast<TreeTemplate<Node>*>(expectedMapping)->getNodes();
  // double branchLength;
  // Node* node;
  // VVDouble expectedDwellingTimes;
  // expectedDwellingTimes.clear();
  // expectedDwellingTimes.resize(nodes.size(), VDouble(states.size()));
  // for (size_t s = 0; s < states.size(); ++s)
  // {
  //   alpha->setIndex(states[s], 1); // set the reward of the state as 1 and the reward for the rest of the states as 0
  //   for (size_t m = 0; m < states.size(); ++m)
  //   {
  //     if (m != s)
  //     {
  //       alpha->setIndex(states[m], 0); //Note@Laurent (Julien 17/06/20): can you chack my correction there and above? I changed s/m to states[s] and states[m], is that correct?
  //     }
  //   }
  //   unique_ptr<Reward> reward(new DecompositionReward(model, alpha));
  //   unique_ptr<ProbabilisticRewardMapping> mapping(RewardMappingTools::computeRewardVectors(*drtl, ids, *reward, false));
  //   for (size_t n = 0; n < nodes.size(); ++n)
  //   {
  //     node = nodes[n];
  //     if (node->hasFather()) // for any node except to the root
  //     {
  //       expectedDwellingTimes[static_cast<size_t>(node->getId())][s] = mapping->getReward(node->getId(), 0); //Note@Laurent (Julien 17/06/20): what is nodeId is negative?
  //     }
  //   }
  // }

  // // standardize expected dwelling itmes, if needed, and update the mapping accorgingly
  // double sumOfDwellingTimes;
  // bool updateBranch;
  // nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;
  // for (size_t n = 0; n < nodes.size(); ++n)
  // {
  //   node = nodes[n];
  //   if (node->hasFather()) // for any node except to the root
  //   {
  //     branchLength = node->getDistanceToFather();
  //     sumOfDwellingTimes = 0;
  //     updateBranch = true;
  //     for (size_t s = 0; s < states.size(); ++s)
  //     {
  //       if (expectedDwellingTimes[static_cast<size_t>(node->getId())][s] == 0) //Note@Laurent (Julien 17/06/20): what is nodeId is negative?

  //       {
  //         updateBranch = false;
  //       }
  //       sumOfDwellingTimes = sumOfDwellingTimes + expectedDwellingTimes[static_cast<size_t>(node->getId())][s]; //Note@Laurent (Julien 17/06/20): what is nodeId is negative?

  //     }

  //     if (branchLength < 0.00001) // branch length is 0 -> no need to update mapping on the branch
  //     {
  //       node->setDistanceToFather(branchLength);
  //       updateBranch = false;
  //     }
  //     else
  //     {
  //       if (sumOfDwellingTimes != branchLength)
  //       {
  //         for (size_t s = 0; s < states.size(); ++s)
  //         {
  //           expectedDwellingTimes[static_cast<size_t>(node->getId())][s] =  branchLength * (expectedDwellingTimes[static_cast<size_t>(node->getId())][s]) / sumOfDwellingTimes;
  //         }
  //       }
  //     }

  //     if (updateBranch)
  //     {
  //       updateBranchByDwellingTimes(node, expectedDwellingTimes[static_cast<size_t>(node->getId())], posteriorProbabilities, divMethod);
  //     }
  //   }
  // }
  // nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;

  // /* free the resources */
  // delete alpha;
  // delete rDist;
  // delete tlModel;
  // delete drtl;

  return expectedMapping;
}
/******************************************************************************/

int StochasticMapping::getNodeState(const PhyloNode* node) const
{
  return (dynamic_cast<const BppInteger*>(node->getProperty(STATE)))->getValue();
}

/******************************************************************************/

void StochasticMapping::setNodeState(PhyloNode* node, size_t state)
{
  BppInteger stateProperty(static_cast<int>(state));
  node->setProperty(STATE, stateProperty);
}

/******************************************************************************/

void StochasticMapping::giveNamesToInternalNodes(PhyloTree& tree)
{
  auto nodes = tree.getAllInnerNodes();
  for (auto& node:nodes)
  {
    if (!node->hasName())
      node->setName("_baseInternal_" + TextTools::toString(tree.getNodeIndex(node)));
  }
}

/******************************************************************************/

void StochasticMapping::setLeafsStates(std::shared_ptr<PhyloTree> mapping)
{
  // auto leafsStates = likelihood_->getData();
  auto leaves = mapping->getAllLeaves();

  for (auto& leaf: leaves)
  {
    string nodeName = leaf->getName();
//    auto lstates = likelihoods_->getNode(mapping->getNodeIndex(leaf))
//   size_t leafState = static_cast<size_t>(tl_->getAlphabetStateAsInt(leafsStates->getSequence(nodeName).getValue(0)));
//     //note@Laurent (Julien on 17/06/20): I thing the above line is incorrect, in particulat the use of the getAlphabetStateAsInt function. It is supposed to take as input a state index (size_t) and return the corresponding character state as an integer. Here you give as input to the method already a sequence character (integer). In most cases that will still work as the characters states for resolved characters are usually 0..n, and there corresponding states 0..n. But it will fail for models with gaps (character state -1) and Markov modulated models (character states 0..n, but state index 0..k*n)
//     setNodeState(node, leafState);
//   }
  }
}

/******************************************************************************/


void StochasticMapping::computeFractionals()
{
  // // some auxiliiary variables
  // size_t statesNum = tl_->getNumberOfStates();
  // const TransitionModel* model = tl_->getModelForSite(0, 0); // this calls assumes that all the sites and all the branches are assoiacted with the same node
  // const SiteContainer* leafsStates = tl_->getData();
  // TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(baseTree_);
  // vector<Node*> nodes = ttree->getNodes();

  // // compute the fractional probabilities according to Felsenstein prunnig algorithm: for each node nodes[i] and state s compute: P(Data[leafs under node[i]]|node[i] has state s]
  // for (size_t i = 0; i < nodes.size(); ++i) // traverse the tree in post-order
  // {
  //   int nodeId = nodes[i]->getId();
  //   string nodeName = nodes[i]->getName();
  //   if (nodes[i]->isLeaf()) // if the node is a leaf, set the fractional probability of its state to 1, and the rest ot 0
  //   {
  //     size_t leafState = static_cast<int>(tl_->getAlphabetStateAsInt(leafsStates->getSequence(nodeName).getValue(0)));
  //     for (size_t nodeState = 0; nodeState < statesNum; ++nodeState)
  //     {
  //       if (nodeState != leafState)
  //       {
  //         fractionalProbabilities_[nodeId][nodeState] = 0;
  //       }
  //       else
  //       {
  //         fractionalProbabilities_[nodeId][nodeState] = 1;
  //       }
  //     }
  //   }
  //   else                // if the node is internal, follow the Felesenstein computation rule to compute the fractional probability
  //   {
  //     for (size_t nodeState = 0; nodeState < statesNum; ++nodeState)
  //     {
  //       double fullProb = 1;
  //       for (size_t j = 0; j < (nodes[i]->getNumberOfSons()); ++j) // for each son of the node, sum over the probabilities of all its assignments given its father's state (i.e, nodeState)
  //       {
  //         double sonProb = 0;
  //         double bl = nodes[i]->getSon(j)->getDistanceToFather();
  //         for (size_t sonState = 0; sonState < statesNum; ++sonState)
  //         {
  //           sonProb += model->Pij_t(nodeState, sonState, bl) * fractionalProbabilities_[nodes[i]->getSon(j)->getId()][sonState];
  //         }
  //         fullProb *= sonProb;
  //       }
  //       fractionalProbabilities_[nodeId][nodeState] = fullProb;
  //     }
  //   }
  // }
}

/******************************************************************************/

void StochasticMapping::ComputeConditionals()
{
  // // some auxiliiary variables
  // size_t statesNum = tl_->getNumberOfStates();
  // VDouble rootProbabilities = tl_->getRootFrequencies(0);
  // const TransitionModel* model = tl_->getModelForSite(0, 0);
  // TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(baseTree_);
  // vector<Node*> nodes = ttree->getNodes(); // here - how many nodes? only 4! what happened?

  // /* compute the fractional probabilities for each node and state */
  // fractionalProbabilities_.clear();
  // fractionalProbabilities_.resize(nodes.size(), VDouble(statesNum));
  // computeFractionals();

  // /*  compute the conditional probabilities: for each combination of nodes son, father, compute Pr(son recieves sonState | father has fatherState) */
  // ConditionalProbabilities_.clear();
  // ConditionalProbabilities_.resize((nodes.size()));

  // for (size_t i = 0; i < nodes.size(); ++i)
  // {
  //   if (!nodes[i]->isLeaf() || !nodes[i]->hasFather())  // the second condition will catch the root even if it has a single child (in which case, isLeaf() returns true)
  //   {
  //     int nodeId = nodes[i]->getId();
  //     ConditionalProbabilities_[nodes[i]->getId()].resize(statesNum, VDouble(statesNum)); // use aux variable numOfStates
  //     if (!(nodes[i]->hasFather()))  // if the node is the root -> set the conditional probability to be same for all "fatherStates"
  //     {
  //       double sum = 0.0;
  //       for (size_t sonState = 0; sonState < statesNum; ++sonState)
  //       {
  //         double stateConditionalNominator = fractionalProbabilities_[nodeId][sonState] * rootProbabilities[sonState];
  //         for (size_t fatherState = 0; fatherState < statesNum; ++fatherState)
  //         {
  //           ConditionalProbabilities_[nodeId][fatherState][sonState] = stateConditionalNominator;
  //         }
  //         sum += stateConditionalNominator;
  //       }
  //       for (size_t fatherState = 0; fatherState < statesNum; ++fatherState)
  //       {
  //         for (size_t sonState = 0; sonState < statesNum; ++sonState)
  //         {
  //           ConditionalProbabilities_[nodeId][fatherState][sonState] /= sum;
  //         }
  //       }
  //     }
  //     else                            // else -> follfow equation (10) from the paper to compute the consitional assingment probabilities given the ones of his father
  //     {
  //       for (size_t fatherState = 0; fatherState < statesNum; ++fatherState)
  //       {
  //         double sum = 0.0;
  //         for (size_t sonState = 0; sonState < statesNum; ++sonState)
  //         {
  //           double stateConditionalNominator = fractionalProbabilities_[nodes[i]->getId()][sonState] * model->Pij_t(fatherState, sonState, nodes[i]->getDistanceToFather());
  //           ConditionalProbabilities_[nodeId][fatherState][sonState] = stateConditionalNominator;
  //           sum += stateConditionalNominator;
  //         }
  //         for (size_t sonState = 0; sonState < statesNum; ++sonState)
  //         {
  //           ConditionalProbabilities_[nodeId][fatherState][sonState] /= sum;
  //         }
  //       }
  //     }
  //   }
  // }
}

/******************************************************************************/

void StochasticMapping::computeStatesFrequencies(VVDouble& ancestralStatesFrequencies, vector<shared_ptr<PhyloTree> >& mappings)
{
  // // some auxiliiary variables
  // size_t statesNum = tl_->getNumberOfStates();
  // const SiteContainer* leafsStates = tl_->getData();
  // TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(baseTree_);
  // vector<Node*> nodes = ttree->getNodes();

  // // compute the node assignment probabilities based on their frequency in the mappings
  // for (size_t i = 0; i < nodes.size(); ++i)
  // {
  //   Node* node = nodes[i];
  //   int nodeId = node->getId();
  //   string nodeName = node->getName();
  //   // in leafs - don't iterate to save time, as the frequency of a state is either 0 or 1 based on the known character data
  //   if (node->isLeaf())
  //   {
  //     size_t leafState = static_cast<int>(tl_->getAlphabetStateAsInt(leafsStates->getSequence(nodeName).getValue(0)));
  //     for (size_t nodeState = 0; nodeState < statesNum; ++nodeState)
  //     {
  //       if (nodeState != leafState)
  //       {
  //         ancestralStatesFrequencies[nodeId][nodeState] = 0;
  //       }
  //       else
  //       {
  //         ancestralStatesFrequencies[nodeId][nodeState] = 1;
  //       }
  //     }
  //   }
  //   else
  //   {
  //     // else, go over all the mappings and collect the number of states assignment per state
  //     fill(ancestralStatesFrequencies[nodeId].begin(), ancestralStatesFrequencies[nodeId].end(), 0); // reset all the values to 0
  //     for (size_t h = 0; h < mappings.size(); ++h)
  //     {
  //       Node* nodeInMapping = dynamic_cast<TreeTemplate<Node>*>(mappings[h])->getNode(nodeName);
  //       ancestralStatesFrequencies[nodeId][static_cast<size_t>(getNodeState(nodeInMapping))]++; //Note@Laurent (Julien 17/06/20): assuming node state is positive, is that so?
  //     }
  //     // now divide the vector entries by the number of mappings
  //     for (size_t nodeState = 0; nodeState < statesNum; ++nodeState)
  //     {
  //       ancestralStatesFrequencies[nodeId][nodeState] = ancestralStatesFrequencies[nodeId][nodeState] / static_cast<int>(mappings.size());
  //     }
  //   }
  // }
}


/******************************************************************************/

size_t StochasticMapping::sampleState(const VDouble& distibution)
{
  size_t state = 0;        // the default state is 0
  // double prob = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);

  // for (size_t i = 0; i < distibution.size(); ++i)
  // {
  //   prob -= distibution[i];
  //   if (prob < 0)  // if the the sampled probability is smaller than the probability to choose state i -> set state to be i
  //   {
  //     state = i;
  //     break;
  //   }
  // }
  return state;
}

/******************************************************************************/

void StochasticMapping::sampleAncestrals(shared_ptr<PhyloTree> mapping)
{
  // TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(mapping);
  // PreOrderTreeIterator* treeIt = new PreOrderTreeIterator(*ttree);
  // for (Node* node = treeIt->begin(); node != treeIt->end(); node = treeIt->next())
  // {
  //   if (!node->isLeaf())
  //   {
  //     int nodeId = node->getId();
  //     if (!node->hasFather())
  //     {
  //       size_t rootState = sampleState(ConditionalProbabilities_[nodeId][0]); // set father state to 0 (all the entries in the fatherState level are the same anyway)
  //       setNodeState(node, rootState);
  //     }
  //     else
  //     {
  //       int fatherState = getNodeState(node->getFather());
  //       size_t sonState = sampleState(ConditionalProbabilities_[nodeId][fatherState]);
  //       setNodeState(node, sonState); // in the case of a leaf, the assigned state must be sampled
  //     }
  //   }
  // }
  // delete(treeIt);
}

/******************************************************************************/

void StochasticMapping::sampleMutationsGivenAncestrals(shared_ptr<PhyloTree> mapping)
{
//   TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(mapping);
//   vector<Node*> nodes = ttree->getNodes();
//   nodesCounter_ = nodes.size() - 1; // intialize nodesCounter_ according ot the number of nodes in the base tree
//   for (size_t i = 0; i < nodes.size(); ++i)
//   {
//     Node* son = nodes[i];
//     if (son->hasFather())
//     {
//       sampleMutationsGivenAncestralsPerBranch(son);
//     }
//   }
}

/******************************************************************************/

void StochasticMapping::updateBranchMapping(PhyloNode* son, const MutationPath& branchMapping)
{
  // const vector<size_t> states = branchMapping.getStates();
  // const VDouble times = branchMapping.getTimes();
  // Node* curNode = son;
  // Node* nextNode;
  // int eventsNum = static_cast<int>(branchMapping.getNumberOfEvents());

  // if (eventsNum == 0) // if there are no events >-return nothing
  //   return;
  // else
  // {
  //   for (int i = eventsNum - 1; i > -1; --i) // add a new node to represent the transition
  //   {
  //     nodesCounter_ = nodesCounter_ + 1;
  //     const string name = "_mappingInternal" + TextTools::toString(nodesCounter_) + "_";
  //     nextNode = new Node(static_cast<int>(nodesCounter_), name);
  //     setNodeState(nextNode, states[i]);
  //     nextNode->setDistanceToFather(times[i]);

  //     // set the father to no longer be the father of curNode
  //     Node* originalFather = curNode->getFather();
  //     originalFather->removeSon(curNode); // also removes originalFather as the father of curNode
  //     // set nextNode to be the new father of curNode
  //     curNode->setFather(nextNode); // also adds curNode to the sons of nextNode
  //     // set curNode's original father ot be the father of nextNode
  //     nextNode->setFather(originalFather); // also adds nextNode to the sons of originalFather - or not? make sure this is the father at all times
  //     curNode = nextNode;
  //   }
  //   return;
  // }
}

/******************************************************************************/

void StochasticMapping::sampleMutationsGivenAncestralsPerBranch(PhyloNode* son, size_t maxIterNum)
{
  // Node* father = son->getFather();
  // double branchLength = son->getDistanceToFather();
  // size_t fatherState = getNodeState(father);
  // size_t sonState = getNodeState(son);

  // /* simulate mapping on a branch until you manage to finish at the son's state */
  // for (size_t i = 0; i < maxIterNum; ++i)
  // {
  //   double disFromNode = 0.0;
  //   size_t curState = fatherState;
  //   MutationPath tryMapping(mappingParameters_->getSubstitutionModel()->getAlphabet(), fatherState, branchLength);

  //   double timeTillChange;
  //   // if the father's state is not the same as the son's state -> use the correction corresponding to equation (11) in the paper
  //   if (fatherState != sonState)
  //   {   // sample timeTillChange conditional on it being smaller than branchLength
  //     double u = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
  //     double waitingTimeParam = -1 * mappingParameters_->getSubstitutionModel()->Qij(fatherState, fatherState); // get the parameter for the exoponential distribution to draw the waiting time from
  //     double tmp = u * (1.0 - exp(branchLength * -waitingTimeParam));
  //     timeTillChange =  -log(1.0 - tmp) / waitingTimeParam;
  //     assert (timeTillChange < branchLength);
  //   }
  //   else
  //   {
  //     timeTillChange = mappingParameters_->getTimeBeforeNextMutationEvent(fatherState); // draw the time until a transition from exponential distribution with the rate of leaving fatherState
  //   }

  //   while (disFromNode + timeTillChange < branchLength)  // a jump occured but not passed the whole branch ->
  //   {
  //     tryMapping.addEvent(curState, timeTillChange);                                        // add the current state and time to branch history

  //     disFromNode += timeTillChange;
  //     timeTillChange = mappingParameters_->getTimeBeforeNextMutationEvent(curState);        // draw the time until a transition from exponential distribution with the rate of leaving curState
  //     curState = mappingParameters_->mutate(curState);                                      // draw the state to transition to after from initial state curState based on the relative tranistion rates distribution (see MutationProcess.cpp line 50)
  //   }
  //   // the last jump passed the length of the branch -> finish the simulation and check if it's sucessfull (i.e, mapping is finished at the son's state)
  //   if (curState != sonState) // if the simulation failed, try again
  //   {
  //     continue;
  //   }
  //   else                      // if the simulation was sucessfully, add it to the build mapping
  //   {
  //     double timeOfJump = branchLength - disFromNode;
  //     son->setDistanceToFather(timeOfJump);
  //     updateBranchMapping(son, tryMapping);     // add the successfull simulation to the build mapping
  //     return;
  //   }
  // }
  // // if all simulations failed -> throw an exception
  // throw Exception("could not produce simulations with father = " + TextTools::toString(fatherState) + " son " + TextTools::toString(sonState) + " branch length = " + TextTools::toString(branchLength));
}

/******************************************************************************/

void StochasticMapping::updateBranchByDwellingTimes(PhyloNode* node, VDouble& dwellingTimes, VVDouble& ancestralStatesFrequencies, size_t divMethod)
{
  // /* first, convert the dwelling times vector to a mutation path of the branch */
  // size_t statesNum = tl_->getNumberOfStates();
  // size_t sonState = getNodeState(node);
  // size_t fatherState = getNodeState(node->getFather());
  // double branchLength = node->getDistanceToFather();
  // double Pf = 1;
  // double Ps = 1;
  // double shareOfFather = 0;
  // double shareOfSon = 0;
  // MutationPath branchMapping(mappingParameters_->getSubstitutionModel()->getAlphabet(), fatherState, branchLength);
  // // set the first event with the dwelling time that matches the state of the father
  // if (fatherState == sonState)
  // {
  //   if (divMethod == 0)
  //   {
  //     if (node->hasFather())
  //     {
  //       Pf = ancestralStatesFrequencies[node->getFather()->getId()][fatherState];
  //     }
  //     Ps = 1;
  //     if (!node->isLeaf())
  //     {
  //       Ps = ancestralStatesFrequencies[node->getId()][sonState];
  //     }
  //     shareOfFather = Pf / (Pf + Ps);
  //     branchMapping.addEvent(fatherState, dwellingTimes[fatherState] * shareOfFather);
  //   }
  //   else
  //   {
  //     branchMapping.addEvent(fatherState, 0);
  //   }
  // }
  // else
  // {
  //   branchMapping.addEvent(fatherState, dwellingTimes[fatherState]);
  // }
  // // set all events except for the one entering the son
  // for (size_t state = 0; state < statesNum; ++state)
  // {
  //   if (state != fatherState && state != sonState && dwellingTimes[state] > 0) // if the state matches an event which is not the first or the last -> add it
  //   {
  //     branchMapping.addEvent(state, dwellingTimes[state]);
  //   }
  // }
  // // change the length of the branch whose bottom node is the son according to the dwelling time of the relevant state
  // if (fatherState == sonState)
  // {
  //   if (divMethod == 0)
  //   {
  //     shareOfSon = 1 - shareOfFather;
  //     node->setDistanceToFather(dwellingTimes[sonState] * shareOfSon);
  //   }
  //   else
  //   {
  //     node->setDistanceToFather(dwellingTimes[sonState]);
  //   }
  // }
  // else
  // {
  //   node->setDistanceToFather(dwellingTimes[sonState]);
  // }

  // /* secondly, update the expected history with the dwelling times-based mutation path */
  // updateBranchMapping(node, branchMapping);
}
