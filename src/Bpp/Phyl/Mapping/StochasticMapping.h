//
// File: StochasticMapping.h
// Created by: Keren Halabi
// Created on: June 2018
//

/*
   Copyright or © or Copr. CNRS, (November 16, 2004)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef ___STOCHASTIC_MAPPING_H
#define ___STOCHASTIC_MAPPING_H

#include "../Likelihood/TreeLikelihood.h"
#include "../Simulation/MutationProcess.h"

// From the STL:
#include <iostream>
#include <iomanip>
#include <map>

using namespace std;

typedef vector<vector<vector<double> > > VVVDouble;
typedef vector<vector<double> > VVDouble;
typedef vector<double> VDouble;

/* class for reprenting the framework of Stochastic mapping
   A StochasticMapping instance can be used to sample histories of state transitions along a tree, given a substitution model and the states at the tip taxa
   For more information, see: Nielsen, Rasmus. "Mapping mutations on phylogenies." Systematic biology 51.5 (2002): 729-739.‏ */

namespace bpp
{
class StochasticMapping
{
protected:
  const SimpleMutationProcess* mappingParameters_; // this instance will hold the parameters required for the sotchastic mapping procedure, and be used to generate stochastic mappings
  Tree* baseTree_;                                 // this is the base tree, which will act as the skeleton of each induced mapping in the procedure
  const TreeLikelihood* tl_;                       // the tree likelihood instance is used for computing the the conditional sampling probabilities of the ancestral states as well as the root assignment probabilities
  VVDouble fractionalProbabilities_;               // vector that holds the fractional probabilities per state per node in the tree, based on which the conditional and posterior probabilities are computed
  VVVDouble ConditionalProbabilities_;             // vector that holds the conditionl states assignment probabilities of the nodes in the tree (node*father_states*son_states)
  size_t nodesCounter_;                            // counter of nodes hat allows adding unique names to the generated nodes while breaking branching in a mapping
  size_t numOfMappings_;                           // the number of stochastic mappings to generate
  map<int,size_t> nodeIdToIndex_;

public:
  /* constructors and destructors */

  explicit StochasticMapping(const TreeLikelihood* tl, size_t numOfMappings = 10000); // it is a good general practice to use "explicit" keyword on constructors with a single argument: https://stackoverflow.com/questions/121162/what-does-the-explicit-keyword-mean

  ~StochasticMapping();

  StochasticMapping(const StochasticMapping& sm) : // must pass sm by repference to avoid infinitie recusion in the copy construcor
    mappingParameters_(sm.mappingParameters_), baseTree_(0), tl_(sm.tl_), fractionalProbabilities_(sm.fractionalProbabilities_), ConditionalProbabilities_(sm.ConditionalProbabilities_), nodesCounter_(0), numOfMappings_(sm.numOfMappings_), nodeIdToIndex_(sm.nodeIdToIndex_)
  { baseTree_ = sm.baseTree_->clone(); } // the tree must be cloned so that instead of copying the pointer to the tree, a new tree with a new pointer will be created

  /**
   * @brief Assignment operator
   */
  StochasticMapping& operator=(const StochasticMapping& sm)
  {
    tl_ = sm.tl_;
    baseTree_ = sm.baseTree_->clone();
    mappingParameters_ = sm.mappingParameters_;
    fractionalProbabilities_ = sm.fractionalProbabilities_;
    ConditionalProbabilities_ = sm.ConditionalProbabilities_;
    numOfMappings_ = sm.numOfMappings_;
    nodeIdToIndex_ = sm.nodeIdToIndex_;
    return *this;
  }

  /**
   * @cloning function used by the copy constructor of ./Likelihood/JointLikelihoodFunction/h
   */
  StochasticMapping* clone() const { return new StochasticMapping(*this); }


  /* generates a stochastic mappings based on the sampling parameters
   * @param     Number of histories to sample
   * @param mappings          Vector of Tree pointers that will hold the sampled stochastic mappings.
   *                          Note that this function generates new tree instances, that must be deleted by the calling function.
   */
  void generateStochasticMapping(vector<Tree*>& mappings);

  /* creates a single expected (i.e, average) history based on a given set of mappings
   * steps correspond to Nielsen, Rasmus. "Mapping mutations on phylogenies." Systematic biology 51.5 (2002): 729-739.‏
   * the function assumes that there is only one site to simulate history for */
  /* @param mappings          A vector of stochastic mappings to average
     @param divMethod         The method used in the case that the son and father share the same state (either divide the wdelling time of the staed state by 2 for  two transitions (method 0) or allocate the entire dwelling time to be adjacent to the son(method 1))
   */
  Tree* generateExpectedMapping(const vector<Tree*>& mappings, size_t divMethod = 0);

  /* creates a single expected (i.e, average) history based the rewards prvided by te algorithm of Minin and Suchard (2008)
   * the function assumes that there is only one site to simulate history for */
  /* @param divMethod         The method used in the case that the son and father share the same state (either divide the wdelling time of the staed state by 2 for  two transitions (method 0) or allocate the entire dwelling time to be adjacent to the son(method 1))
   */
  Tree* generateAnalyticExpectedMapping(size_t divMethod = 0);

  /* extracts the state of a node in a mapping
   * @param node              The node to get the state of
   * @return                  Node state is int
   */
  static size_t getNodeState(const Node* node);

  /* sets the state of a node in a mapping
   * @param node               The node to get the state of
   * @param state              The state that needs to be assigned to the node
   */
  static void setNodeState(Node* node, size_t state);

  /* compute the ancestral frequenceis of character states of all the nodes based on the mappings
  * @param                     A vector of the posterior probabilities probabilities to fill in (node**state combinaion in each entry)
  * @param                     A vector of mappings to base the frequencies on
  */
  void computeStatesFrequencies(VVDouble& ancestralStatesFreuquencies, const vector<Tree*>& mappings);


private:
  /* adds names to the internal nodes, in case of absence.
   * @param tree               The tree whose nodes should be edited if needed.
   */
  void giveNamesToInternalNodes(Tree* tree);

  /* returns the possible model states assigned to a leaf based on its character state
   * @param node - the leaf of interest
   */
  vector<size_t> getLeafModelStates(Node* node);

  /* returns the possible model states assignments of the leafs as properties of thier nodes instances
   * @param mapping - the tree to sets the properties in
   */
  map<int,vector<size_t>> setLeafsStates(Tree* mapping);

  /* compute the fractional probabilities of all the nodes assignements
   * @param                     A vector of the fractional probabilities probabilities to fill in (node**state combinaion in each entry)
   */
  void computeFractionals();

  /* compute the conditional probabilities of all the nodes assignments
   * @param rootProbabilities  The root frequencies
   */
  void ComputeConditionals();

  /* auxiliary function that samples a state based on a given discrete distribution
   * @param distibution       The distribution to sample states based on
   */
  size_t sampleState(const VDouble& distibution); // k: best by ref

  /* samples ancestral states based on the conditional probabilities at each node in the base (user input) tree and the root assignment probabilities. States will be updated as nodes properties
   * @param mapping               The tree whose nodes names should be updated according to their assigned states.
   */
  void sampleAncestrals(Tree* mapping, map<int,vector<size_t>> leafIdToStates);

  /* set ancestral states in the expected history based on the conditional probabilities at each node in the base (user input) tree and the root assignment probabilities. States will be updated as nodes properties
   * @param expectedMapping           The expected mapping instance whose nodes names should be updated according to their assigned states.
   * @param posteriorProbabilities    Vector of posterior assignment proabilities to inner node to decide on assignments
   */
  void setExpectedAncestrals(Tree* expectedMapping, VVDouble& posteriorProbabilities);

  /* simulates mutations on phylogeny based the sampled ancestrals, tips data, and the simulation parameters
   * @param mapping               The tree that should be edited according to the sampled history
   */
  void sampleMutationsGivenAncestrals(Tree* mapping);

  /* adds a branch mapping to the mapping in a tree format by repeatedly braking branches and adding internal nodes with single children
   * @param son                   The node at the bottom of the branch
   * @param branchMapping         The branchMapping of transitions in a MutationProcess format
   */
  void updateBranchMapping(Node* son, const MutationPath& branchMapping);

  /* sample mutations based on the stochastic mapping parameters, the source and destination state, and the branch length, and updates the simulated history along the branch in the input tree, on the fly
   * @param son                   Node of interest
   * @param maxIterNum            Maximal number of imulation trials
   */
  void sampleMutationsGivenAncestralsPerBranch(Node* son, size_t maxIterNum = 10000);

  /* converts a vector of dwelling times to a mutation path and then updates the bracnh stemming from the given node */
  /* @param node                      The node at the bottom of the branch
   * @param dwellingTimes             A vector of dwelling times where the value at each entry i corresponds to the dwelling time under the i'th state
   *                                  Note that this function generates a new tree instance, that must be deleted by the calling function.
   * @param posteriorProbabilities    Posterior probaibitlies to divide the time spent in a shared state between father and son into two transitions
     @param divMethod                 The method used in the case that the son and father share the same state (either divide the wdelling time of the staed state by 2 for  two transitions (method 0) or allocate the entire dwelling time to be adjacent to the son(method 1))
   */
  void updateBranchByDwellingTimes(Node* node, VDouble& dwellingTimes, VVDouble& posteriorProbabilities, size_t divMethod = 0);
};
}

#endif// ___STOCHASTIC_MAPPING_H