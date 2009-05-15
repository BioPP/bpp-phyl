//
// File: NonHomogeneousSequenceSimulator.cpp
//       (previously SequenceSimulator.cpp, then HomogeneousSequenceSimulator.cpp)
// Created by: Julien Dutheil
//             Bastien Boussau
// Created on: Wed Feb  4 16:30:51 2004
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

#include "NonHomogeneousSequenceSimulator.h"
#include "SubstitutionModelSetTools.h"

// From Utils:
#include <Utils/ApplicationTools.h>

// From SeqLib:
#include <Seq/VectorSiteContainer.h>

// From NumCalc:
#include <NumCalc/VectorTools.h>
#include <NumCalc/MatrixTools.h>

using namespace bpp;

/******************************************************************************/

NonHomogeneousSequenceSimulator::NonHomogeneousSequenceSimulator(
  const SubstitutionModelSet * modelSet,
  const DiscreteDistribution * rate,
  const TreeTemplate<Node> * tree) throw(Exception):
  _modelSet(modelSet),
  _rate(rate),
  _templateTree(tree),
  _tree(*tree),
  _ownModelSet(false),
  _continuousRates(false)
{
  if(!modelSet->isFullySetUpFor(*tree))
    throw Exception("NonHomogeneousSequenceSimulator(constructor). Model set is not fully specified.");
  init();
}

/******************************************************************************/

NonHomogeneousSequenceSimulator::NonHomogeneousSequenceSimulator(
  const SubstitutionModel * model,
  const DiscreteDistribution * rate,
  const TreeTemplate<Node> * tree):
  _modelSet(NULL),
  _rate(rate),
  _templateTree(tree),
  _tree(*tree),
  _ownModelSet(true),
  _continuousRates(false)
{
  FullFrequenciesSet* fSet = new FullFrequenciesSet(model->getAlphabet(), "anc");
  fSet->setFrequencies(model->getFrequencies());
  _modelSet = SubstitutionModelSetTools::createHomogeneousModelSet(dynamic_cast<SubstitutionModel*>(model->clone()), fSet, _templateTree);
  init();
}
  
/******************************************************************************/

void NonHomogeneousSequenceSimulator::init()
{
  _leaves   = _tree.getLeaves();
  _alphabet = _modelSet->getAlphabet();
  _seqNames.resize(_leaves.size());
  for(unsigned int i = 0; i < _seqNames.size(); i++)
    _seqNames[i] = _leaves[i]->getName();
  // Initialize cumulative pxy:
  vector<SNode *> nodes = _tree.getNodes();
  nodes.pop_back(); //remove root
  _nbNodes = nodes.size();
  _nbClasses = _rate->getNumberOfCategories();
  _nbStates = _modelSet->getNumberOfStates();
  
  for(unsigned int i = 0; i < nodes.size(); i++)
  {
    SNode * node = nodes[i];
    node->getInfos().model = _modelSet->getModelForNode(node->getId());
    double d = node->getDistanceToFather();
    VVVdouble * _cumpxy_node = &node->getInfos().cumpxy;
    _cumpxy_node->resize(_nbClasses);
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      VVdouble * _cumpxy_node_c = & (* _cumpxy_node)[c];
      _cumpxy_node_c->resize(_nbStates);
      RowMatrix<double> P = node->getInfos().model->getPij_t(d * _rate->getCategory(c));
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        Vdouble * _cumpxy_node_c_x = & (* _cumpxy_node_c)[x];
        _cumpxy_node_c_x->resize(_nbStates);
        (* _cumpxy_node_c_x)[0] = P(x, 0);
        for(unsigned int y = 1; y < _nbStates; y++)
        {
          (* _cumpxy_node_c_x)[y] = (* _cumpxy_node_c_x)[y - 1] + P(x, y);
        }
      }
    }
  }
}

/******************************************************************************/

Site * NonHomogeneousSequenceSimulator::simulate() const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  int initialState = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  vector<double> freqs = _modelSet->getRootFrequencies();
  for(unsigned int i = 0; i < _nbStates; i++)
  {
    cumprob += freqs[i];
    if(r <= cumprob)
    {
      initialState = (int)i;
      break;
    }
  }
  return simulate(initialState);
}

/******************************************************************************/
  
Site * NonHomogeneousSequenceSimulator::simulate(int initialState) const
{
  if(_continuousRates)
  {
    // Draw a random rate:
    double rate = _rate->randC();
    // Make this state evolve:
    return simulate(initialState, rate);
  }
  else
  {
    // Draw a random rate:
    unsigned int rateClass = (unsigned int)RandomTools::giveIntRandomNumberBetweenZeroAndEntry(_rate->getNumberOfCategories());
    // Make this state evolve:
    return simulate(initialState, rateClass);
  }
}

/******************************************************************************/

Site * NonHomogeneousSequenceSimulator::simulate(int initialState, unsigned int rateClass) const
{
  // Launch recursion:
  SNode* root = _tree.getRootNode();
  root->getInfos().state = initialState;
  for(unsigned int i = 0; i < root->getNumberOfSons(); i++)
  {
    evolveInternal(root->getSon(i), rateClass);
  }
  // Now create a Site object:
  Vint site(_leaves.size());
  for(unsigned int i = 0; i < _leaves.size(); i++)
  {
    site[i] = _leaves[i]->getInfos().model->getState(_leaves[i]->getInfos().state);
  }
  return new Site(site, _alphabet);
}

/******************************************************************************/
  
Site * NonHomogeneousSequenceSimulator::simulate(int initialState, double rate) const
{
  // Launch recursion:
  SNode* root = _tree.getRootNode();
  root->getInfos().state = initialState;
  for(unsigned int i = 0; i < root->getNumberOfSons(); i++)
  {
    evolveInternal(root->getSon(i), rate);
  }
  // Now create a Site object:
  Vint site(_leaves.size());
  for(unsigned int i = 0; i < _leaves.size(); i++)
  {
    site[i] = _leaves[i]->getInfos().model->getState(_leaves[i]->getInfos().state);
  }
  return new Site(site, _alphabet);
}

/******************************************************************************/

Site * NonHomogeneousSequenceSimulator::simulate(double rate) const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  int initialState = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;  
  vector<double> freqs = _modelSet->getRootFrequencies();
  for(unsigned int i = 0; i < _nbStates; i++)
  {
    cumprob += freqs[i];
    if(r <= cumprob)
    {
      initialState = (int)i;
      break;
    }
  }
  // Make this state evolve:
  return simulate(initialState, rate);
}

/******************************************************************************/
    
SiteContainer * NonHomogeneousSequenceSimulator::simulate(unsigned int numberOfSites) const
{
  Vint * initialStates = new Vint(numberOfSites, 0);
  for(unsigned int j = 0; j < numberOfSites; j++)
  { 
    double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
    double cumprob = 0;  
    vector<double> freqs = _modelSet->getRootFrequencies();
    for(unsigned int i = 0; i < _nbStates; i++)
    {
      cumprob += freqs[i]; 
      if(r <= cumprob)
      {
        (* initialStates)[j] = (int)i;
        break;
      }
    }
  }
  if(_continuousRates)
  {
    VectorSiteContainer * sites = new VectorSiteContainer(_seqNames.size(), _alphabet);
    sites->setSequencesNames(_seqNames);
    for(unsigned int j = 0; j < numberOfSites; j++)
    {
      Site * site = simulate();
      site->setPosition(j);
      sites->addSite(*site);
      delete site;
    }
    return sites;
  }
  else
  { 
    //More efficient to do site this way:
    // Draw random rates:
    vector<unsigned int> rateClasses(numberOfSites);
    unsigned int nCat = _rate->getNumberOfCategories();
    for(unsigned int j = 0; j < numberOfSites; j++)
      rateClasses[j] = (unsigned int)RandomTools::giveIntRandomNumberBetweenZeroAndEntry(nCat);
    // Make these states evolve:
    SiteContainer * sites = multipleEvolve(* initialStates, rateClasses);
    return sites;
  }
}

/******************************************************************************/

RASiteSimulationResult * NonHomogeneousSequenceSimulator::dSimulate() const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  int initialState = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;  
  vector<double> freqs = _modelSet->getRootFrequencies();
  for(unsigned int i = 0; i < _nbStates; i++)
  {
    cumprob += freqs[i];
    if(r <= cumprob)
    {
      initialState = i;
      break;
    }
  }
  
  return dSimulate(initialState);
}

/******************************************************************************/

RASiteSimulationResult * NonHomogeneousSequenceSimulator::dSimulate(int initialState) const
{
  // Draw a random rate:
  if(_continuousRates)
  {
    double rate = _rate->randC();
    return dSimulate(initialState, rate);
  }
  else
  {
    unsigned int rateClass = (unsigned int)RandomTools::giveIntRandomNumberBetweenZeroAndEntry(_rate->getNumberOfCategories());
    return dSimulate(initialState, rateClass);
    //NB: this is more efficient than dSimulate(initialState, _rDist->rand())
  }
}

/******************************************************************************/

RASiteSimulationResult * NonHomogeneousSequenceSimulator::dSimulate(int initialState, double rate) const
{
  // Make this state evolve:
  RASiteSimulationResult * hssr = new RASiteSimulationResult(_templateTree, _modelSet->getAlphabet(), initialState, rate);
  dEvolve(initialState, rate, *hssr);
  return hssr;
}

/******************************************************************************/

RASiteSimulationResult * NonHomogeneousSequenceSimulator::dSimulate(int initialState, unsigned int rateClass) const
{
  return dSimulate(initialState, _rate->getCategory(rateClass));
}

/******************************************************************************/

RASiteSimulationResult * NonHomogeneousSequenceSimulator::dSimulate(double rate) const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  int initialState = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  vector <double> freqs =  _modelSet->getRootFrequencies();
  for(unsigned int i = 0; i < _nbStates; i++)
  {
    cumprob += freqs[i];
    if(r <= cumprob)
    {
      initialState = i;
      break;
    }
  }
  return dSimulate(initialState, rate);
}

/******************************************************************************/

int NonHomogeneousSequenceSimulator::evolve(const SNode* node, int initialState, unsigned int rateClass) const
{
  const Vdouble * _cumpxy_node_c_x = & node->getInfos().cumpxy[rateClass][initialState];
  double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  for(int y = 0; y < (int)_nbStates; y++)
  {
    if(rand < (* _cumpxy_node_c_x)[y]) return y;
  }
  cerr << "DEBUG: This message should never happen! (HomogeneousSequenceSimulator::evolve)" << endl;
  cout << "   rand = " << rand << endl;
  return -1;
}

/******************************************************************************/

int NonHomogeneousSequenceSimulator::evolve(const SNode* node, int initialState, double rate) const
{
  double cumpxy = 0;
  double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double l = rate * node->getDistanceToFather();
  const SubstitutionModel* model = node->getInfos().model;
  for(int y = 0; y < (int)_nbStates; y++)
  {
    cumpxy += model->Pij_t(initialState, y, l);
    if(rand < cumpxy) return y;
  }
  cerr << "DEBUG: This message should never happen! (NonHomogeneousSequenceSimulator::evolve)" << endl;
  cout << "  rand = " << rand << endl;
  cout << "cumpxy = " << cumpxy << endl;
  MatrixTools::print(model->getPij_t(l));
  return -1;
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::multipleEvolve(const SNode* node, const Vint & initialStates, const vector<unsigned int> & rateClasses, Vint & finalStates) const
{
  const VVVdouble * _cumpxy_node = & node->getInfos().cumpxy;
  for(unsigned int i = 0; i < initialStates.size(); i++)
  {
    const Vdouble * _cumpxy_node_c_x = & (* _cumpxy_node)[rateClasses[i]][initialStates[i]];
    double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
    for(unsigned int y = 0; y < _nbStates; y++)
    {
      if(rand < (* _cumpxy_node_c_x)[y])
      {
        finalStates[i] = (int)y;
        break;
      }
    }
  }
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::evolveInternal(SNode* node, unsigned int rateClass) const
{
  if(!node->hasFather())
  { 
    cerr << "DEBUG: NonHomogeneousSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  node->getInfos().state = evolve(node, node->getFather()->getInfos().state, rateClass);
  for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    evolveInternal(node->getSon(i), rateClass);
  }
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::evolveInternal(SNode* node, double rate) const
{
  if(!node->hasFather())
  { 
    cerr << "DEBUG: NonHomogeneousSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  node->getInfos().state = evolve(node, node->getFather()->getInfos().state, rate);
  for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    evolveInternal(node->getSon(i), rate);
  }
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::multipleEvolveInternal(SNode * node, const vector<unsigned int> & rateClasses) const
{
  if(!node->hasFather())
  { 
    cerr << "DEBUG: NonHomogeneousSequenceSimulator::multipleEvolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  const vector<int> * initialStates = &node->getFather()->getInfos().states;
  unsigned int n = initialStates->size();
  node->getInfos().states.resize(n); //allocation.
  multipleEvolve(node, node->getFather()->getInfos().states, rateClasses, node->getInfos().states);
  for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    multipleEvolveInternal(node->getSon(i), rateClasses);
  }
}

/******************************************************************************/

SiteContainer * NonHomogeneousSequenceSimulator::multipleEvolve(const Vint & initialStates, const vector<unsigned int> & rateClasses) const
{
  // Launch recursion:
  SNode * root = _tree.getRootNode();
  root->getInfos().states = initialStates;
  for(unsigned int i = 0; i < root->getNumberOfSons(); i++)
  {
    multipleEvolveInternal(root->getSon(i), rateClasses);
  }
  // Now create a SiteContainer object:
  AlignedSequenceContainer * sites = new AlignedSequenceContainer(_alphabet);
  unsigned int n = _leaves.size();
  unsigned int nbSites = initialStates.size();
  const SubstitutionModel * model = NULL;
  for(unsigned int i = 0; i < n; i++)
  {
    vector<int> content(nbSites);
    vector<int> * states = &_leaves[i]->getInfos().states;
    model = _leaves[i]->getInfos().model;
    for(unsigned int j = 0; j < nbSites; j++)
    {
      content[j] = model->getState((*states)[j]);
    }
    sites->addSequence(Sequence(_leaves[i]->getName(), content, _alphabet), false);
  }
  return sites;
}

/******************************************************************************/
  
void NonHomogeneousSequenceSimulator::dEvolve(int initialState, double rate, RASiteSimulationResult & rassr) const
{
  // Launch recursion:
  SNode * root = _tree.getRootNode();
  root->getInfos().state = initialState;
  for(unsigned int i = 0; i < root->getNumberOfSons(); i++)
  {
    dEvolveInternal(root->getSon(i), rate, rassr);
  }
}

/******************************************************************************/

void NonHomogeneousSequenceSimulator::dEvolveInternal(SNode * node, double rate, RASiteSimulationResult & rassr) const
{
  if(!node->hasFather())
  { 
    cerr << "DEBUG: NonHomogeneousSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  SimpleMutationProcess* process = new SimpleMutationProcess(node->getInfos().model);
  MutationPath mp = process->detailedEvolve(node->getFather()->getInfos().state, node->getDistanceToFather() * rate);
  node->getInfos().state = mp.getFinalState();

  // Now append infos in rassr:
  rassr.addNode(node->getId(), mp);

  // Now jump to son nodes:
  for(unsigned int i = 0; i < node->getNumberOfSons(); i++)
  {
    dEvolveInternal(node->getSon(i), rate, rassr);
  }
}

/******************************************************************************/

