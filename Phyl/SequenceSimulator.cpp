//
// File: SequenceSimulator.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Feb  4 16:30:51 2004
//

#include "SequenceSimulator.h"

// From SeqLib:
#include <Seq/VectorSiteContainer.h>

/******************************************************************************/

HomogeneousSequenceSimulator::HomogeneousSequenceSimulator(
	const MutationProcess * process,
	const DiscreteDistribution * rate,
	const Tree * tree):
	_process(process),
	_rate(rate),
	_tree(tree)
{
	_leaves   = _tree    -> getLeaves();
	_model    = _process -> getSubstitutionModel();
	_alphabet = _model   -> getAlphabet();
	_seqNames.resize(_leaves.size());
	for(unsigned int i = 0; i < _seqNames.size(); i++) _seqNames[i] = _leaves[i] -> getName();
	_hssr = NULL;
}

/******************************************************************************/

void HomogeneousSequenceSimulator::evolveInternal(const Node * node, double rate) const
{
	if(!node -> hasFather()) { cerr << "DEBUG: HomogeneousSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl; return; }
	int initialState = _states[node -> getFather()];
	int   finalState = _process -> evolve(initialState, node -> getDistanceToFather() * rate);
	_states[node] = finalState;
	for(unsigned int i = 0; i < node -> getNumberOfSons(); i++) {
		evolveInternal(node -> getSon(i), rate);
	}
}

/******************************************************************************/

void HomogeneousSequenceSimulator::preciseEvolveInternal(const Node * node, double rate) const
{
	if(!node -> hasFather()) { cerr << "DEBUG: HomogeneousSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl; return; }
	int initialState = _states[node -> getFather()];
	MutationPath mp = _process -> detailedEvolve(initialState, node -> getDistanceToFather() * rate);
	_states[node] = mp.getFinalState();

	// Now append infos in _ssr:
	_hssr -> addNode(node, mp);

	// Now jump to son nodes:
	for(unsigned int i = 0; i < node -> getNumberOfSons(); i++) {
		preciseEvolveInternal(node -> getSon(i), rate);
	}
}
		
/******************************************************************************/
	
Site * HomogeneousSequenceSimulator::evolve(int initialState, double rate) const
{
	// Launch recursion:
	const Node * root = _tree -> getRootNode();
	_states[root] = initialState;
	for(unsigned int i = 0; i < root -> getNumberOfSons(); i++) {
		evolveInternal(root -> getSon(i), rate);
	}
	// Now create a Site object:
	Vint site(_leaves.size());
	for(unsigned int i = 0; i < _leaves.size(); i++) {
		site[i] = _states[_leaves[i]];
	}
	return new Site(site, _alphabet);
}

/******************************************************************************/
	
void HomogeneousSequenceSimulator::preciseEvolve(int initialState, double rate) const
{
	// Launch recursion:
	const Node * root = _tree -> getRootNode();
	_states[root] = initialState;
	for(unsigned int i = 0; i < root -> getNumberOfSons(); i++) {
		preciseEvolveInternal(root -> getSon(i), rate);
	}
}

/******************************************************************************/

Site * HomogeneousSequenceSimulator::simulate() const
{
	// Draw an initial state randomly according to equilibrum frequencies:
	int initialState = 0;
	double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
	double cumprob = 0;
	for(unsigned int i = 0; i < _alphabet -> getSize(); i++)	{
		cumprob += _model -> freq(i);
		if(r <= cumprob) {
			initialState = i;
			break;
		}
	}
	// Draw a random rate:
	double rate = _rate -> rand();
	// Make this state evolve:
	return evolve(initialState, rate);
}

/******************************************************************************/

SiteContainer * HomogeneousSequenceSimulator::simulate(unsigned int numberOfSites) const
{
	vector<const Site *> vs(numberOfSites);
	for(unsigned int i = 0; i < numberOfSites; i++) {
		Site * s = simulate();
		s -> setPosition((int)i);
		vs[i] = s;
	}
	SiteContainer * sites = new VectorSiteContainer(vs, _alphabet);
	sites -> setSequencesNames(_seqNames, false);
	// Freeing memory:
	for(unsigned int i = 0; i < numberOfSites; i++) delete vs[i];
	return sites;
}

/******************************************************************************/

SequenceSimulationResult * HomogeneousSequenceSimulator::preciseSimulate() const
{
	// Draw an initial state randomly according to equilibrum frequencies:
	int initialState = 0;
	double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
	double cumprob = 0;
	for(unsigned int i = 0; i < _alphabet -> getSize(); i++)	{
		cumprob += _model -> freq(i);
		if(r <= cumprob) {
			initialState = i;
			break;
		}
	}
	// Draw a random rate:
	double rate = _rate -> rand();
	// Make this state evolve:
	_hssr = new HomogeneousSequenceSimulationResult(_tree, initialState, rate);
	preciseEvolve(initialState, rate);
	return _hssr;
}

/******************************************************************************/

const MutationProcess * 
HomogeneousSequenceSimulator::getMutationProcess() const { return _process; }

const SubstitutionModel * 
HomogeneousSequenceSimulator::getSubstitutionModel() const {
	return _process -> getSubstitutionModel();
}

const DiscreteDistribution * 
HomogeneousSequenceSimulator::getRateDistribution() const { return _rate; }

const Tree * 
HomogeneousSequenceSimulator::getTree() const { return _tree; }

/******************************************************************************/
