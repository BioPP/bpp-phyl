//
// File: SequenceSimulator.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Feb  4 16:30:51 2004
//

#include "SequenceSimulator.h"
#include "ApplicationTools.h"

// From SeqLib:
#include <Seq/VectorSiteContainer.h>

// From NumCalc:
#include <NumCalc/VectorTools.h>
using namespace VectorOperators;

/******************************************************************************/

HomogeneousSequenceSimulator::HomogeneousSequenceSimulator(
	const MutationProcess * process,
	const DiscreteDistribution * rate,
	const Tree<Node> * tree,
	bool verbose):
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
	// Initialize cumulative pxy:
	if(verbose) ApplicationTools::displayTask("Initializing probabilities");
	vector<const Node *> nodes = _tree -> getNodes();
	nodes.pop_back(); //remove root
	_nbNodes = nodes.size();
	_nbClasses = _rate -> getNumberOfCategories();
	_nbStates = _alphabet -> getSize();
	for(unsigned int i = 0; i < nodes.size(); i++) {
		const Node * node = nodes[i];
		double d = node -> getDistanceToFather();
		VVVdouble * _cumpxy_node = & _cumpxy[node];
		_cumpxy_node -> resize(_nbClasses);
		for(unsigned int c = 0; c < _nbClasses; c++) {
			VVdouble * _cumpxy_node_c = & (* _cumpxy_node)[c];
			_cumpxy_node_c -> resize(_nbStates);
			Mat P = _model -> getPij_t(d * _rate -> getCategory(c));
			for(unsigned int x = 0; x < _nbStates; x++) {
				Vdouble * _cumpxy_node_c_x = & (* _cumpxy_node_c)[x];
				_cumpxy_node_c_x -> resize(_nbStates);
				(* _cumpxy_node_c_x)[0] = P(x, 0);
				for(unsigned int y = 1; y < _nbStates; y++) {
					(* _cumpxy_node_c_x)[y] = (* _cumpxy_node_c_x)[y - 1] + P(x, y);
				}
			}
		}
	}
	if(verbose) ApplicationTools::displayTaskDone();
}

/******************************************************************************/

int HomogeneousSequenceSimulator::evolve(int initialState, const Node * node, int rateClass) const
{
	Vdouble * _cumpxy_node_c_x = & _cumpxy[node][rateClass][initialState];
	double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
	for(unsigned int y = 0; y < _nbStates; y++) {
		if(rand < (* _cumpxy_node_c_x)[y]) return y;
	}
}

/******************************************************************************/

void HomogeneousSequenceSimulator::multipleEvolve(Vint & initialStates, const Node * node, Vint & rateClasses, Vint & finalStates) const
{
	VVVdouble * _cumpxy_node = & _cumpxy[node];
	for(unsigned int i = 0; i < initialStates.size(); i++) {
		Vdouble * _cumpxy_node_c_x = & (* _cumpxy_node)[rateClasses[i]][initialStates[i]];
		double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
		for(unsigned int y = 0; y < _nbStates; y++) {
			if(rand < (* _cumpxy_node_c_x)[y]) {
				finalStates[i] = y;
				break;
			}
		}
	}
}

/******************************************************************************/

void HomogeneousSequenceSimulator::evolveInternal(const Node * node, int rateClass) const
{
	if(!node -> hasFather()) { cerr << "DEBUG: HomogeneousSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl; return; }
	int initialState = _states[node -> getFather()];
	//int   finalState = _process -> evolve(initialState, node -> getDistanceToFather() * _rate -> getCategory[rateClass]);
	int   finalState = evolve(initialState, node, rateClass);
	_states[node] = finalState;
	for(unsigned int i = 0; i < node -> getNumberOfSons(); i++) {
		evolveInternal(node -> getSon(i), rateClass);
	}
}

/******************************************************************************/

void HomogeneousSequenceSimulator::multipleEvolveInternal(const Node * node, Vint & rateClasses) const
{
	if(!node -> hasFather()) { cerr << "DEBUG: HomogeneousSequenceSimulator::multipleEvolveInternal. Forbidden call of method on root node." << endl; return; }
	Vint * initialStates = _multipleStates[node -> getFather()];
	unsigned int n = initialStates -> size();
	Vint * finalStates = new Vint(n);
	//double d = node -> getDistanceToFather();
	//for(unsigned int i = 0; i < n; i++) {
	//	(* finalStates)[i] = _process -> evolve((* initialStates)[i], rates[i] * d);
	//}
	multipleEvolve(*initialStates, node, rateClasses, *finalStates);
	delete _multipleStates[node];
	_multipleStates[node] = finalStates;
	for(unsigned int i = 0; i < node -> getNumberOfSons(); i++) {
		multipleEvolveInternal(node -> getSon(i), rateClasses);
	}
}


/******************************************************************************/

void HomogeneousSequenceSimulator::preciseEvolveInternal(const Node * node, int rateClass) const
{
	if(!node -> hasFather()) { cerr << "DEBUG: HomogeneousSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl; return; }
	int initialState = _states[node -> getFather()];
	MutationPath mp = _process -> detailedEvolve(initialState, node -> getDistanceToFather() * _rate -> getCategory(rateClass));
	_states[node] = mp.getFinalState();

	// Now append infos in _ssr:
	_hssr -> addNode(node, mp);

	// Now jump to son nodes:
	for(unsigned int i = 0; i < node -> getNumberOfSons(); i++) {
		preciseEvolveInternal(node -> getSon(i), rateClass);
	}
}
		
/******************************************************************************/
	
Site * HomogeneousSequenceSimulator::evolve(int initialState, int rateClass) const
{
	// Launch recursion:
	const Node * root = _tree -> getRootNode();
	_states[root] = initialState;
	for(unsigned int i = 0; i < root -> getNumberOfSons(); i++) {
		evolveInternal(root -> getSon(i), rateClass);
	}
	// Now create a Site object:
	Vint site(_leaves.size());
	for(unsigned int i = 0; i < _leaves.size(); i++) {
		site[i] = _states[_leaves[i]];
	}
	return new Site(site, _alphabet);
}

/******************************************************************************/
	
SiteContainer * HomogeneousSequenceSimulator::multipleEvolve(Vint & initialStates, Vint & rateClasses) const
{
	// Launch recursion:
	const Node * root = _tree -> getRootNode();
	_multipleStates[root] = & initialStates;
	for(unsigned int i = 0; i < root -> getNumberOfSons(); i++) {
		multipleEvolveInternal(root -> getSon(i), rateClasses);
	}
	// Now create a Site object:
	AlignedSequenceContainer * sites = new AlignedSequenceContainer(_alphabet);
	unsigned int n = _leaves.size();
	for(unsigned int i = 0; i < n; i++) {
		sites -> addSequence(Sequence(_leaves[i] -> getName(), * _multipleStates[_leaves[i]], _alphabet), false);
	}
	return sites;
}

/******************************************************************************/
	
void HomogeneousSequenceSimulator::preciseEvolve(int initialState, int rateClass) const
{
	// Launch recursion:
	const Node * root = _tree -> getRootNode();
	_states[root] = initialState;
	for(unsigned int i = 0; i < root -> getNumberOfSons(); i++) {
		preciseEvolveInternal(root -> getSon(i), rateClass);
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
	int rateClass = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(_rate -> getNumberOfCategories());
	// Make this state evolve:
	return evolve(initialState, rateClass);
}

/******************************************************************************/

#include <Seq/Mase.h>

SiteContainer * HomogeneousSequenceSimulator::simulate(unsigned int numberOfSites) const
{
	// really unefficient!
	//vector<const Site *> vs(numberOfSites);
	//for(unsigned int i = 0; i < numberOfSites; i++) {
	//	Site * s = simulate();
	//	s -> setPosition((int)i);
	//	vs[i] = s;
	//}
	//SiteContainer * sites = new VectorSiteContainer(vs, _alphabet);
	//sites -> setSequencesNames(_seqNames, false);
	// Freeing memory:
	//for(unsigned int i = 0; i < numberOfSites; i++) delete vs[i];
	//return sites;
	// Draw an initial state randomly according to equilibrum frequencies:
	Vint * initialStates = new Vint(numberOfSites, 0);
	for(unsigned int j = 0; j < numberOfSites; j++) { 
		double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
		double cumprob = 0;
		for(unsigned int i = 0; i < _alphabet -> getSize(); i++)	{
			cumprob += _model -> freq(i);
			if(r <= cumprob) {
				(* initialStates)[j] = i;
				break;
			}
		}
	}
	// Draw random rates:
	Vint rateClasses(numberOfSites);
	unsigned int nCat = _rate -> getNumberOfCategories();
	for(unsigned int j = 0; j < numberOfSites; j++) rateClasses[j] = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(nCat);
	// Make these states evolve:
	SiteContainer * sites = multipleEvolve(* initialStates, rateClasses);
	//Mase mase;
	//mase.write(string("tmp.mase"), *sites, true);
	return sites;
	//return multipleEvolve(* initialStates, rateClasses);
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
	int rateClass = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(_rate -> getNumberOfCategories());
	// Make this state evolve:
	_hssr = new HomogeneousSequenceSimulationResult(_tree, initialState, rateClass);
	preciseEvolve(initialState, rateClass);
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

const Tree<Node> * 
HomogeneousSequenceSimulator::getTree() const { return _tree; }

/******************************************************************************/
