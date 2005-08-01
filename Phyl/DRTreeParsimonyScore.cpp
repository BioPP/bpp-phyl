//
// File: DRTreeParsimonyScore.cpp
// Created by: Julien Dutheil
// Created on: Thu Jul 28 17:25 2005
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

#include "DRTreeParsimonyScore.h"
#include "ApplicationTools.h"
#include "PatternTools.h"

/******************************************************************************/

void DRTreeParsimonyData::init(const SiteContainer & sites) throw (Exception)
{
	_nbStates = sites.getAlphabet() -> getSize();
 	_nbSites  = sites.getNumberOfSites();
	Pattern pattern = PatternTools::countSites(sites);
	_shrunkData = PatternTools::getSites(pattern, sites.getAlphabet());
	_rootWeights = PatternTools::getWeights(pattern);
	_rootPatternLinks = PatternTools::getIndices(pattern);
	_nbDistinctSites = _shrunkData -> getNumberOfSites();
		
	//Init data:
	// Clone data for more efficiency on sequences access:
	const SiteContainer * sequences = new AlignedSequenceContainer(* _shrunkData);
	init(_tree -> getRootNode(), *sequences);
	delete sequences;

	// Now initialize root arrays:
	_rootBitsets.resize(_nbDistinctSites);
	_rootScores.resize(_nbDistinctSites);
}

void DRTreeParsimonyData::init(const Node * node, const SiteContainer & sites) throw (Exception)
{
	const Alphabet * alphabet = sites.getAlphabet();
	if(node -> isLeaf()) {
		DRTreeParsimonyLeafData * leafData    = & _leafData[node];
		vector<Bitset> * leafData_bitsets     = & leafData -> getBitsetsArray();
		leafData -> setNode(node);
		
		leafData_bitsets -> resize(_nbDistinctSites);
		
		for(unsigned int i = 0; i < _nbDistinctSites; i++) {
			Bitset * leafData_bitsets_i = & (* leafData_bitsets)[i];
			for(unsigned int s = 0; s < _nbStates; s++) {
				//Leaves bitset are set to 1 if the char correspond to the site in the sequence,
				//otherwise value set to 0:
				try {
					int state = sites.getSequence(node -> getName()) -> getValue(i);
					vector<int> states = alphabet -> getAlias(state);
					for(unsigned int j = 0; j < states.size(); j++)
						if(s == states[j]) (* leafData_bitsets_i)[s].flip();
				} catch (SequenceNotFoundException & snfe) {
					throw SequenceNotFoundException("DRTreeParsimonyData:init(node, sites). Leaf name in tree not found in site container: ", (node -> getName()));
				}
			}
		}

	} else {
		DRTreeParsimonyNodeData * nodeData = & _nodeData[node];
		nodeData -> setNode(node);
	
		int nbSons = node -> getNumberOfSons();
	
		for(int n = (node -> hasFather() ? -1 : 0); n < nbSons; n++) {
			const Node * neighbor = (* node)[n];
			vector<Bitset> * neighborData_bitsets       = & nodeData -> getBitsetsArrayForNeighbor(neighbor);
			vector<unsigned int> * neighborData_scores  = & nodeData -> getScoresArrayForNeighbor(neighbor);
		
			neighborData_bitsets -> resize(_nbDistinctSites);
			neighborData_scores  -> resize(_nbDistinctSites);
		}
	}

	// We initialize each son node:
	unsigned int nbSonNodes = node -> getNumberOfSons();
	for(unsigned int l = 0; l < nbSonNodes; l++) {
		//For each son node,
		init(node -> getSon(l), sites);
	}
}

/******************************************************************************/

DRTreeParsimonyScore::DRTreeParsimonyScore(
		TreeTemplate<Node> & tree,
		const SiteContainer & data,
		bool verbose)
	throw (Exception) :
	AbstractTreeParsimonyScore(tree, data, verbose)
{
	if(verbose) ApplicationTools::message << "Double-Recursive Tree Parsimony Score" << endl;	
	_parsimonyData = new DRTreeParsimonyData(*_tree);
	
	if(verbose) ApplicationTools::displayTask("Initializing data structure");
	_parsimonyData -> init(data);
	_nbDistinctSites = _parsimonyData -> getNumberOfDistinctSites();
	_rootScores.resize(_nbDistinctSites);
	_rootBitsets.resize(_nbDistinctSites);
	computeScores();
	if(verbose) ApplicationTools::displayTaskDone();
	if(verbose) ApplicationTools::displayResult("Number of distinct sites",
			TextTools::toString(_nbDistinctSites));
}

/******************************************************************************/

DRTreeParsimonyScore::~DRTreeParsimonyScore() { delete _parsimonyData; }

/******************************************************************************/

void DRTreeParsimonyScore::computeScores()
{
	computeScoresPostorder(_tree -> getRootNode());
	computeScoresPreorder(_tree -> getRootNode());
	computeScoresForNode(_parsimonyData -> getNodeData(_tree -> getRootNode()), _rootBitsets, _rootScores);
}

void DRTreeParsimonyScore::computeScoresPostorder(const Node * node)
{
	if(node -> isLeaf()) return;
	DRTreeParsimonyNodeData * pData = & _parsimonyData -> getNodeData(node);
	for(unsigned int k = 0; k < node -> getNumberOfSons(); k++) {
		const Node * son = node -> getSon(k);
		computeScoresPostorder(son);
		vector<Bitset> * bitsets      = & pData -> getBitsetsArrayForNeighbor(son);
		vector<unsigned int> * scores = & pData -> getScoresArrayForNeighbor(son);		
		if(son -> isLeaf()) {
			// son has no NodeData associated, must use LeafData instead
			vector<Bitset> * sonBitsets = & _parsimonyData -> getLeafData(son).getBitsetsArray();
			for(unsigned int i = 0; i < sonBitsets -> size(); i++) {
				(*bitsets)[i] = (*sonBitsets)[i];
				(*scores)[i]  = 0; 
			}
		} else {
			computeScoresPostorderForNode(
					_parsimonyData -> getNodeData(son), 
					*bitsets,
					*scores);		
		}
	}
}

void DRTreeParsimonyScore::computeScoresPostorderForNode(const DRTreeParsimonyNodeData & pData, vector<Bitset> & rBitsets, vector<unsigned int> & rScores)
{
	const Node * node   = pData.getNode();
	unsigned int nbSons = node -> getNumberOfSons();
	if(nbSons != 2) throw Exception("DRTreeParsimonyScore::computeScoresPostorderForNode. Unimplemented for multifurcations.");
	unsigned int nbPos  = rBitsets.size();
	
	const Node * node1 = node -> getSon(0);
	const Node * node2 = node -> getSon(1);
	const vector<Bitset> * son1Bitsets = & pData.getBitsetsArrayForNeighbor(node1);
	const vector<Bitset> * son2Bitsets = & pData.getBitsetsArrayForNeighbor(node2);
	const vector<unsigned int> * son1Scores = & pData.getScoresArrayForNeighbor(node1);
	const vector<unsigned int> * son2Scores = & pData.getScoresArrayForNeighbor(node2);
	for(unsigned int i = 0; i < nbPos; i++) {
		rBitsets[i] = (*son1Bitsets)[i] & (*son2Bitsets)[i];
		rScores[i]  = (*son1Scores)[i] + (*son2Scores)[i]; 
		if(rBitsets[i] == 0) {
			rBitsets[i] = (*son1Bitsets)[i] | (*son2Bitsets)[i];
			rScores[i]++;
		}
	}
}

void DRTreeParsimonyScore::computeScoresPreorder(const Node * node)
{
	if(!node -> hasFather() || node -> getNumberOfSons() == 0) return;
	DRTreeParsimonyNodeData * pData = & _parsimonyData -> getNodeData(node);
	const Node * father = node -> getFather();
	vector<Bitset> * bitsets      = &pData -> getBitsetsArrayForNeighbor(father);
	vector<unsigned int> * scores = &pData -> getScoresArrayForNeighbor(father);		
	if(father -> isLeaf()) { // Means that the tree is rooted by a leaf... dunno if we must allow that! Let it be for now.
		// son has no NodeData associated, must use LeafData instead
		vector<Bitset> * sonBitsets = & _parsimonyData -> getLeafData(father).getBitsetsArray();
		for(unsigned int i = 0; i < sonBitsets -> size(); i++) {
			(*bitsets)[i] = (*sonBitsets)[i];
			(*scores)[i]  = 0; 
		}
	} else {
		computeScoresPreorderForNode(
				_parsimonyData -> getNodeData(father),
				node,
				*bitsets,
				*scores);		
	}
	for(unsigned int k = 0; k < node -> getNumberOfSons(); k++) computeScoresPreorder(node -> getSon(k));
}

void DRTreeParsimonyScore::computeScoresPreorderForNode(const DRTreeParsimonyNodeData & pData, const Node * source, vector<Bitset> & rBitsets, vector<unsigned int> & rScores)
{
	const Node * node   = pData.getNode();
	unsigned int nbSons = node -> getNumberOfSons();
	if(nbSons != 2) throw Exception("DRTreeParsimonyScore::computeScoresPreorderForNode. Unimplemented for multifurcations.");
	unsigned int nbPos  = rBitsets.size();

	const Node * node1 = node -> getSon(0);
	const Node * node2 = node -> getSon(1);
	if(node1 == source) node1 = node -> getFather();
	else                node2 = node -> getFather(); // node2 == source
	const vector<Bitset> * son1Bitsets = & pData.getBitsetsArrayForNeighbor(node1);
	const vector<Bitset> * son2Bitsets = & pData.getBitsetsArrayForNeighbor(node2);
	const vector<unsigned int> * son1Scores = & pData.getScoresArrayForNeighbor(node1);
	const vector<unsigned int> * son2Scores = & pData.getScoresArrayForNeighbor(node2);
	for(unsigned int i = 0; i < nbPos; i++) {
		rBitsets[i] = (*son1Bitsets)[i] & (*son2Bitsets)[i];
		rScores[i]  = (*son1Scores)[i] + (*son2Scores)[i]; 
		if(rBitsets[i] == 0) {
			rBitsets[i] = (*son1Bitsets)[i] | (*son2Bitsets)[i];
			rScores[i]++;
		}
	}
}

void DRTreeParsimonyScore::computeScoresForNode(const DRTreeParsimonyNodeData & pData, vector<Bitset> & rBitsets, vector<unsigned int> & rScores)
{
	const Node * node   = pData.getNode();
	unsigned int nbNeighbors = node -> degree();
	if(nbNeighbors != 3) throw Exception("DRTreeParsimonyScore::computeScoresForNode. Unimplemented for multifurcations.");
	vector<const Node *> neighbors = node -> getNeighbors();
	unsigned int nbPos  = rBitsets.size();
	
	const Node * node1 = neighbors[0];
	const Node * node2 = neighbors[1];
	const Node * node3 = neighbors[2];
	const vector<Bitset> * son1Bitsets = & pData.getBitsetsArrayForNeighbor(node1);
	const vector<Bitset> * son2Bitsets = & pData.getBitsetsArrayForNeighbor(node2);
	const vector<Bitset> * son3Bitsets = & pData.getBitsetsArrayForNeighbor(node3);
	const vector<unsigned int> * son1Scores = & pData.getScoresArrayForNeighbor(node1);
	const vector<unsigned int> * son2Scores = & pData.getScoresArrayForNeighbor(node2);
	const vector<unsigned int> * son3Scores = & pData.getScoresArrayForNeighbor(node3);
	for(unsigned int i = 0; i < nbPos; i++) {
		rBitsets[i] = (*son1Bitsets)[i] & (*son2Bitsets)[i];
		rScores[i]  = (*son1Scores)[i] + (*son2Scores)[i]; 
		if(rBitsets[i] == 0) {
			rBitsets[i] = (*son1Bitsets)[i] | (*son2Bitsets)[i];
			rScores[i] += 1;
		}
		rBitsets[i] &= (*son3Bitsets)[i];
		rScores[i] += (*son3Scores)[i]; 
		if(rBitsets[i] == 0) {
			rBitsets[i] |= (*son3Bitsets)[i];
			rScores[i] += 1;
		}
	}
}

/******************************************************************************/

unsigned int DRTreeParsimonyScore::getScore() const
{
	unsigned int score = 0;
	for(unsigned int i = 0; i < _nbDistinctSites; i++) {
		score += _rootScores[i] * _parsimonyData -> getWeight(i);
	}
	return score;
}

/******************************************************************************/

unsigned int DRTreeParsimonyScore::getScoreForSite(unsigned int site) const
{
	return _rootScores[_parsimonyData -> getRootArrayPosition(site)];
}
	
/******************************************************************************/

