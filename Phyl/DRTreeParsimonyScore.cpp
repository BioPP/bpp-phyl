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
#include "PatternTools.h"
#include "SitePatterns.h"
#include "TreeTools.h" //Needed for NNIs

// From Utils:
#include <Utils/ApplicationTools.h>

// From NumCalc:
#include <NumCalc/VectorTools.h>
using namespace VectorFunctions;

// From SeqLib:
#include <Seq/AlignedSequenceContainer.h>

/******************************************************************************/

void DRTreeParsimonyData::init(const SiteContainer & sites) throw (Exception)
{
	_nbStates = sites.getAlphabet() -> getSize();
 	_nbSites  = sites.getNumberOfSites();
	SitePatterns pattern(sites);
	_shrunkData = pattern.getSites();
	_rootWeights = pattern.getWeights();
	_rootPatternLinks = pattern.getIndices();
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
		const Sequence * seq;
		try {
			seq = sites.getSequence(node -> getName());
		} catch (SequenceNotFoundException & snfe) {
			throw SequenceNotFoundException("DRTreeParsimonyData:init(node, sites). Leaf name in tree not found in site container: ", (node -> getName()));
		}
		DRTreeParsimonyLeafData * leafData    = & _leafData[node];
		vector<Bitset> * leafData_bitsets     = & leafData -> getBitsetsArray();
		leafData -> setNode(node);
		
		leafData_bitsets -> resize(_nbDistinctSites);
		
		for(unsigned int i = 0; i < _nbDistinctSites; i++) {
			Bitset * leafData_bitsets_i = & (* leafData_bitsets)[i];
			for(unsigned int s = 0; s < _nbStates; s++) {
				//Leaves bitset are set to 1 if the char correspond to the site in the sequence,
				//otherwise value set to 0:
				int state = seq -> getValue(i);
				vector<int> states = alphabet -> getAlias(state);
				for(unsigned int j = 0; j < states.size(); j++)
					if((int)s == states[j]) (* leafData_bitsets_i)[s].flip();
			}
		}

	} else {
		DRTreeParsimonyNodeData * nodeData = & _nodeData[node];
		nodeData -> setNode(node);
		nodeData -> eraseNeighborArrays();
	
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

void DRTreeParsimonyData::reInit() throw (Exception)
{
  reInit(_tree -> getRootNode());
}

void DRTreeParsimonyData::reInit(const Node * node) throw (Exception)
{
	if(node -> isLeaf()) {
		return;
	} else {
		DRTreeParsimonyNodeData * nodeData = & _nodeData[node];
		nodeData -> setNode(node);
		nodeData -> eraseNeighborArrays();
	
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
		reInit(node -> getSon(l));
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
	//First initialize the vectors from input:
	const Node * node = pData.getNode();
	const Node * source = node -> getFather();
	vector<const Node *> neighbors = node -> getNeighbors();
	unsigned int nbNeighbors = node -> degree();
	vector< const vector<Bitset>       *> iBitsets; 
	vector< const vector<unsigned int> *> iScores;
	for(unsigned int k = 0; k < nbNeighbors; k++) {
		const Node * n = neighbors[k];
		if(n != source) {
			iBitsets.push_back(& pData.getBitsetsArrayForNeighbor(n));
			iScores.push_back(& pData.getScoresArrayForNeighbor(n));
		}
	}
	//Then call the general method on these arrays:
	computeScoresFromArrays(iBitsets, iScores, rBitsets, rScores);
}

void DRTreeParsimonyScore::computeScoresPreorder(const Node * node)
{
	if(node -> getNumberOfSons() == 0) return;
	DRTreeParsimonyNodeData * pData = & _parsimonyData -> getNodeData(node);
	if(node -> hasFather()) {
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
	}
	// Recurse call:
	for(unsigned int k = 0; k < node -> getNumberOfSons(); k++) computeScoresPreorder(node -> getSon(k));
}

void DRTreeParsimonyScore::computeScoresPreorderForNode(const DRTreeParsimonyNodeData & pData, const Node * source, vector<Bitset> & rBitsets, vector<unsigned int> & rScores)
{
	//First initialize the vectors from input:
	const Node * node = pData.getNode();
	vector<const Node *> neighbors = node -> getNeighbors();
	unsigned int nbNeighbors = node -> degree();
	vector< const vector<Bitset>       *> iBitsets; 
	vector< const vector<unsigned int> *> iScores;
	for(unsigned int k = 0; k < nbNeighbors; k++) {
		const Node * n = neighbors[k];
		if(n != source) {
			iBitsets.push_back(& pData.getBitsetsArrayForNeighbor(n));
			iScores.push_back(& pData.getScoresArrayForNeighbor(n));
		}
	}
	//Then call the general method on these arrays:
	computeScoresFromArrays(iBitsets, iScores, rBitsets, rScores);
}

void DRTreeParsimonyScore::computeScoresForNode(const DRTreeParsimonyNodeData & pData, vector<Bitset> & rBitsets, vector<unsigned int> & rScores)
{
	const Node * node   = pData.getNode();
	unsigned int nbNeighbors = node -> degree();
	vector<const Node *> neighbors = node -> getNeighbors();
	//First initialize the vectors fro input:
	vector< const vector<Bitset>       *> iBitsets(nbNeighbors); 
	vector< const vector<unsigned int> *> iScores(nbNeighbors);
	for(unsigned int k = 0; k < nbNeighbors; k++) {
		const Node * n = neighbors[k];
		iBitsets[k] =  & pData.getBitsetsArrayForNeighbor(n);
		iScores [k] =  & pData.getScoresArrayForNeighbor(n);
	}
	//Then call the general method on these arrays:
	computeScoresFromArrays(iBitsets, iScores, rBitsets, rScores);
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

void DRTreeParsimonyScore::computeScoresFromArrays(
				const vector< const vector<Bitset>       *> & iBitsets,
				const vector< const vector<unsigned int> *> & iScores,
				vector<Bitset> & oBitsets,
				vector<unsigned int> & oScores)
{
	unsigned int nbPos  = oBitsets.size();
	unsigned int nbNodes = iBitsets.size();
	if(iScores.size() != nbNodes)
		throw Exception("DRTreeParsimonyScore::computeScores(); Error, input arrays must have the same length.");
	if(nbNodes < 1)
		throw Exception("DRTreeParsimonyScore::computeScores(); Error, input arrays must have a size >= 1.");
	const vector<Bitset> * bitsets0 = iBitsets[0];
	const vector<unsigned int> * scores0 = iScores[0];
	for(unsigned int i = 0; i < nbPos; i++) {
		oBitsets[i] = (*bitsets0)[i];
		oScores[i]  = (*scores0)[i]; 
	}
	for(unsigned int k = 1; k < nbNodes; k++) {
		const vector<Bitset> * bitsetsk = iBitsets[k];
		const vector<unsigned int> * scoresk = iScores[k];
		for(unsigned int i = 0; i < nbPos; i++) {
			Bitset bs = oBitsets[i] & (*bitsetsk)[i];
			oScores[i] += (*scoresk)[i]; 
			if(bs == 0) {
				bs = oBitsets[i] | (*bitsetsk)[i];
				oScores[i] += 1;
			}
			oBitsets[i] = bs;
		}
	}
}

/******************************************************************************/

double DRTreeParsimonyScore::testNNI(const Node * parent, const Node * son) const throw (NodeException)
{
	if(!parent->hasFather())       throw NodeException("DRTreeParsimonyScore::testNNI(). Node 'parent' must not be the root node.", parent);
	if(!son   ->hasFather())       throw NodeException("DRTreeParsimonyScore::testNNI(). Node 'son' must not be the root node.", parent);
	if(son->getFather() != parent) throw NodeException("DRTreeParsimonyScore::testNNI(). Node 'son' must be a son of node 'parent'.", son);
	
	const Node * grandFather = parent->getFather();
	//From here: Bifurcation assumed.
	//In case of multifurcation, an arbitrary uncle is chosen.
	//If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
	unsigned int parentPosition = grandFather->getSonPosition(*parent);
	const Node * uncle = grandFather->getSon(parentPosition > 1 ? parentPosition -1 : 1 - parentPosition);
	
	//Retrieving arrays of interest:
	const DRTreeParsimonyNodeData * parentData = & _parsimonyData->getNodeData(parent);
	const vector<Bitset>          * sonBitsets = & parentData->getBitsetsArrayForNeighbor(son); 
	const vector<unsigned int>    * sonScores  = & parentData->getScoresArrayForNeighbor(son); 
	vector<const Node *> parentNeighbors = TreeTools::getRemainingNeighbors(parent, grandFather, son);
	unsigned int nbParentNeighbors = parentNeighbors.size();
	vector< const vector<Bitset>       *> parentBitsets(nbParentNeighbors);
	vector< const vector<unsigned int> *> parentScores(nbParentNeighbors);
	for(unsigned int k = 0; k < nbParentNeighbors; k++) {
		const Node * n = parentNeighbors[k]; // This neighbor
		parentBitsets[k] = & parentData->getBitsetsArrayForNeighbor(n); 
		parentScores[k] = & parentData->getScoresArrayForNeighbor(n); 
	}
	
	const DRTreeParsimonyNodeData * grandFatherData = & _parsimonyData->getNodeData(grandFather);
	const vector<Bitset>          * uncleBitsets = & grandFatherData->getBitsetsArrayForNeighbor(uncle); 
	const vector<unsigned int>    * uncleScores  = & grandFatherData->getScoresArrayForNeighbor(uncle); 
	vector<const Node *> grandFatherNeighbors = TreeTools::getRemainingNeighbors(grandFather, parent, uncle);
	unsigned int nbGrandFatherNeighbors = grandFatherNeighbors.size();
	vector< const vector<Bitset>       *> grandFatherBitsets(nbGrandFatherNeighbors);
	vector< const vector<unsigned int> *> grandFatherScores(nbGrandFatherNeighbors);
	for(unsigned int k = 0; k < nbGrandFatherNeighbors; k++) {
		const Node * n = grandFatherNeighbors[k]; // This neighbor
		grandFatherBitsets[k] = & grandFatherData->getBitsetsArrayForNeighbor(n); 
		grandFatherScores[k] = & grandFatherData->getScoresArrayForNeighbor(n); 
	}
	
	//Compute arrays and scores for grand-father node:
	grandFatherBitsets.push_back(sonBitsets);
	grandFatherScores.push_back(sonScores);
	//Init arrays:
	vector<Bitset> gfBitsets(sonBitsets->size());//All arrays supposed to have the same size!
	vector<unsigned int> gfScores(sonScores->size());
	//Fill arrays:
	computeScoresFromArrays(grandFatherBitsets, grandFatherScores, gfBitsets, gfScores);

	//Now computes arrays and scores for parent node:
	parentBitsets.push_back(uncleBitsets);
	parentScores.push_back(uncleScores);
	parentBitsets.push_back(&gfBitsets);
	parentScores.push_back(&gfScores);
	//Init arrays:
	vector<Bitset> pBitsets(sonBitsets->size());//All arrays supposed to have the same size!
	vector<unsigned int> pScores(sonScores->size());
	//Fill arrays:
	computeScoresFromArrays(parentBitsets, parentScores, pBitsets, pScores);

	//Final computation:
	unsigned int score = 0;
	for(unsigned int i = 0; i < _nbDistinctSites; i++) {
		score += pScores[i] * _parsimonyData -> getWeight(i);
	}
	return (double)score - (double)getScore();
}

/******************************************************************************/

void DRTreeParsimonyScore::doNNI(Node * parent, Node * son) throw (NodeException)
{
	if(!parent->hasFather())       throw NodeException("DRTreeParsimonyScore::doNNI(). Node 'parent' must not be the root node.", parent);
	if(!son   ->hasFather())       throw NodeException("DRTreeParsimonyScore::doNNI(). Node 'son' must not be the root node.", parent);
	if(son->getFather() != parent) throw NodeException("DRTreeParsimonyScore::doNNI(). Node 'son' must be a son of node 'parent'.", son);
	Node * grandFather = parent->getFather();
	//From here: Bifurcation assumed.
	//In case of multifurcation, an arbitrary uncle is chosen.
	//If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
	unsigned int parentPosition = grandFather->getSonPosition(*parent);
	Node * uncle = grandFather->getSon(parentPosition > 1 ? parentPosition -1 : 1 - parentPosition);
	//Swap nodes:
	parent->removeSon(*son);
	grandFather->removeSon(*uncle);
	parent->addSon(*uncle);
	grandFather->addSon(*son);
}

/******************************************************************************/

