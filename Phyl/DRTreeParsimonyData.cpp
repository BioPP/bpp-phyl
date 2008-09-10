//
// File: DRTreeParsimonyData.cpp
// Created by: Julien Dutheil
// Created on: Tue Jan O9 17:38 2007
// From file: DRHTreeParsimonyScore.cpp
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

#include "DRTreeParsimonyData.h"
#include "SitePatterns.h"

// From SeqLib:
#include <Seq/AlignedSequenceContainer.h>

using namespace bpp;

/******************************************************************************/

DRTreeParsimonyData::DRTreeParsimonyData(const DRTreeParsimonyData & data):
  AbstractTreeParsimonyData(data), 
  _nodeData(data._nodeData),
  _leafData(data._leafData),
  _rootBitsets(data._rootBitsets),
  _rootScores(data._rootScores),
  _shrunkData(NULL),
  _nbSites(data._nbSites),
  _nbStates(data._nbStates),
  _nbDistinctSites(data._nbDistinctSites)
{
  _shrunkData      = dynamic_cast<SiteContainer *>(data._shrunkData->clone());
}
    
/******************************************************************************/

DRTreeParsimonyData & DRTreeParsimonyData::operator=(const DRTreeParsimonyData & data)
{
  AbstractTreeParsimonyData::operator=(data);
  _nodeData        = data._nodeData;
  _leafData        = data._leafData;
  _rootBitsets     = data._rootBitsets;
  _rootScores      = data._rootScores;
  if(_shrunkData) delete _shrunkData;
  _shrunkData      = dynamic_cast<SiteContainer *>(data._shrunkData->clone());
  _nbSites         = data._nbSites;
  _nbStates        = data._nbStates;
  _nbDistinctSites = data._nbDistinctSites;
  return *this;
}

/******************************************************************************/

void DRTreeParsimonyData::init(const SiteContainer & sites) throw (Exception)
{
	_nbStates = sites.getAlphabet()->getSize();
 	_nbSites  = sites.getNumberOfSites();
	SitePatterns pattern(&sites);
	_shrunkData = pattern.getSites();
	_rootWeights = pattern.getWeights();
	_rootPatternLinks = pattern.getIndices();
	_nbDistinctSites = _shrunkData->getNumberOfSites();
		
	//Init data:
	// Clone data for more efficiency on sequences access:
	const SiteContainer * sequences = new AlignedSequenceContainer(* _shrunkData);
	init(_tree->getRootNode(), *sequences);
	delete sequences;

	// Now initialize root arrays:
	_rootBitsets.resize(_nbDistinctSites);
	_rootScores.resize(_nbDistinctSites);
}

/******************************************************************************/

void DRTreeParsimonyData::init(const Node * node, const SiteContainer & sites) throw (Exception)
{
	const Alphabet * alphabet = sites.getAlphabet();
	if(node->isLeaf())
  {
		const Sequence * seq;
		try {
			seq = sites.getSequence(node->getName());
		} catch (SequenceNotFoundException & snfe) {
			throw SequenceNotFoundException("DRTreeParsimonyData:init(node, sites). Leaf name in tree not found in site container: ", (node -> getName()));
		}
		DRTreeParsimonyLeafData * leafData    = & _leafData[node->getId()];
		vector<Bitset> * leafData_bitsets     = & leafData->getBitsetsArray();
		leafData->setNode(*node);
		
		leafData_bitsets->resize(_nbDistinctSites);
		
		for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
			Bitset * leafData_bitsets_i = & (* leafData_bitsets)[i];
			for(unsigned int s = 0; s < _nbStates; s++)
      {
				//Leaves bitset are set to 1 if the char correspond to the site in the sequence,
				//otherwise value set to 0:
				int state = seq->getValue(i);
				vector<int> states = alphabet->getAlias(state);
				for(unsigned int j = 0; j < states.size(); j++)
					if((int)s == states[j]) (* leafData_bitsets_i)[s].flip();
			}
		}

	}
  else
  {
		DRTreeParsimonyNodeData * nodeData = & _nodeData[node->getId()];
		nodeData->setNode(*node);
		nodeData->eraseNeighborArrays();
	
		int nbSons = node->getNumberOfSons();
	
		for(int n = (node->hasFather() ? -1 : 0); n < nbSons; n++)
    {
			const Node * neighbor = (* node)[n];
			vector<Bitset> * neighborData_bitsets       = & nodeData->getBitsetsArrayForNeighbor(neighbor->getId());
			vector<unsigned int> * neighborData_scores  = & nodeData->getScoresArrayForNeighbor(neighbor->getId());
		
			neighborData_bitsets->resize(_nbDistinctSites);
			neighborData_scores->resize(_nbDistinctSites);
		}
	}

	// We initialize each son node:
	unsigned int nbSonNodes = node->getNumberOfSons();
	for(unsigned int l = 0; l < nbSonNodes; l++)
  {
		//For each son node,
		init(node->getSon(l), sites);
	}
}

/******************************************************************************/

void DRTreeParsimonyData::reInit() throw (Exception)
{
  reInit(_tree->getRootNode());
}

/******************************************************************************/

void DRTreeParsimonyData::reInit(const Node * node) throw (Exception)
{
	if(node->isLeaf())
  {
		return;
	}
  else
  {
		DRTreeParsimonyNodeData * nodeData = & _nodeData[node->getId()];
		nodeData->setNode(*node);
		nodeData->eraseNeighborArrays();
	
		int nbSons = node->getNumberOfSons();
	
		for(int n = (node->hasFather() ? -1 : 0); n < nbSons; n++)
    {
			const Node * neighbor = (* node)[n];
			vector<Bitset> * neighborData_bitsets       = & nodeData->getBitsetsArrayForNeighbor(neighbor->getId());
			vector<unsigned int> * neighborData_scores  = & nodeData->getScoresArrayForNeighbor(neighbor->getId());
		
			neighborData_bitsets->resize(_nbDistinctSites);
			neighborData_scores->resize(_nbDistinctSites);
		}
	}

	// We initialize each son node:
	unsigned int nbSonNodes = node->getNumberOfSons();
	for(unsigned int l = 0; l < nbSonNodes; l++)
  {
		//For each son node,
		reInit(node->getSon(l));
	}
}

/******************************************************************************/

