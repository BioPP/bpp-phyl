//
// File: DRASDRTreeLikelihoodData.cpp
// Created by: Julien Dutheil
// Created on: Sat Dec 30 14:20 2006
// From file DRHomogeneousTreeLikelihood.cpp
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

#include "DRASDRTreeLikelihoodData.h"
#include "PatternTools.h"

//From SeqLib:
#include <Seq/SiteTools.h>

using namespace bpp;

/******************************************************************************/

void DRASDRTreeLikelihoodData::initLikelihoods(const SiteContainer & sites, const SubstitutionModel & model) throw (Exception)
{  
  if(sites.getNumberOfSequences() == 1) throw Exception("Error, only 1 sequence!");
  if(sites.getNumberOfSequences() == 0) throw Exception("Error, no sequence!");
  if(sites.getAlphabet()->getAlphabetType()
      != model.getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("DRASDRTreeLikelihoodData::initLikelihoods. Data and model must have the same alphabet type.",
        sites.getAlphabet(),
        model.getAlphabet());
  _alphabet = sites.getAlphabet();
  _nbStates = model.getNumberOfStates();
  _nbSites  = sites.getNumberOfSites();
  
  SitePatterns pattern(&sites);
  if(_shrunkData != NULL) delete _shrunkData;
  _shrunkData = pattern.getSites();
  _rootWeights = pattern.getWeights();
  _rootPatternLinks = pattern.getIndices();
  _nbDistinctSites = _shrunkData->getNumberOfSites();
  
  //Init data:
  // Clone data for more efficiency on sequences access:
  const SiteContainer * sequences = new AlignedSequenceContainer(* _shrunkData);
  initLikelihoods(_tree->getRootNode(), * sequences, model);
  delete sequences;

  // Now initialize root likelihoods and derivatives:
  _rootLikelihoods.resize(_nbDistinctSites);
  _rootLikelihoodsS.resize(_nbDistinctSites);
  _rootLikelihoodsSR.resize(_nbDistinctSites);
  for(unsigned int i = 0; i < _nbDistinctSites; i++)
  {
    VVdouble * _rootLikelihoods_i = & _rootLikelihoods[i];
    Vdouble * _rootLikelihoodsS_i = & _rootLikelihoodsS[i];
    _rootLikelihoods_i->resize(_nbClasses);
    _rootLikelihoodsS_i->resize(_nbClasses);
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      Vdouble * _rootLikelihoods_i_c = & (* _rootLikelihoods_i)[c];
      _rootLikelihoods_i_c->resize(_nbStates);
      for(unsigned int x = 0; x < _nbStates; x++)
      {
        (* _rootLikelihoods_i_c)[x] = 1.;
      }
    }
  }
}

/******************************************************************************/

void DRASDRTreeLikelihoodData::initLikelihoods(const Node * node, const SiteContainer & sites, const SubstitutionModel & model) throw (Exception)
{
  if(node->isLeaf())
  {
    // Init leaves likelihoods:
    const Sequence * seq;
    try
    {
      seq = sites.getSequence(node->getName());
    }
    catch (SequenceNotFoundException & snfe)
    {
      throw SequenceNotFoundException("DRASDRTreeLikelihoodData::initlikelihoods. Leaf name in tree not found in site container: ", (node->getName()));
    }
    DRASDRTreeLikelihoodLeafData * leafData = & _leafData[node->getId()];
    VVdouble * leavesLikelihoods_leaf = & leafData->getLikelihoodArray();
    leafData->setNode(*node);
    leavesLikelihoods_leaf->resize(_nbDistinctSites);
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      Vdouble * leavesLikelihoods_leaf_i = & (* leavesLikelihoods_leaf)[i];
      leavesLikelihoods_leaf_i->resize(_nbStates);
      int state = seq->getValue(i);
      double test = 0.;
      for(unsigned int s = 0; s < _nbStates; s++)
      {
        //Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
        //otherwise value set to 0:
        ( * leavesLikelihoods_leaf_i)[s] = model.getInitValue(s, state);
        test += ( * leavesLikelihoods_leaf_i)[s];
      }
      if(test < 0.000001) cerr << "WARNING!!! Likelihood will be 0 for this site." << endl;
    }
  }

  // We initialize each son node first:
  unsigned int nbSonNodes = node->getNumberOfSons();
  for(unsigned int l = 0; l < nbSonNodes; l++)
  {
    //For each son node,
    initLikelihoods(node->getSon(l), sites, model);
  }

  //Initialize likelihood vector:
  DRASDRTreeLikelihoodNodeData * nodeData = & _nodeData[node->getId()];
  map<int, VVVdouble> * _likelihoods_node = & nodeData->getLikelihoodArrays();
  nodeData->setNode(*node);
  
  int nbSons = node->getNumberOfSons();
  
  for(int n = (node->hasFather() ? -1 : 0); n < nbSons; n++)
  {
    const Node * neighbor = (* node)[n];
    VVVdouble * _likelihoods_node_neighbor = & (* _likelihoods_node)[neighbor->getId()];
    
    _likelihoods_node_neighbor->resize(_nbDistinctSites);

    if(neighbor->isLeaf())
    {
      VVdouble * _leavesLikelihoods_leaf = & _leafData[neighbor->getId()].getLikelihoodArray();
      for(unsigned int i = 0; i < _nbDistinctSites; i++)
      {
        Vdouble  * _leavesLikelihoods_leaf_i = & (* _leavesLikelihoods_leaf)[i];
        VVdouble * _likelihoods_node_neighbor_i = & (* _likelihoods_node_neighbor)[i];
        _likelihoods_node_neighbor_i->resize(_nbClasses);
        for(unsigned int c = 0; c < _nbClasses; c++)
        {
          Vdouble * _likelihoods_node_neighbor_i_c = & (* _likelihoods_node_neighbor_i)[c];
          _likelihoods_node_neighbor_i_c->resize(_nbStates);
          for(unsigned int s = 0; s < _nbStates; s++)
          {
            (* _likelihoods_node_neighbor_i_c)[s] = (* _leavesLikelihoods_leaf_i)[s];
          }
        }
      }
    }
    else
    {
      for(unsigned int i = 0; i < _nbDistinctSites; i++)
      {
        VVdouble * _likelihoods_node_neighbor_i = & (* _likelihoods_node_neighbor)[i];
        _likelihoods_node_neighbor_i->resize(_nbClasses);
        for(unsigned int c = 0; c < _nbClasses; c++)
        {
          Vdouble * _likelihoods_node_neighbor_i_c = & (* _likelihoods_node_neighbor_i)[c];
          _likelihoods_node_neighbor_i_c->resize(_nbStates);
          for(unsigned int s = 0; s < _nbStates; s++)
          {
            (* _likelihoods_node_neighbor_i_c)[s] = 1.; //All likelihoods are initialized to 1.
          }
        }
      }
    }
  }

  // Initialize d and d2 likelihoods:
  Vdouble * _dLikelihoods_node = & nodeData->getDLikelihoodArray();
  Vdouble * _d2Likelihoods_node = & nodeData->getD2LikelihoodArray();
  _dLikelihoods_node->resize(_nbDistinctSites);
  _d2Likelihoods_node->resize(_nbDistinctSites);   
}

/******************************************************************************/

void DRASDRTreeLikelihoodData::reInit() throw (Exception)
{
  reInit(_tree->getRootNode());
}

void DRASDRTreeLikelihoodData::reInit(const Node * node) throw (Exception)
{
	if(node->isLeaf())
  {
    DRASDRTreeLikelihoodLeafData * leafData = & _leafData[node->getId()];
	  leafData->setNode(*node);
  }

  DRASDRTreeLikelihoodNodeData * nodeData = & _nodeData[node->getId()];
	nodeData->setNode(*node);
	nodeData->eraseNeighborArrays();
	
	int nbSons = node->getNumberOfSons();
	
	for(int n = (node->hasFather() ? -1 : 0); n < nbSons; n++)
  {
		const Node * neighbor = (* node)[n];
		VVVdouble *array = & nodeData->getLikelihoodArrayForNeighbor(neighbor->getId());
		
		array->resize(_nbDistinctSites);
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      VVdouble * array_i = & (* array)[i];
      array_i->resize(_nbClasses);
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        Vdouble * array_i_c = & (* array_i)[c];
        array_i_c->resize(_nbStates);
        for(unsigned int s = 0; s < _nbStates; s++)
        {
          (* array_i_c)[s] = 1.; //All likelihoods are initialized to 1.
        }
      }
    }
	}

	// We re-initialize each son node:
	unsigned int nbSonNodes = node->getNumberOfSons();
	for(unsigned int l = 0; l < nbSonNodes; l++)
  {
		//For each son node,
		reInit(node->getSon(l));
	}

  nodeData->getDLikelihoodArray().resize(_nbDistinctSites);
  nodeData->getD2LikelihoodArray().resize(_nbDistinctSites);
}

/******************************************************************************/

