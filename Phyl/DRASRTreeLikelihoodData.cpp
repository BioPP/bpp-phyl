//
// File: DRASRTreeLikelihoodData.cpp
// Created by: Julien Dutheil
// Created on: Sat Dec 30 14:20 2006
// From file HomogeneousTreeLikelihood.cpp
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

#include "DRASRTreeLikelihoodData.h"
#include "PatternTools.h"

//From SeqLib:
#include <Seq/SiteTools.h>
#include <Seq/AlignedSequenceContainer.h>
#include <Seq/SequenceContainerTools.h>
#include <Seq/VectorSiteContainer.h>

using namespace bpp;

/******************************************************************************/

void DRASRTreeLikelihoodData::initLikelihoods(const SiteContainer & sites, const SubstitutionModel & model) throw (Exception)
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
  if(_shrunkData != NULL) delete _shrunkData;
  SitePatterns * patterns;
  if(_usePatterns)
  {
    patterns = initLikelihoodsWithPatterns(_tree->getRootNode(), sites, model);
    _shrunkData = patterns->getSites();
    _rootWeights = patterns->getWeights();
    _rootPatternLinks = patterns->getIndices();
    _nbDistinctSites = _shrunkData->getNumberOfSites();
  }
  else
  {
    patterns = new SitePatterns(&sites);
    _shrunkData = patterns->getSites();
    _rootWeights = patterns->getWeights();
    _rootPatternLinks = patterns->getIndices();
    _nbDistinctSites = _shrunkData->getNumberOfSites();
    initLikelihoods(_tree->getRootNode(), *_shrunkData, model);
  }
  delete patterns;
}

/******************************************************************************/

void DRASRTreeLikelihoodData::initLikelihoods(const Node * node, const SiteContainer & sequences, const SubstitutionModel & model) throw (Exception)
{
  //Initialize likelihood vector:
  DRASRTreeLikelihoodNodeData * nodeData = & _nodeData[node->getId()];
  nodeData->setNode(*node);
  VVVdouble * _likelihoods_node = & nodeData->getLikelihoodArray();
  VVVdouble * _dLikelihoods_node = & nodeData->getDLikelihoodArray();
  VVVdouble * _d2Likelihoods_node = & nodeData->getD2LikelihoodArray();
  
  _likelihoods_node->resize(_nbDistinctSites);
  _dLikelihoods_node->resize(_nbDistinctSites);
  _d2Likelihoods_node->resize(_nbDistinctSites);

  for(unsigned int i = 0; i < _nbDistinctSites; i++)
  {
    VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i];
    VVdouble * _dLikelihoods_node_i = & (* _dLikelihoods_node)[i];
    VVdouble * _d2Likelihoods_node_i = & (* _d2Likelihoods_node)[i];
    _likelihoods_node_i->resize(_nbClasses);
    _dLikelihoods_node_i->resize(_nbClasses);
    _d2Likelihoods_node_i->resize(_nbClasses);
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c];
      Vdouble * _dLikelihoods_node_i_c = & (* _dLikelihoods_node_i)[c];
      Vdouble * _d2Likelihoods_node_i_c = & (* _d2Likelihoods_node_i)[c];
      _likelihoods_node_i_c->resize(_nbStates);
      _dLikelihoods_node_i_c->resize(_nbStates);
      _d2Likelihoods_node_i_c->resize(_nbStates);
      for(unsigned int s = 0; s < _nbStates; s++)
      {
        (* _likelihoods_node_i_c)[s] = 1; //All likelihoods are initialized to 1.
        (* _dLikelihoods_node_i_c)[s] = 0; //All dLikelihoods are initialized to 0.
        (* _d2Likelihoods_node_i_c)[s] = 0; //All d2Likelihoods are initialized to 0.
      }
    }
  }

  //Now initialize likelihood values and pointers:
  
  if(node->isLeaf())
  {
    const Sequence * seq;
    try
    {
      seq = sequences.getSequence(node->getName());
    }
    catch (SequenceNotFoundException snfe)
    {
      throw SequenceNotFoundException("DRASRTreeLikelihoodData::initTreelikelihoods. Leaf name in tree not found in site conainer: ", (node->getName()));
    }  
    for(unsigned int i = 0; i < _nbDistinctSites; i++)
    {
      VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i]; 
      int state = seq->getValue(i);
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c]; 
        double test = 0.;
        for(unsigned int s = 0; s < _nbStates; s++)
        {
          //Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
          //otherwise value set to 0:
          //cout << "i=" << i << "\tc=" << c << "\ts=" << s << endl;
          (* _likelihoods_node_i_c)[s] = model.getInitValue(s, state);
          test += (* _likelihoods_node_i_c)[s];
        }
        if(test < 0.000001) cerr << "WARNING!!! Likelihood will be 0 for this site." << endl;
      }
    }
  }
  else
  {
    //'node' is an internal node.
    map<int, vector<unsigned int> > * _patternLinks_node = & _patternLinks[node->getId()];
    unsigned int nbSonNodes = node->getNumberOfSons();
    for(unsigned int l = 0; l < nbSonNodes; l++)
    {
      //For each son node,
      const Node * son = (* node)[l];
      initLikelihoods(son, sequences, model);
      vector<unsigned int> * _patternLinks_node_son = & (* _patternLinks_node)[son->getId()];

      //Init map:
      _patternLinks_node_son->resize(_nbDistinctSites);

      for(unsigned int i = 0; i < _nbDistinctSites; i++)
      {
        (* _patternLinks_node_son)[i] = i;
      }
    }
  }
}

/******************************************************************************/

SitePatterns * DRASRTreeLikelihoodData::initLikelihoodsWithPatterns(const Node * node, const SiteContainer & sequences, const SubstitutionModel & model) throw (Exception)
{
  SiteContainer * tmp = PatternTools::getSequenceSubset(sequences, * node);
  SitePatterns * patterns = new SitePatterns(tmp, true);
  SiteContainer * subSequences = patterns->getSites();

  unsigned int nbSites = subSequences->getNumberOfSites();
  
  //Initialize likelihood vector:
  DRASRTreeLikelihoodNodeData * nodeData = & _nodeData[node->getId()];
  nodeData->setNode(*node);
  VVVdouble * _likelihoods_node = & nodeData->getLikelihoodArray();
  VVVdouble * _dLikelihoods_node = & nodeData->getDLikelihoodArray();
  VVVdouble * _d2Likelihoods_node = & nodeData->getD2LikelihoodArray();
  _likelihoods_node->resize(nbSites);
  _dLikelihoods_node->resize(nbSites);
  _d2Likelihoods_node->resize(nbSites);

  for(unsigned int i = 0; i < nbSites; i++)
  {
    VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i];
    VVdouble * _dLikelihoods_node_i = & (* _dLikelihoods_node)[i];
    VVdouble * _d2Likelihoods_node_i = & (* _d2Likelihoods_node)[i];
    _likelihoods_node_i->resize(_nbClasses);
    _dLikelihoods_node_i->resize(_nbClasses);
    _d2Likelihoods_node_i->resize(_nbClasses);
    for(unsigned int c = 0; c < _nbClasses; c++)
    {
      Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c];
      Vdouble * _dLikelihoods_node_i_c = & (* _dLikelihoods_node_i)[c];
      Vdouble * _d2Likelihoods_node_i_c = & (* _d2Likelihoods_node_i)[c];
      _likelihoods_node_i_c->resize(_nbStates);
      _dLikelihoods_node_i_c->resize(_nbStates);
      _d2Likelihoods_node_i_c->resize(_nbStates);
      for(unsigned int s = 0; s < _nbStates; s++)
      {
        (* _likelihoods_node_i_c)[s] = 1; //All likelihoods are initialized to 1.
        (* _dLikelihoods_node_i_c)[s] = 0; //All dLikelihoods are initialized to 0.
        (* _d2Likelihoods_node_i_c)[s] = 0; //All d2Likelihoods are initialized to 0.
      }
    }
  }

  //Now initialize likelihood values and pointers:
    
  if(node->isLeaf())
  {
    const Sequence * seq;
    try
    {
      seq = subSequences->getSequence(node->getName());
    }
    catch (SequenceNotFoundException snfe)
    {
      throw SequenceNotFoundException("HomogeneousTreeLikelihood::initTreelikelihoodsWithPatterns. Leaf name in tree not found in site conainer: ", (node->getName()));
    }  
    for(unsigned int i = 0; i < nbSites; i++)
    {
      VVdouble * _likelihoods_node_i = & (* _likelihoods_node)[i];
      int state = seq->getValue(i);
      for(unsigned int c = 0; c < _nbClasses; c++)
      {
        Vdouble * _likelihoods_node_i_c = & (* _likelihoods_node_i)[c];
        double test = 0.;
        for(unsigned int s = 0; s < _nbStates; s++)
        {
          //Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
          //otherwise value set to 0:
          //cout << "i=" << i << "\tc=" << c << "\ts=" << s << endl;
          (* _likelihoods_node_i_c)[s] = model.getInitValue(s, state);
          test += (* _likelihoods_node_i_c)[s];
        }
        if(test < 0.000001) cerr << "WARNING!!! Likelihood will be 0 for this site." << endl;
      }
    }
  }
  else
  {
    //'node' is an internal node.
    map<int, vector<unsigned int> > * _patternLinks_node = & _patternLinks[node->getId()];
    
    //Now initialize pattern links:
    unsigned int nbSonNodes = node->getNumberOfSons();
    for(unsigned int l = 0; l < nbSonNodes; l++)
    {
      //For each son node,
      const Node * son = (* node)[l];

      vector<unsigned int> * _patternLinks_node_son = & (* _patternLinks_node)[son->getId()];
      
      //Initialize subtree 'l' and retrieves corresponding subSequences:
      SitePatterns * subPatterns = initLikelihoodsWithPatterns(son, *subSequences, model);
      (* _patternLinks_node_son) = subPatterns->getIndices();
      delete subPatterns;
    }
  }
  delete subSequences;
  return patterns;
}

/******************************************************************************/

