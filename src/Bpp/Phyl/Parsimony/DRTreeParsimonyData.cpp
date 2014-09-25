//
// File: DRTreeParsimonyData.cpp
// Created by: Julien Dutheil
// Created on: Tue Jan O9 17:38 2007
// From file: DRHTreeParsimonyScore.cpp
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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
#include "../SitePatterns.h"

// From SeqLib:
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

DRTreeParsimonyData::DRTreeParsimonyData(const DRTreeParsimonyData& data) :
  AbstractTreeParsimonyData(data),
  nodeData_(data.nodeData_),
  leafData_(data.leafData_),
  rootBitsets_(data.rootBitsets_),
  rootScores_(data.rootScores_),
  shrunkData_(0),
  nbSites_(data.nbSites_),
  nbStates_(data.nbStates_),
  nbDistinctSites_(data.nbDistinctSites_)
{
  if (data.shrunkData_)
    shrunkData_ = dynamic_cast<SiteContainer*>(data.shrunkData_->clone());
  else
    shrunkData_ = 0;
}

/******************************************************************************/

DRTreeParsimonyData& DRTreeParsimonyData::operator=(const DRTreeParsimonyData& data)
{
  AbstractTreeParsimonyData::operator=(data);
  nodeData_        = data.nodeData_;
  leafData_        = data.leafData_;
  rootBitsets_     = data.rootBitsets_;
  rootScores_      = data.rootScores_;
  if (shrunkData_) delete shrunkData_;
  if (data.shrunkData_)
    shrunkData_ = dynamic_cast<SiteContainer*>(data.shrunkData_->clone());
  else
    shrunkData_ = 0;
  nbSites_         = data.nbSites_;
  nbStates_        = data.nbStates_;
  nbDistinctSites_ = data.nbDistinctSites_;
  return *this;
}

/******************************************************************************/
void DRTreeParsimonyData::init(const SiteContainer& sites, const StateMap& stateMap) throw (Exception)
{
  nbStates_         = stateMap.getNumberOfModelStates();
  nbSites_          = sites.getNumberOfSites();
  SitePatterns pattern(&sites);
  shrunkData_       = pattern.getSites();
  rootWeights_      = pattern.getWeights();
  rootPatternLinks_ = pattern.getIndices();
  nbDistinctSites_  = shrunkData_->getNumberOfSites();

  // Init data:
  // Clone data for more efficiency on sequences access:
  const SiteContainer* sequences = new AlignedSequenceContainer(*shrunkData_);
  init(getTreeP_()->getRootNode(), *sequences, stateMap);
  delete sequences;

  // Now initialize root arrays:
  rootBitsets_.resize(nbDistinctSites_);
  rootScores_.resize(nbDistinctSites_);
}

/******************************************************************************/
void DRTreeParsimonyData::init(const Node* node, const SiteContainer& sites, const StateMap& stateMap) throw (Exception)
{
  const Alphabet* alphabet = sites.getAlphabet();
  if (node->isLeaf())
  {
    const Sequence* seq;
    try
    {
      seq = &sites.getSequence(node->getName());
    }
    catch (SequenceNotFoundException& snfe)
    {
      throw SequenceNotFoundException("DRTreeParsimonyData:init(node, sites). Leaf name in tree not found in site container: ", (node->getName()));
    }
    DRTreeParsimonyLeafData* leafData    = &leafData_[node->getId()];
    vector<Bitset>* leafData_bitsets     = &leafData->getBitsetsArray();
    leafData->setNode(node);

    leafData_bitsets->resize(nbDistinctSites_);

    for (unsigned int i = 0; i < nbDistinctSites_; i++)
    {
      Bitset* leafData_bitsets_i = &(*leafData_bitsets)[i];
      for (unsigned int s = 0; s < nbStates_; s++)
      {
        // Leaves bitset are set to 1 if the char correspond to the site in the sequence,
        // otherwise value set to 0:
        int state = seq->getValue(i);
        vector<int> states = alphabet->getAlias(state);
        for (size_t j = 0; j < states.size(); j++)
        {
          if (stateMap.getAlphabetStateAsInt(s) == states[j])
            (*leafData_bitsets_i)[s].flip();
        }
      }
    }
  }
  else
  {
    DRTreeParsimonyNodeData* nodeData = &nodeData_[node->getId()];
    nodeData->setNode(node);
    nodeData->eraseNeighborArrays();

    int nbSons = static_cast<int>(node->getNumberOfSons());

    for (int n = (node->hasFather() ? -1 : 0); n < nbSons; n++)
    {
      const Node* neighbor = (*node)[n];
      vector<Bitset>* neighborData_bitsets       = &nodeData->getBitsetsArrayForNeighbor(neighbor->getId());
      vector<unsigned int>* neighborData_scores  = &nodeData->getScoresArrayForNeighbor(neighbor->getId());

      neighborData_bitsets->resize(nbDistinctSites_);
      neighborData_scores->resize(nbDistinctSites_);
    }
  }

  // We initialize each son node:
  size_t nbSonNodes = node->getNumberOfSons();
  for (unsigned int l = 0; l < nbSonNodes; l++)
  {
    // For each son node,
    init(node->getSon(l), sites, stateMap);
  }
}

/******************************************************************************/
void DRTreeParsimonyData::reInit() throw (Exception)
{
  reInit(getTreeP_()->getRootNode());
}

/******************************************************************************/
void DRTreeParsimonyData::reInit(const Node* node) throw (Exception)
{
  if (node->isLeaf())
  {
    return;
  }
  else
  {
    DRTreeParsimonyNodeData* nodeData = &nodeData_[node->getId()];
    nodeData->setNode(node);
    nodeData->eraseNeighborArrays();

    int nbSons = static_cast<int>(node->getNumberOfSons());

    for (int n = (node->hasFather() ? -1 : 0); n < nbSons; n++)
    {
      const Node* neighbor = (*node)[n];
      vector<Bitset>* neighborData_bitsets       = &nodeData->getBitsetsArrayForNeighbor(neighbor->getId());
      vector<unsigned int>* neighborData_scores  = &nodeData->getScoresArrayForNeighbor(neighbor->getId());

      neighborData_bitsets->resize(nbDistinctSites_);
      neighborData_scores->resize(nbDistinctSites_);
    }
  }

  // We initialize each son node:
  size_t nbSonNodes = node->getNumberOfSons();
  for (unsigned int l = 0; l < nbSonNodes; l++)
  {
    // For each son node,
    reInit(node->getSon(l));
  }
}

/******************************************************************************/

