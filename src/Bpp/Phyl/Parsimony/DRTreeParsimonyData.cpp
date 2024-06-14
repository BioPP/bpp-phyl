// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../SitePatterns.h"
#include "DRTreeParsimonyData.h"

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
  shrunkData_(nullptr),
  nbSites_(data.nbSites_),
  nbStates_(data.nbStates_),
  nbDistinctSites_(data.nbDistinctSites_)
{
  if (data.shrunkData_)
    shrunkData_.reset(data.shrunkData_->clone());
}

/******************************************************************************/

DRTreeParsimonyData& DRTreeParsimonyData::operator=(const DRTreeParsimonyData& data)
{
  AbstractTreeParsimonyData::operator=(data);
  nodeData_        = data.nodeData_;
  leafData_        = data.leafData_;
  rootBitsets_     = data.rootBitsets_;
  rootScores_      = data.rootScores_;
  if (data.shrunkData_)
    shrunkData_.reset(data.shrunkData_->clone());
  else
    shrunkData_.reset(nullptr);
  nbSites_         = data.nbSites_;
  nbStates_        = data.nbStates_;
  nbDistinctSites_ = data.nbDistinctSites_;
  return *this;
}

/******************************************************************************/

void DRTreeParsimonyData::init(
    shared_ptr<const SiteContainerInterface> sites,
    shared_ptr<const StateMapInterface> stateMap)
{
  nbStates_         = stateMap->getNumberOfModelStates();
  nbSites_          = sites->getNumberOfSites();

  SitePatterns pattern(*sites);

  auto tmp = pattern.getSites();
  shrunkData_.reset(dynamic_cast<SiteContainerInterface*>(tmp.release()));

  if (shrunkData_ == nullptr)
    throw Exception("DRTreeParsimonyData::init : Data must be plain alignments.");

  rootWeights_      = pattern.getWeights();

  rootPatternLinks_.resize(size_t(pattern.getIndices().size()));
  SitePatterns::IndicesType::Map(&rootPatternLinks_[0], pattern.getIndices().size()) = pattern.getIndices();
  nbDistinctSites_  = shrunkData_->getNumberOfSites();

  // Init data:
  // Clone data for more efficiency on sequences access:
  auto sequences = make_shared<AlignedSequenceContainer>(*shrunkData_);
  init_(tree().getRootNode(), sequences, stateMap);

  // Now initialize root arrays:
  rootBitsets_.resize(nbDistinctSites_);
  rootScores_.resize(nbDistinctSites_);
}

/******************************************************************************/

void DRTreeParsimonyData::init_(
    const Node* node,
    shared_ptr<const SiteContainerInterface> sites,
    shared_ptr<const StateMapInterface> stateMap)
{
  auto alphabet = sites->getAlphabet();
  if (node->isLeaf())
  {
    const Sequence* seq;
    try
    {
      seq = &sites->sequence(node->getName());
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
          if (stateMap->getAlphabetStateAsInt(s) == states[j])
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
    init_(node->getSon(l), sites, stateMap);
  }
}

/******************************************************************************/

void DRTreeParsimonyData::reInit()
{
  reInit_(tree().getRootNode());
}

/******************************************************************************/

void DRTreeParsimonyData::reInit_(const Node* node)
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
    reInit_(node->getSon(l));
  }
}

/******************************************************************************/
