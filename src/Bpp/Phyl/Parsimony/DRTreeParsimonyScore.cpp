//
// File: DRTreeParsimonyScore.cpp
// Created by: Julien Dutheil
// Created on: Thu Jul 28 17:25 2005
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

#include "DRTreeParsimonyScore.h"
#include "../PatternTools.h"
#include "../Tree/TreeTemplateTools.h" // Needed for NNIs

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/VectorTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

DRTreeParsimonyScore::DRTreeParsimonyScore(
  const Tree& tree,
  const SiteContainer& data,
  bool verbose,
  bool includeGaps)
throw (Exception) :
  AbstractTreeParsimonyScore(tree, data, verbose, includeGaps),
  parsimonyData_(new DRTreeParsimonyData(getTreeP_())),
  nbDistinctSites_()
{
  init_(data, verbose);
}

DRTreeParsimonyScore::DRTreeParsimonyScore(
  const Tree& tree,
  const SiteContainer& data,
  const StateMap* statesMap,
  bool verbose)
throw (Exception) :
  AbstractTreeParsimonyScore(tree, data, statesMap, verbose),
  parsimonyData_(new DRTreeParsimonyData(getTreeP_())),
  nbDistinctSites_()
{
  init_(data, verbose);
}

void DRTreeParsimonyScore::init_(const SiteContainer& data, bool verbose)
{
  if (verbose)
    ApplicationTools::displayTask("Initializing data structure");
  parsimonyData_->init(data, getStateMap());
  nbDistinctSites_ = parsimonyData_->getNumberOfDistinctSites();
  computeScores();
  if (verbose)
    ApplicationTools::displayTaskDone();
  if (verbose)
    ApplicationTools::displayResult("Number of distinct sites",
                                    TextTools::toString(nbDistinctSites_));
}

/******************************************************************************/

DRTreeParsimonyScore::DRTreeParsimonyScore(const DRTreeParsimonyScore& tp) :
  AbstractTreeParsimonyScore(tp),
  parsimonyData_(dynamic_cast<DRTreeParsimonyData*>(tp.parsimonyData_->clone())),
  nbDistinctSites_(tp.nbDistinctSites_)
{
  parsimonyData_->setTree(getTreeP_());
}

/******************************************************************************/

DRTreeParsimonyScore& DRTreeParsimonyScore::operator=(const DRTreeParsimonyScore& tp)
{
  AbstractTreeParsimonyScore::operator=(tp);
  parsimonyData_ = dynamic_cast<DRTreeParsimonyData*>(tp.parsimonyData_->clone());
  parsimonyData_->setTree(getTreeP_());
  nbDistinctSites_ = tp.nbDistinctSites_;
  return *this;
}

/******************************************************************************/

DRTreeParsimonyScore::~DRTreeParsimonyScore()
{
  delete parsimonyData_;
}

/******************************************************************************/
void DRTreeParsimonyScore::computeScores()
{
  computeScoresPostorder(getTreeP_()->getRootNode());
  computeScoresPreorder(getTreeP_()->getRootNode());
  computeScoresForNode(
    parsimonyData_->getNodeData(getTree().getRootId()),
    parsimonyData_->getRootBitsets(),
    parsimonyData_->getRootScores());
}

void DRTreeParsimonyScore::computeScoresPostorder(const Node* node)
{
  if (node->isLeaf()) return;
  DRTreeParsimonyNodeData* pData = &parsimonyData_->getNodeData(node->getId());
  for (unsigned int k = 0; k < node->getNumberOfSons(); k++)
  {
    const Node* son = node->getSon(k);
    computeScoresPostorder(son);
    vector<Bitset>* bitsets      = &pData->getBitsetsArrayForNeighbor(son->getId());
    vector<unsigned int>* scores = &pData->getScoresArrayForNeighbor(son->getId());
    if (son->isLeaf())
    {
      // son has no NodeData associated, must use LeafData instead
      vector<Bitset>* sonBitsets = &parsimonyData_->getLeafData(son->getId()).getBitsetsArray();
      for (unsigned int i = 0; i < sonBitsets->size(); i++)
      {
        (*bitsets)[i] = (*sonBitsets)[i];
        (*scores)[i]  = 0;
      }
    }
    else
    {
      computeScoresPostorderForNode(
        parsimonyData_->getNodeData(son->getId()),
        *bitsets,
        *scores);
    }
  }
}

void DRTreeParsimonyScore::computeScoresPostorderForNode(const DRTreeParsimonyNodeData& pData, vector<Bitset>& rBitsets, vector<unsigned int>& rScores)
{
  // First initialize the vectors from input:
  const Node* node = pData.getNode();
  const Node* source = node->getFather();
  vector<const Node*> neighbors = node->getNeighbors();
  size_t nbNeighbors = node->degree();
  vector< const vector<Bitset>*> iBitsets;
  vector< const vector<unsigned int>*> iScores;
  for (unsigned int k = 0; k < nbNeighbors; k++)
  {
    const Node* n = neighbors[k];
    if (n != source)
    {
      iBitsets.push_back(&pData.getBitsetsArrayForNeighbor(n->getId()));
      iScores.push_back(&pData.getScoresArrayForNeighbor(n->getId()));
    }
  }
  // Then call the general method on these arrays:
  computeScoresFromArrays(iBitsets, iScores, rBitsets, rScores);
}

void DRTreeParsimonyScore::computeScoresPreorder(const Node* node)
{
  if (node->getNumberOfSons() == 0) return;
  DRTreeParsimonyNodeData* pData = &parsimonyData_->getNodeData(node->getId());
  if (node->hasFather())
  {
    const Node* father = node->getFather();
    vector<Bitset>* bitsets      = &pData->getBitsetsArrayForNeighbor(father->getId());
    vector<unsigned int>* scores = &pData->getScoresArrayForNeighbor(father->getId());
    if (father->isLeaf())
    { // Means that the tree is rooted by a leaf... dunno if we must allow that! Let it be for now.
      // son has no NodeData associated, must use LeafData instead
      vector<Bitset>* sonBitsets = &parsimonyData_->getLeafData(father->getId()).getBitsetsArray();
      for (unsigned int i = 0; i < sonBitsets->size(); i++)
      {
        (*bitsets)[i] = (*sonBitsets)[i];
        (*scores)[i]  = 0;
      }
    }
    else
    {
      computeScoresPreorderForNode(
        parsimonyData_->getNodeData(father->getId()),
        node,
        *bitsets,
        *scores);
    }
  }
  // Recurse call:
  for (unsigned int k = 0; k < node->getNumberOfSons(); k++)
  {
    computeScoresPreorder(node->getSon(k));
  }
}

void DRTreeParsimonyScore::computeScoresPreorderForNode(const DRTreeParsimonyNodeData& pData, const Node* source, std::vector<Bitset>& rBitsets, std::vector<unsigned int>& rScores)
{
  // First initialize the vectors from input:
  const Node* node = pData.getNode();
  vector<const Node*> neighbors = node->getNeighbors();
  size_t nbNeighbors = node->degree();
  vector< const vector<Bitset>*> iBitsets;
  vector< const vector<unsigned int>*> iScores;
  for (unsigned int k = 0; k < nbNeighbors; k++)
  {
    const Node* n = neighbors[k];
    if (n != source)
    {
      iBitsets.push_back(&pData.getBitsetsArrayForNeighbor(n->getId()));
      iScores.push_back(&pData.getScoresArrayForNeighbor(n->getId()));
    }
  }
  // Then call the general method on these arrays:
  computeScoresFromArrays(iBitsets, iScores, rBitsets, rScores);
}

void DRTreeParsimonyScore::computeScoresForNode(const DRTreeParsimonyNodeData& pData, std::vector<Bitset>& rBitsets, std::vector<unsigned int>& rScores)
{
  const Node* node = pData.getNode();
  size_t nbNeighbors = node->degree();
  vector<const Node*> neighbors = node->getNeighbors();
  // First initialize the vectors fro input:
  vector< const vector<Bitset>*> iBitsets(nbNeighbors);
  vector< const vector<unsigned int>*> iScores(nbNeighbors);
  for (unsigned int k = 0; k < nbNeighbors; k++)
  {
    const Node* n = neighbors[k];
    iBitsets[k] =  &pData.getBitsetsArrayForNeighbor(n->getId());
    iScores [k] =  &pData.getScoresArrayForNeighbor(n->getId());
  }
  // Then call the general method on these arrays:
  computeScoresFromArrays(iBitsets, iScores, rBitsets, rScores);
}

/******************************************************************************/
unsigned int DRTreeParsimonyScore::getScore() const
{
  unsigned int score = 0;
  for (unsigned int i = 0; i < nbDistinctSites_; i++)
  {
    score += parsimonyData_->getRootScore(i) * parsimonyData_->getWeight(i);
  }
  return score;
}

/******************************************************************************/
unsigned int DRTreeParsimonyScore::getScoreForSite(size_t site) const
{
  return parsimonyData_->getRootScore(parsimonyData_->getRootArrayPosition(site));
}

/******************************************************************************/
void DRTreeParsimonyScore::computeScoresFromArrays(
  const vector< const vector<Bitset>*>& iBitsets,
  const vector< const vector<unsigned int>*>& iScores,
  vector<Bitset>& oBitsets,
  vector<unsigned int>& oScores)
{
  size_t nbPos  = oBitsets.size();
  size_t nbNodes = iBitsets.size();
  if (iScores.size() != nbNodes)
    throw Exception("DRTreeParsimonyScore::computeScores(); Error, input arrays must have the same length.");
  if (nbNodes < 1)
    throw Exception("DRTreeParsimonyScore::computeScores(); Error, input arrays must have a size >= 1.");
  const vector<Bitset>* bitsets0 = iBitsets[0];
  const vector<unsigned int>* scores0 = iScores[0];
  for (size_t i = 0; i < nbPos; i++)
  {
    oBitsets[i] = (*bitsets0)[i];
    oScores[i]  = (*scores0)[i];
  }
  for (size_t k = 1; k < nbNodes; k++)
  {
    const vector<Bitset>* bitsetsk = iBitsets[k];
    const vector<unsigned int>* scoresk = iScores[k];
    for (unsigned int i = 0; i < nbPos; i++)
    {
      Bitset bs = oBitsets[i] & (*bitsetsk)[i];
      oScores[i] += (*scoresk)[i];
      if (bs == 0)
      {
        bs = oBitsets[i] | (*bitsetsk)[i];
        oScores[i] += 1;
      }
      oBitsets[i] = bs;
    }
  }
}

/******************************************************************************/
double DRTreeParsimonyScore::testNNI(int nodeId) const throw (NodeException)
{
  const Node* son = getTreeP_()->getNode(nodeId);
  if (!son->hasFather()) throw NodePException("DRTreeParsimonyScore::testNNI(). Node 'son' must not be the root node.", son);
  const Node* parent = son->getFather();
  if (!parent->hasFather()) throw NodePException("DRTreeParsimonyScore::testNNI(). Node 'parent' must not be the root node.", parent);
  const Node* grandFather = parent->getFather();
  // From here: Bifurcation assumed.
  // In case of multifurcation, an arbitrary uncle is chosen.
  // If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
  size_t parentPosition = grandFather->getSonPosition(parent);
  const Node* uncle = grandFather->getSon(parentPosition > 1 ? parentPosition - 1 : 1 - parentPosition);

  // Retrieving arrays of interest:
  const DRTreeParsimonyNodeData* parentData = &parsimonyData_->getNodeData(parent->getId());
  const vector<Bitset>* sonBitsets = &parentData->getBitsetsArrayForNeighbor(son->getId());
  const vector<unsigned int>* sonScores  = &parentData->getScoresArrayForNeighbor(son->getId());
  vector<const Node*> parentNeighbors = TreeTemplateTools::getRemainingNeighbors(parent, grandFather, son);
  size_t nbParentNeighbors = parentNeighbors.size();
  vector< const vector<Bitset>*> parentBitsets(nbParentNeighbors);
  vector< const vector<unsigned int>*> parentScores(nbParentNeighbors);
  for (unsigned int k = 0; k < nbParentNeighbors; k++)
  {
    const Node* n = parentNeighbors[k]; // This neighbor
    parentBitsets[k] = &parentData->getBitsetsArrayForNeighbor(n->getId());
    parentScores[k] = &parentData->getScoresArrayForNeighbor(n->getId());
  }

  const DRTreeParsimonyNodeData* grandFatherData = &parsimonyData_->getNodeData(grandFather->getId());
  const vector<Bitset>* uncleBitsets = &grandFatherData->getBitsetsArrayForNeighbor(uncle->getId());
  const vector<unsigned int>* uncleScores  = &grandFatherData->getScoresArrayForNeighbor(uncle->getId());
  vector<const Node*> grandFatherNeighbors = TreeTemplateTools::getRemainingNeighbors(grandFather, parent, uncle);
  size_t nbGrandFatherNeighbors = grandFatherNeighbors.size();
  vector< const vector<Bitset>*> grandFatherBitsets(nbGrandFatherNeighbors);
  vector< const vector<unsigned int>*> grandFatherScores(nbGrandFatherNeighbors);
  for (unsigned int k = 0; k < nbGrandFatherNeighbors; k++)
  {
    const Node* n = grandFatherNeighbors[k]; // This neighbor
    grandFatherBitsets[k] = &grandFatherData->getBitsetsArrayForNeighbor(n->getId());
    grandFatherScores[k] = &grandFatherData->getScoresArrayForNeighbor(n->getId());
  }

  // Compute arrays and scores for grand-father node:
  grandFatherBitsets.push_back(sonBitsets);
  grandFatherScores.push_back(sonScores);
  // Init arrays:
  vector<Bitset> gfBitsets(sonBitsets->size()); // All arrays supposed to have the same size!
  vector<unsigned int> gfScores(sonScores->size());
  // Fill arrays:
  computeScoresFromArrays(grandFatherBitsets, grandFatherScores, gfBitsets, gfScores);

  // Now computes arrays and scores for parent node:
  parentBitsets.push_back(uncleBitsets);
  parentScores.push_back(uncleScores);
  parentBitsets.push_back(&gfBitsets);
  parentScores.push_back(&gfScores);
  // Init arrays:
  vector<Bitset> pBitsets(sonBitsets->size()); // All arrays supposed to have the same size!
  vector<unsigned int> pScores(sonScores->size());
  // Fill arrays:
  computeScoresFromArrays(parentBitsets, parentScores, pBitsets, pScores);

  // Final computation:
  unsigned int score = 0;
  for (unsigned int i = 0; i < nbDistinctSites_; i++)
  {
    score += pScores[i] * parsimonyData_->getWeight(i);
  }
  return (double)score - (double)getScore();
}

/******************************************************************************/
void DRTreeParsimonyScore::doNNI(int nodeId) throw (NodeException)
{
  Node* son = getTreeP_()->getNode(nodeId);
  if (!son->hasFather()) throw NodePException("DRTreeParsimonyScore::doNNI(). Node 'son' must not be the root node.", son);
  Node* parent = son->getFather();
  if (!parent->hasFather()) throw NodePException("DRTreeParsimonyScore::doNNI(). Node 'parent' must not be the root node.", parent);
  Node* grandFather = parent->getFather();
  // From here: Bifurcation assumed.
  // In case of multifurcation, an arbitrary uncle is chosen.
  // If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
  size_t parentPosition = grandFather->getSonPosition(parent);
  Node* uncle = grandFather->getSon(parentPosition > 1 ? parentPosition - 1 : 1 - parentPosition);
  // Swap nodes:
  parent->removeSon(son);
  grandFather->removeSon(uncle);
  parent->addSon(uncle);
  grandFather->addSon(son);
}

/******************************************************************************/

