// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../../PatternTools.h"
#include "DRASDRTreeLikelihoodData.h"

// From SeqLib:
#include <Bpp/Seq/SiteTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

void DRASDRTreeLikelihoodData::initLikelihoods(
    const AlignmentDataInterface& sites,
    const TransitionModelInterface& model)
{
  if (sites.getNumberOfSequences() == 1)
    throw Exception("Error, only 1 sequence!");
  if (sites.getNumberOfSequences() == 0)
    throw Exception("Error, no sequence!");
  if (sites.getAlphabet()->getAlphabetType()
      != model.getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("DRASDRTreeLikelihoodData::initLikelihoods. Data and model must have the same alphabet type.",
          sites.getAlphabet(),
          model.getAlphabet());
  alphabet_ = sites.getAlphabet();
  nbStates_ = model.getNumberOfStates();
  nbSites_  = sites.getNumberOfSites();

  SitePatterns pattern(sites);
  shrunkData_       = pattern.getSites();
  rootWeights_      = pattern.getWeights();

  rootPatternLinks_.resize(static_cast<size_t>(pattern.getIndices().size()));
  SitePatterns::IndicesType::Map(&rootPatternLinks_[0], pattern.getIndices().size()) = pattern.getIndices();
  nbDistinctSites_  = shrunkData_->getNumberOfSites();

  // Init data:
  initLikelihoods(tree_->getRootNode(), *shrunkData_, model);

  // Now initialize root likelihoods and derivatives:
  rootLikelihoods_.resize(nbDistinctSites_);
  rootLikelihoodsS_.resize(nbDistinctSites_);
  rootLikelihoodsSR_.resize(nbDistinctSites_);
  for (size_t i = 0; i < nbDistinctSites_; ++i)
  {
    VVdouble* rootLikelihoods_i_ = &rootLikelihoods_[i];
    Vdouble* rootLikelihoodsS_i_ = &rootLikelihoodsS_[i];
    rootLikelihoods_i_->resize(nbClasses_);
    rootLikelihoodsS_i_->resize(nbClasses_);
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* rootLikelihoods_i_c_ = &(*rootLikelihoods_i_)[c];
      rootLikelihoods_i_c_->resize(nbStates_);
      for (size_t x = 0; x < nbStates_; x++)
      {
        (*rootLikelihoods_i_c_)[x] = 1.;
      }
    }
  }
}

/******************************************************************************/

void DRASDRTreeLikelihoodData::initLikelihoods(
    const Node* node,
    const AlignmentDataInterface& sites,
    const TransitionModelInterface& model)
{
  if (node->isLeaf())
  {
    // Init leaves likelihoods:
    size_t posSeq;
    try
    {
      posSeq = sites.getSequencePosition(node->getName());
    }
    catch (SequenceNotFoundException& snfe)
    {
      throw SequenceNotFoundException("DRASDRTreeLikelihoodData::initlikelihoods. Leaf name in tree not found in site container: ", (node->getName()));
    }
    DRASDRTreeLikelihoodLeafData* leafData = &leafData_[node->getId()];
    VVdouble* leavesLikelihoods_leaf = &leafData->getLikelihoodArray();
    leafData->setNode(node);
    leavesLikelihoods_leaf->resize(nbDistinctSites_);
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      Vdouble* leavesLikelihoods_leaf_i = &(*leavesLikelihoods_leaf)[i];
      leavesLikelihoods_leaf_i->resize(nbStates_);
      double test = 0.;

      for (size_t s = 0; s < nbStates_; s++)
      {
        // Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
        // otherwise value set to 0:
        ( *leavesLikelihoods_leaf_i)[s] = sites.getStateValueAt(i, posSeq, model.getAlphabetStateAsInt(s));
        test += ( *leavesLikelihoods_leaf_i)[s];
      }
      if (test < 0.000001)
        std::cerr << "WARNING!!! Likelihood will be 0 for site " << i << std::endl;
    }
  }

  // We initialize each son node first:
  size_t nbSonNodes = node->getNumberOfSons();
  for (size_t l = 0; l < nbSonNodes; l++)
  {
    // For each son node,
    initLikelihoods(node->getSon(l), sites, model);
  }

  // Initialize likelihood vector:
  DRASDRTreeLikelihoodNodeData* nodeData = &nodeData_[node->getId()];
  std::map<int, VVVdouble>* likelihoods_node_ = &nodeData->getLikelihoodArrays();
  nodeData->setNode(node);

  int nbSons = static_cast<int>(node->getNumberOfSons());

  for (int n = (node->hasFather() ? -1 : 0); n < nbSons; n++)
  {
    const Node* neighbor = (*node)[n];
    VVVdouble* likelihoods_node_neighbor_ = &(*likelihoods_node_)[neighbor->getId()];

    likelihoods_node_neighbor_->resize(nbDistinctSites_);

    if (neighbor->isLeaf())
    {
      VVdouble* leavesLikelihoods_leaf_ = &leafData_[neighbor->getId()].getLikelihoodArray();
      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        Vdouble* leavesLikelihoods_leaf_i_ = &(*leavesLikelihoods_leaf_)[i];
        VVdouble* likelihoods_node_neighbor_i_ = &(*likelihoods_node_neighbor_)[i];
        likelihoods_node_neighbor_i_->resize(nbClasses_);
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* likelihoods_node_neighbor_i_c_ = &(*likelihoods_node_neighbor_i_)[c];
          likelihoods_node_neighbor_i_c_->resize(nbStates_);
          for (size_t s = 0; s < nbStates_; s++)
          {
            (*likelihoods_node_neighbor_i_c_)[s] = (*leavesLikelihoods_leaf_i_)[s];
          }
        }
      }
    }
    else
    {
      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        VVdouble* likelihoods_node_neighbor_i_ = &(*likelihoods_node_neighbor_)[i];
        likelihoods_node_neighbor_i_->resize(nbClasses_);
        for (size_t c = 0; c < nbClasses_; c++)
        {
          Vdouble* likelihoods_node_neighbor_i_c_ = &(*likelihoods_node_neighbor_i_)[c];
          likelihoods_node_neighbor_i_c_->resize(nbStates_);
          for (size_t s = 0; s < nbStates_; s++)
          {
            (*likelihoods_node_neighbor_i_c_)[s] = 1.; // All likelihoods are initialized to 1.
          }
        }
      }
    }
  }

  // Initialize d and d2 likelihoods:
  Vdouble* dLikelihoods_node_ = &nodeData->getDLikelihoodArray();
  Vdouble* d2Likelihoods_node_ = &nodeData->getD2LikelihoodArray();
  dLikelihoods_node_->resize(nbDistinctSites_);
  d2Likelihoods_node_->resize(nbDistinctSites_);
}

/******************************************************************************/

void DRASDRTreeLikelihoodData::reInit()
{
  reInit(tree_->getRootNode());
}

void DRASDRTreeLikelihoodData::reInit(const Node* node)
{
  if (node->isLeaf())
  {
    DRASDRTreeLikelihoodLeafData* leafData = &leafData_[node->getId()];
    leafData->setNode(node);
  }

  DRASDRTreeLikelihoodNodeData* nodeData = &nodeData_[node->getId()];
  nodeData->setNode(node);
  nodeData->eraseNeighborArrays();

  int nbSons = static_cast<int>(node->getNumberOfSons());

  for (int n = (node->hasFather() ? -1 : 0); n < nbSons; n++)
  {
    const Node* neighbor = (*node)[n];
    VVVdouble* array = &nodeData->getLikelihoodArrayForNeighbor(neighbor->getId());

    array->resize(nbDistinctSites_);
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble* array_i = &(*array)[i];
      array_i->resize(nbClasses_);
      for (size_t c = 0; c < nbClasses_; c++)
      {
        Vdouble* array_i_c = &(*array_i)[c];
        array_i_c->resize(nbStates_);
        for (size_t s = 0; s < nbStates_; s++)
        {
          (*array_i_c)[s] = 1.; // All likelihoods are initialized to 1.
        }
      }
    }
  }

  // We re-initialize each son node:
  size_t nbSonNodes = node->getNumberOfSons();
  for (size_t l = 0; l < nbSonNodes; l++)
  {
    // For each son node,
    reInit(node->getSon(l));
  }

  nodeData->getDLikelihoodArray().resize(nbDistinctSites_);
  nodeData->getD2LikelihoodArray().resize(nbDistinctSites_);
}

/******************************************************************************/
