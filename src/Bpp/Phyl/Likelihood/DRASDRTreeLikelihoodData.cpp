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
#include "../PatternTools.h"

// From SeqLib:
#include <Bpp/Seq/SiteTools.h>

using namespace bpp;

/******************************************************************************/

void DRASDRTreeLikelihoodData::initLikelihoods(const SiteContainer& sites, const SubstitutionModel& model) throw (Exception)
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

  SitePatterns pattern(&sites);
  if (shrunkData_)
    delete shrunkData_;
  shrunkData_       = pattern.getSites();
  rootWeights_      = pattern.getWeights();
  rootPatternLinks_ = pattern.getIndices();
  nbDistinctSites_  = shrunkData_->getNumberOfSites();

  // Init data:
  // Clone data for more efficiency on sequences access:
  const SiteContainer* sequences = new AlignedSequenceContainer(*shrunkData_);
  initLikelihoods(tree_->getRootNode(), *sequences, model);
  delete sequences;

  // Now initialize root likelihoods and derivatives:
  rootLikelihoods_.resize(nbDistinctSites_);
  rootLikelihoodsS_.resize(nbDistinctSites_);
  rootLikelihoodsSR_.resize(nbDistinctSites_);
  for (size_t i = 0; i < nbDistinctSites_; i++)
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

void DRASDRTreeLikelihoodData::initLikelihoods(const Node* node, const SiteContainer& sites, const SubstitutionModel& model) throw (Exception)
{
  if (node->isLeaf())
  {
    // Init leaves likelihoods:
    const Sequence* seq;
    try
    {
      seq = &sites.getSequence(node->getName());
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
      int state = seq->getValue(i);
      double test = 0.;
      for (size_t s = 0; s < nbStates_; s++)
      {
        // Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
        // otherwise value set to 0:
        ( *leavesLikelihoods_leaf_i)[s] = model.getInitValue(s, state);
        test += ( *leavesLikelihoods_leaf_i)[s];
      }
      if (test < 0.000001)
        std::cerr << "WARNING!!! Likelihood will be 0 for this site." << std::endl;
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

void DRASDRTreeLikelihoodData::reInit() throw (Exception)
{
  reInit(tree_->getRootNode());
}

void DRASDRTreeLikelihoodData::reInit(const Node* node) throw (Exception)
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

