//
// File: DoubleRecursiveTreeLikelihoodData.cpp
// Created by: Julien Dutheil
// Created on: Sat Dec 30 14:20 2006
// From file HomogeneousTreeLikelihood.cpp
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

#include "DoubleRecursiveTreeLikelihoodData.h"
#include "../PatternTools.h"

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

// From the STL:
#include <memory>

using namespace bpp;
using namespace newlik;
using namespace std;

/******************************************************************************/

void DoubleRecursiveTreeLikelihoodData::initLikelihoods(const SiteContainer& sites, const SubstitutionProcess& process)
throw (Exception)
{
  if (sites.getNumberOfSequences() == 1) throw Exception("DoubleRecursiveTreeLikelihoodData::initLikelihoods. Only 1 sequence in data set.");
  if (sites.getNumberOfSequences() == 0) throw Exception("DoubleRecursiveTreeLikelihoodData::initLikelihoods. No sequence in data set.");
  if (!process.isCompatibleWith(sites))
    throw Exception("DoubleRecursiveTreeLikelihoodData::initLikelihoods. Data and model are not compatible.");
  alphabet_ = sites.getAlphabet();
  nbStates_ = process.getNumberOfStates();
  nbSites_  = sites.getNumberOfSites();
  auto_ptr<SitePatterns> patterns;

  patterns.reset(new SitePatterns(&sites));
  shrunkData_.reset(patterns->getSites());
  rootWeights_      = patterns->getWeights();
  rootPatternLinks_ = patterns->getIndices();
  nbDistinctSites_  = shrunkData_->getNumberOfSites();

  
  // ??? A remettre???
  
  // Init data:
  // Clone data for more efficiency on sequences access:
  // const SiteContainer* sequences = new AlignedSequenceContainer(*shrunkData_);
  // initLikelihoods_(tree_->getRootNode(), *sequences, process);
  // delete sequences;
  
  initLikelihoods_(process.getTree().getRootNode(), *shrunkData_, process);

  // Now initialize root likelihoods and derivatives:
  rootLikelihoods_.resize(nbClasses_);
  rootLikelihoodsS_.resize(nbClasses_);
  rootLikelihoodsSC_.resize(nbDistinctSites_);
  for (size_t c = 0; c < nbClasses_; c++)
  {
    VVdouble* rootLikelihoods_c_ = &rootLikelihoods_[c];
    Vdouble* rootLikelihoodsS_c_ = &rootLikelihoodsS_[c];
    rootLikelihoods_c_->resize(nbDistinctSites_);
    rootLikelihoodsS_c_->resize(nbDistinctSites_);
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      Vdouble* rootLikelihoods_c_i_ = &(*rootLikelihoods_c_)[i];
      rootLikelihoods_c_i_->resize(nbStates_);
      for (size_t x = 0; x < nbStates_; x++)
      {
        (*rootLikelihoods_c_i_)[x] = 1.;
      }
    }
  }
}


/******************************************************************************/

void DoubleRecursiveTreeLikelihoodData::initLikelihoods_(const Node* node, const SiteContainer& sequences, const SubstitutionProcess& process) throw (Exception)
{
  if (node->isLeaf())
  {
    // Init leaves likelihoods:
    const Sequence* seq;
    try
    {
      seq = &sequences.getSequence(node->getName());
    }
    catch (SequenceNotFoundException& snfe)
    {
      throw SequenceNotFoundException("DoubleRecursiveTreeLikelihoodData::initlikelihoods_. Leaf name in tree not found in site container:", (node->getName()));
    }
    DoubleRecursiveTreeLikelihoodLeafData* leafData = &leafData_[node->getId()];
    VVdouble* leavesLikelihoods_leaf = &leafData->getLikelihoodArray();
    leafData->setNodeId(node->getId());
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
        ( *leavesLikelihoods_leaf_i)[s] = process.getInitValue(s, state);
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
    initLikelihoods_(node->getSon(l), sequences, process);
  }

  // Initialize likelihood vector:
  DoubleRecursiveTreeLikelihoodNodeData* nodeData = &nodeData_[node->getId()];

  std::map<int, VVVdouble>* likelihoods_node = &nodeData->getLikelihoodArrays();
  nodeData->setNodeId(node->getId());

  int nbSons = static_cast<int>(node->getNumberOfSons());

  for (int n = (node->hasFather() ? -1 : 0); n < nbSons; n++)
  {
    const Node* neighbor = (*node)[n];
    VVVdouble* likelihoods_node_neighbor = &(*likelihoods_node)[neighbor->getId()];
    
    likelihoods_node_neighbor->resize(nbClasses_);
    for (size_t c = 0; c < nbClasses_; c++)
    {
      VVdouble* likelihoods_node_neighbor_c = &(*likelihoods_node_neighbor)[c];

      likelihoods_node_neighbor_c->resize(nbDistinctSites_);

      if (neighbor->isLeaf())
      {
        VVdouble* leavesLikelihoods_leaf = &leafData_[neighbor->getId()].getLikelihoodArray();

        for (size_t i = 0; i < nbDistinctSites_; i++)
        {
          Vdouble* leavesLikelihoods_leaf_i = &(*leavesLikelihoods_leaf)[i];
          Vdouble* likelihoods_node_neighbor_c_i = &(*likelihoods_node_neighbor_c)[i];
          likelihoods_node_neighbor_c_i->resize(nbStates_);
          for (size_t s = 0; s < nbStates_; s++)
          {
            (*likelihoods_node_neighbor_c_i)[s] = (*leavesLikelihoods_leaf_i)[s];
          }
        }
      }
      else
      {
        for (size_t i = 0; i < nbDistinctSites_; i++)
        {
          Vdouble* likelihoods_node_neighbor_c_i = &(*likelihoods_node_neighbor_c)[i];
          likelihoods_node_neighbor_c_i->resize(nbStates_);
          for (size_t s = 0; s < nbStates_; s++)
          {
            (*likelihoods_node_neighbor_c_i)[s] = 1.; // All likelihoods are initialized to 1.
          }
        }
      }
    }
  }

  // Initialize d and d2 likelihoods:
  VVdouble* dLikelihoods_node = &nodeData->getDLikelihoodArray();
  VVdouble* d2Likelihoods_node = &nodeData->getD2LikelihoodArray();
  dLikelihoods_node->resize(nbClasses_);
  d2Likelihoods_node->resize(nbClasses_);
  for (size_t c = 0; c < nbClasses_; c++)
  {
    (*dLikelihoods_node)[c].resize(nbDistinctSites_); 
    (*d2Likelihoods_node)[c].resize(nbDistinctSites_);
  }
  
}

/******************************************************************************/

void DoubleRecursiveTreeLikelihoodData::reInit() throw (Exception)
{
  reInit(tree_->getRootNode());
}

void DoubleRecursiveTreeLikelihoodData::reInit(const Node* node) throw (Exception)
{
  if (node->isLeaf())
  {
    DoubleRecursiveTreeLikelihoodLeafData* leafData = &leafData_[node->getId()];
    leafData->setNodeId(node->getId());
  }

  DoubleRecursiveTreeLikelihoodNodeData* nodeData = &nodeData_[node->getId()];
  nodeData->setNodeId(node->getId());
  nodeData->eraseNeighborArrays();

  int nbSons = static_cast<int>(node->getNumberOfSons());

  for (int n = (node->hasFather() ? -1 : 0); n < nbSons; n++)
  {
    const Node* neighbor = (*node)[n];
    VVVdouble* array = &nodeData->getLikelihoodArrayForNeighbor(neighbor->getId());
    for (size_t c = 0; c < nbClasses_; c++){
      VVdouble* array_c = &(*array)[c];

      array_c->resize(nbDistinctSites_);
      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        Vdouble* array_c_i = &(*array_c)[i];
        array_c_i->resize(nbStates_);
        for (size_t s = 0; s < nbStates_; s++)
        {
          (*array_c_i)[s] = 1.; // All likelihoods are initialized to 1.
        }
      }
    }
  }


  VVdouble* dArray=&nodeData->getDLikelihoodArray();
  VVdouble* d2Array=&nodeData->getD2LikelihoodArray();
  
  dArray->resize(nbClasses_);
  d2Array->resize(nbClasses_);
  for (size_t c=0;c<nbClasses_;c++){
    (*dArray)[c].resize(nbDistinctSites_);
    (*d2Array)[c].resize(nbDistinctSites_);
  }

  // We re-initialize each son node:
  size_t nbSonNodes = node->getNumberOfSons();
  for (size_t l = 0; l < nbSonNodes; l++)
  {
    // For each son node,
    reInit(node->getSon(l));
  }  
}

