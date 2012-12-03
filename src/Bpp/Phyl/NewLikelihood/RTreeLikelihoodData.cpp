//
// File: RTreeLikelihoodData.cpp
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

#include "RTreeLikelihoodData.h"
#include "../PatternTools.h"

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

// From the STL:
#include <memory>

using namespace bpp;
using namespace std;

/******************************************************************************/

void RTreeLikelihoodData::initLikelihoods(const SiteContainer& sites, const SubstitutionProcess& process)
throw (Exception)
{
  if (sites.getNumberOfSequences() == 1) throw Exception("Error, only 1 sequence!");
  if (sites.getNumberOfSequences() == 0) throw Exception("Error, no sequence!");
  if (!process.isCompatibleWith(sites))
    throw Exception("RTreeLikelihoodData::initLikelihoods. Data and model are not compatible.");
  alphabet_ = sites.getAlphabet();
  nbStates_ = process.getNumberOfStates();
  nbSites_  = sites.getNumberOfSites();
  if (shrunkData_) delete shrunkData_;
  auto_ptr<SitePatterns> patterns;
  if (usePatterns_)
  {
    patterns.reset(initLikelihoodsWithPatterns_(tree_->getRootNode(), sites, process));
    shrunkData_       = patterns->getSites();
    rootWeights_      = patterns->getWeights();
    rootPatternLinks_ = patterns->getIndices();
    nbDistinctSites_  = shrunkData_->getNumberOfSites();
  }
  else
  {
    patterns.reset(new SitePatterns(&sites));
    shrunkData_       = patterns->getSites();
    rootWeights_      = patterns->getWeights();
    rootPatternLinks_ = patterns->getIndices();
    nbDistinctSites_  = shrunkData_->getNumberOfSites();
    initLikelihoodsWithoutPatterns_(tree_->getRootNode(), *shrunkData_, process);
  }
}

/******************************************************************************/

void RTreeLikelihoodData::initLikelihoodsWithoutPatterns_(const Node* node, const SiteContainer& sequences, const SubstitutionProcess& process) throw (Exception)
{
  // Initialize likelihood vector:
  RTreeLikelihoodNodeData* nodeData = &nodeData_[node->getId()];
  nodeData->setNode(node);
  VVVdouble* likelihoods_node = &nodeData->getLikelihoodArray();
  VVVdouble* dLikelihoods_node = &nodeData->getDLikelihoodArray();
  VVVdouble* d2Likelihoods_node = &nodeData->getD2LikelihoodArray();

  likelihoods_node->resize(nbDistinctSites_);
  dLikelihoods_node->resize(nbDistinctSites_);
  d2Likelihoods_node->resize(nbDistinctSites_);

  for (unsigned int i = 0; i < nbDistinctSites_; i++)
  {
    VVdouble* likelihoods_node_i = &(*likelihoods_node)[i];
    VVdouble* dLikelihoods_node_i = &(*dLikelihoods_node)[i];
    VVdouble* d2Likelihoods_node_i = &(*d2Likelihoods_node)[i];
    likelihoods_node_i->resize(nbClasses_);
    dLikelihoods_node_i->resize(nbClasses_);
    d2Likelihoods_node_i->resize(nbClasses_);
    for (unsigned int c = 0; c < nbClasses_; c++)
    {
      Vdouble* likelihoods_node_i_c = &(*likelihoods_node_i)[c];
      Vdouble* dLikelihoods_node_i_c = &(*dLikelihoods_node_i)[c];
      Vdouble* d2Likelihoods_node_i_c = &(*d2Likelihoods_node_i)[c];
      likelihoods_node_i_c->resize(nbStates_);
      dLikelihoods_node_i_c->resize(nbStates_);
      d2Likelihoods_node_i_c->resize(nbStates_);
      for (unsigned int s = 0; s < nbStates_; s++)
      {
        (*likelihoods_node_i_c)[s] = 1; // All likelihoods are initialized to 1.
        (*dLikelihoods_node_i_c)[s] = 0; // All dLikelihoods are initialized to 0.
        (*d2Likelihoods_node_i_c)[s] = 0; // All d2Likelihoods are initialized to 0.
      }
    }
  }

  // Now initialize likelihood values and pointers:

  if (node->isLeaf())
  {
    const Sequence* seq;
    try
    {
      seq = &sequences.getSequence(node->getName());
    }
    catch (SequenceNotFoundException snfe)
    {
      throw SequenceNotFoundException("RTreeLikelihoodData::initTreelikelihoods. Leaf name in tree not found in site conainer: ", (node->getName()));
    }
    for (unsigned int i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble* likelihoods_node_i = &(*likelihoods_node)[i];
      int state = seq->getValue(i);
      for (unsigned int c = 0; c < nbClasses_; c++)
      {
        Vdouble* likelihoods_node_i_c = &(*likelihoods_node_i)[c];
        double test = 0.;
        for (unsigned int s = 0; s < nbStates_; s++)
        {
          // Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
          // otherwise value set to 0:
          // cout << "i=" << i << "\tc=" << c << "\ts=" << s << endl;
          (*likelihoods_node_i_c)[s] = process.getInitValue(s, state);
          test += (*likelihoods_node_i_c)[s];
        }
        if (test < 0.000001) std::cerr << "WARNING!!! Likelihood will be 0 for this site." << std::endl;
      }
    }
  }
  else
  {
    // 'node' is an internal node.
    std::map<int, std::vector<unsigned int> >* patternLinks_node = &patternLinks_[node->getId()];
    unsigned int nbSonNodes = node->getNumberOfSons();
    for (unsigned int l = 0; l < nbSonNodes; ++l)
    {
      // For each son node,
      const Node* son = (*node)[l];
      initLikelihoodsWithoutPatterns_(son, sequences, process);
      std::vector<unsigned int>* patternLinks_node_son = &(*patternLinks_node)[son->getId()];

      // Init map:
      patternLinks_node_son->resize(nbDistinctSites_);

      for (unsigned int i = 0; i < nbDistinctSites_; ++i)
      {
        (*patternLinks_node_son)[i] = i;
      }
    }
  }
}

/******************************************************************************/

SitePatterns* RTreeLikelihoodData::initLikelihoodsWithPatterns_(const Node* node, const SiteContainer& sequences, const SubstitutionProcess& process) throw (Exception)
{
  auto_ptr<SiteContainer> tmp(PatternTools::getSequenceSubset(sequences, *node));
  auto_ptr<SitePatterns> patterns(new SitePatterns(tmp.get(), false));
  auto_ptr<SiteContainer> subSequences(patterns->getSites());

  unsigned int nbSites = subSequences->getNumberOfSites();

  // Initialize likelihood vector:
  RTreeLikelihoodNodeData* nodeData = &nodeData_[node->getId()];
  nodeData->setNode(node);
  VVVdouble* likelihoods_node = &nodeData->getLikelihoodArray();
  VVVdouble* dLikelihoods_node = &nodeData->getDLikelihoodArray();
  VVVdouble* d2Likelihoods_node = &nodeData->getD2LikelihoodArray();
  likelihoods_node->resize(nbSites);
  dLikelihoods_node->resize(nbSites);
  d2Likelihoods_node->resize(nbSites);

  for (unsigned int i = 0; i < nbSites; ++i)
  {
    VVdouble* likelihoods_node_i = &(*likelihoods_node)[i];
    VVdouble* dLikelihoods_node_i = &(*dLikelihoods_node)[i];
    VVdouble* d2Likelihoods_node_i = &(*d2Likelihoods_node)[i];
    likelihoods_node_i->resize(nbClasses_);
    dLikelihoods_node_i->resize(nbClasses_);
    d2Likelihoods_node_i->resize(nbClasses_);
    for (unsigned int c = 0; c < nbClasses_; ++c)
    {
      Vdouble* likelihoods_node_i_c = &(*likelihoods_node_i)[c];
      Vdouble* dLikelihoods_node_i_c = &(*dLikelihoods_node_i)[c];
      Vdouble* d2Likelihoods_node_i_c = &(*d2Likelihoods_node_i)[c];
      likelihoods_node_i_c->resize(nbStates_);
      dLikelihoods_node_i_c->resize(nbStates_);
      d2Likelihoods_node_i_c->resize(nbStates_);
      for (unsigned int s = 0; s < nbStates_; ++s)
      {
        (*likelihoods_node_i_c)[s] = 1; // All likelihoods are initialized to 1.
        (*dLikelihoods_node_i_c)[s] = 0; // All dLikelihoods are initialized to 0.
        (*d2Likelihoods_node_i_c)[s] = 0; // All d2Likelihoods are initialized to 0.
      }
    }
  }

  // Now initialize likelihood values and pointers:

  if (node->isLeaf())
  {
    const Sequence* seq;
    try
    {
      seq = &subSequences->getSequence(node->getName());
    }
    catch (SequenceNotFoundException snfe)
    {
      throw SequenceNotFoundException("RTreeLikelihoodData::initTreelikelihoodsWithPatterns_. Leaf name in tree not found in site conainer: ", (node->getName()));
    }
    for (unsigned int i = 0; i < nbSites; ++i)
    {
      VVdouble* likelihoods_node_i = &(*likelihoods_node)[i];
      int state = seq->getValue(i);
      for (unsigned int c = 0; c < nbClasses_; ++c)
      {
        Vdouble* likelihoods_node_i_c = &(*likelihoods_node_i)[c];
        double test = 0.;
        for (unsigned int s = 0; s < nbStates_; ++s)
        {
          // Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
          // otherwise value set to 0:
          // cout << "i=" << i << "\tc=" << c << "\ts=" << s << endl;
          (*likelihoods_node_i_c)[s] = process.getInitValue(s, state);
          test += (*likelihoods_node_i_c)[s];
        }
        if (test < 0.000001) std::cerr << "WARNING!!! Likelihood will be 0 for this site." << std::endl;
      }
    }
  }
  else
  {
    // 'node' is an internal node.
    std::map<int, std::vector<unsigned int> >* patternLinks_node = &patternLinks_[node->getId()];

    // Now initialize pattern links:
    unsigned int nbSonNodes = node->getNumberOfSons();
    for (unsigned int l = 0; l < nbSonNodes; ++l)
    {
      // For each son node,
      const Node* son = (*node)[l];

      std::vector<unsigned int>* patternLinks_node_son = &(*patternLinks_node)[son->getId()];

      // Initialize subtree 'l' and retrieves corresponding subSequences:
      auto_ptr<SitePatterns> subPatterns(initLikelihoodsWithPatterns_(son, *subSequences.get(), process));
      (*patternLinks_node_son) = subPatterns->getIndices();
    }
  }
  return patterns.release();
}

/******************************************************************************/

