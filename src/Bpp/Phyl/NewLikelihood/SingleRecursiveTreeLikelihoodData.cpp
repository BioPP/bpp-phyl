//
// File: SingleRecursiveTreeLikelihoodData.cpp
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

#include "SingleRecursiveTreeLikelihoodData.h"
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

void SingleRecursiveTreeLikelihoodData::initLikelihoods(const SiteContainer& sites, const SubstitutionProcess& process)
throw (Exception)
{
  if (sites.getNumberOfSequences() == 1) throw Exception("SingleRecursiveTreeLikelihoodData::initLikelihoods. Only 1 sequence in data set.");
  if (sites.getNumberOfSequences() == 0) throw Exception("SingleRecursiveTreeLikelihoodData::initLikelihoods. No sequence in data set.");
  if (!process.isCompatibleWith(sites))
    throw Exception("SingleRecursiveTreeLikelihoodData::initLikelihoods. Data and model are not compatible.");
  alphabet_ = sites.getAlphabet();
  nbStates_ = process.getNumberOfStates();
  nbSites_  = sites.getNumberOfSites();
  auto_ptr<SitePatterns> patterns;
  if (usePatterns_)
  {
    patterns.reset(initLikelihoodsWithPatterns_(process.getTree().getRootNode(), sites, process));
    shrunkData_.reset(patterns->getSites());
    rootWeights_      = patterns->getWeights();
    rootPatternLinks_ = patterns->getIndices();
    nbDistinctSites_  = shrunkData_->getNumberOfSites();
  }
  else
  {
    patterns.reset(new SitePatterns(&sites));
    shrunkData_.reset(patterns->getSites());
    rootWeights_      = patterns->getWeights();
    rootPatternLinks_ = patterns->getIndices();
    nbDistinctSites_  = shrunkData_->getNumberOfSites();
    initLikelihoodsWithoutPatterns_(process.getTree().getRootNode(), *shrunkData_, process);
  }
}

/******************************************************************************/

void SingleRecursiveTreeLikelihoodData::initLikelihoodsWithoutPatterns_(const Node* node, const SiteContainer& sequences, const SubstitutionProcess& process) throw (Exception)
{
  // Initialize likelihood vector:
  SingleRecursiveTreeLikelihoodNodeData* nodeData = &nodeData_[node->getId()];
  nodeData->setNodeId(node->getId());
  VVVdouble* likelihoods_node = &nodeData->getLikelihoodArray();
  VVVdouble* dLikelihoods_node = &nodeData->getDLikelihoodArray();
  VVVdouble* d2Likelihoods_node = &nodeData->getD2LikelihoodArray();

  likelihoods_node->resize(nbClasses_);
  dLikelihoods_node->resize(nbClasses_);
  d2Likelihoods_node->resize(nbClasses_);

  for (size_t c = 0; c < nbClasses_; c++)
  {
    VVdouble* likelihoods_node_c = &(*likelihoods_node)[c];
    VVdouble* dLikelihoods_node_c = &(*dLikelihoods_node)[c];
    VVdouble* d2Likelihoods_node_c = &(*d2Likelihoods_node)[c];
    likelihoods_node_c->resize(nbDistinctSites_);
    dLikelihoods_node_c->resize(nbDistinctSites_);
    d2Likelihoods_node_c->resize(nbDistinctSites_);
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      Vdouble* likelihoods_node_c_i = &(*likelihoods_node_c)[i];
      Vdouble* dLikelihoods_node_c_i = &(*dLikelihoods_node_c)[i];
      Vdouble* d2Likelihoods_node_c_i = &(*d2Likelihoods_node_c)[i];
      likelihoods_node_c_i->resize(nbStates_);
      dLikelihoods_node_c_i->resize(nbStates_);
      d2Likelihoods_node_c_i->resize(nbStates_);
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*likelihoods_node_c_i)[s] = 1; // All likelihoods are initialized to 1.
        (*dLikelihoods_node_c_i)[s] = 0; // All dLikelihoods are initialized to 0.
        (*d2Likelihoods_node_c_i)[s] = 0; // All d2Likelihoods are initialized to 0.
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
      throw SequenceNotFoundException("SingleRecursiveTreeLikelihoodData::initTreelikelihoods. Leaf name in tree not found in site conainer: ", (node->getName()));
    }
    for (size_t c = 0; c < nbClasses_; c++)
    {
      VVdouble* likelihoods_node_c = &(*likelihoods_node)[c];
      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        int state = seq->getValue(i);
        Vdouble* likelihoods_node_c_i = &(*likelihoods_node_c)[i];
        double test = 0.;
        for (size_t s = 0; s < nbStates_; s++)
        {
          // Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
          // otherwise value set to 0:
          // cout << "i=" << i << "\tc=" << c << "\ts=" << s << endl;
          (*likelihoods_node_c_i)[s] = process.getInitValue(s, state);
          test += (*likelihoods_node_c_i)[s];
        }
        if (test < 0.000001) std::cerr << "WARNING!!! Likelihood will be 0 for this site." << std::endl;
      }
    }
  }
  else
  {
    // 'node' is an internal node.
    std::map<int, std::vector<size_t> >* patternLinks_node = &patternLinks_[node->getId()];
    int nbSonNodes = static_cast<int>(node->getNumberOfSons());
    for (int l = 0; l < nbSonNodes; ++l)
    {
      // For each son node,
      const Node* son = (*node)[l];
      initLikelihoodsWithoutPatterns_(son, sequences, process);
      std::vector<size_t>* patternLinks_node_son = &(*patternLinks_node)[son->getId()];

      // Init map:
      patternLinks_node_son->resize(nbDistinctSites_);

      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        (*patternLinks_node_son)[i] = i;
      }
    }
  }
}

/******************************************************************************/

SitePatterns* SingleRecursiveTreeLikelihoodData::initLikelihoodsWithPatterns_(const Node* node, const SiteContainer& sequences, const SubstitutionProcess& process) throw (Exception)
{
  SiteContainer* tmp = PatternTools::getSequenceSubset(sequences, *node);
  auto_ptr<SitePatterns> patterns(new SitePatterns(tmp, true)); //Important: patterns own tmp, otherwise sizes will not be accessible outside this function.
  auto_ptr<SiteContainer> subSequences(patterns->getSites());

  size_t nbSites = subSequences->getNumberOfSites();

  // Initialize likelihood vector:
  SingleRecursiveTreeLikelihoodNodeData* nodeData = &nodeData_[node->getId()];
  nodeData->setNodeId(node->getId());
  VVVdouble* likelihoods_node = &nodeData->getLikelihoodArray();
  VVVdouble* dLikelihoods_node = &nodeData->getDLikelihoodArray();
  VVVdouble* d2Likelihoods_node = &nodeData->getD2LikelihoodArray();
  likelihoods_node->resize(nbClasses_);
  dLikelihoods_node->resize(nbClasses_);
  d2Likelihoods_node->resize(nbClasses_);

  for (size_t c = 0; c < nbClasses_; ++c)
  {
    VVdouble* likelihoods_node_c = &(*likelihoods_node)[c];
    VVdouble* dLikelihoods_node_c = &(*dLikelihoods_node)[c];
    VVdouble* d2Likelihoods_node_c = &(*d2Likelihoods_node)[c];
    likelihoods_node_c->resize(nbSites);
    dLikelihoods_node_c->resize(nbSites);
    d2Likelihoods_node_c->resize(nbSites);
    for (size_t i = 0; i < nbSites; ++i)
    {
      Vdouble* likelihoods_node_c_i = &(*likelihoods_node_c)[i];
      Vdouble* dLikelihoods_node_c_i = &(*dLikelihoods_node_c)[i];
      Vdouble* d2Likelihoods_node_c_i = &(*d2Likelihoods_node_c)[i];
      likelihoods_node_c_i->resize(nbStates_);
      dLikelihoods_node_c_i->resize(nbStates_);
      d2Likelihoods_node_c_i->resize(nbStates_);
      for (size_t s = 0; s < nbStates_; ++s)
      {
        (*likelihoods_node_c_i)[s] = 1; // All likelihoods are initialized to 1.
        (*dLikelihoods_node_c_i)[s] = 0; // All dLikelihoods are initialized to 0.
        (*d2Likelihoods_node_c_i)[s] = 0; // All d2Likelihoods are initialized to 0.
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
      throw SequenceNotFoundException("SingleRecursiveTreeLikelihoodData::initTreelikelihoodsWithPatterns_. Leaf name in tree not found in site conainer: ", (node->getName()));
    }
    for (size_t c = 0; c < nbClasses_; ++c)
    {
      VVdouble* likelihoods_node_c = &(*likelihoods_node)[c];
      for (size_t i = 0; i < nbSites; ++i)
      {
        int state = seq->getValue(i);
        Vdouble* likelihoods_node_c_i = &(*likelihoods_node_c)[i];
        double test = 0.;
        for (size_t s = 0; s < nbStates_; ++s)
        {
          // Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
          // otherwise value set to 0:
          // cout << "i=" << i << "\tc=" << c << "\ts=" << s << endl;
          (*likelihoods_node_c_i)[s] = process.getInitValue(s, state);
          test += (*likelihoods_node_c_i)[s];
        }
        if (test < 0.000001) std::cerr << "WARNING!!! Likelihood will be 0 for this site." << std::endl;
      }
    }
  }
  else
  {
    // 'node' is an internal node.
    std::map<int, std::vector<size_t> >* patternLinks_node = &patternLinks_[node->getId()];

    // Now initialize pattern links:
    int nbSonNodes = static_cast<int>(node->getNumberOfSons());
    for (int l = 0; l < nbSonNodes; l++)
    {
      // For each son node,
      const Node* son = (*node)[l];

      std::vector<size_t>* patternLinks_node_son = &(*patternLinks_node)[son->getId()];

      // Initialize subtree 'l' and retrieves corresponding subSequences:
      auto_ptr<SitePatterns> subPatterns(initLikelihoodsWithPatterns_(son, *subSequences.get(), process));
      (*patternLinks_node_son) = subPatterns->getIndices();
    }
  }
  return patterns.release();
}

/******************************************************************************/

