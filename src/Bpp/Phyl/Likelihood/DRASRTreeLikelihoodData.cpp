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
#include "../PatternTools.h"

// From SeqLib:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

void DRASRTreeLikelihoodData::initLikelihoods(const SiteContainer& sites, const SubstitutionModel& model)
throw (Exception)
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
  if (shrunkData_)
    delete shrunkData_;
  SitePatterns* patterns;
  if (usePatterns_)
  {
    patterns          = initLikelihoodsWithPatterns(tree_->getRootNode(), sites, model);
    shrunkData_       = patterns->getSites();
    rootWeights_      = patterns->getWeights();
    rootPatternLinks_ = patterns->getIndices();
    nbDistinctSites_  = shrunkData_->getNumberOfSites();
  }
  else
  {
    patterns          = new SitePatterns(&sites);
    shrunkData_       = patterns->getSites();
    rootWeights_      = patterns->getWeights();
    rootPatternLinks_ = patterns->getIndices();
    nbDistinctSites_  = shrunkData_->getNumberOfSites();
    initLikelihoods(tree_->getRootNode(), *shrunkData_, model);
  }
  delete patterns;
}

/******************************************************************************/

void DRASRTreeLikelihoodData::initLikelihoods(const Node* node, const SiteContainer& sequences, const SubstitutionModel& model) throw (Exception)
{
  // Initialize likelihood vector:
  DRASRTreeLikelihoodNodeData* nodeData = &nodeData_[node->getId()];
  nodeData->setNode(node);
  VVVdouble* _likelihoods_node = &nodeData->getLikelihoodArray();
  VVVdouble* _dLikelihoods_node = &nodeData->getDLikelihoodArray();
  VVVdouble* _d2Likelihoods_node = &nodeData->getD2LikelihoodArray();

  _likelihoods_node->resize(nbDistinctSites_);
  _dLikelihoods_node->resize(nbDistinctSites_);
  _d2Likelihoods_node->resize(nbDistinctSites_);

  for (size_t i = 0; i < nbDistinctSites_; i++)
  {
    VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
    VVdouble* _dLikelihoods_node_i = &(*_dLikelihoods_node)[i];
    VVdouble* _d2Likelihoods_node_i = &(*_d2Likelihoods_node)[i];
    _likelihoods_node_i->resize(nbClasses_);
    _dLikelihoods_node_i->resize(nbClasses_);
    _d2Likelihoods_node_i->resize(nbClasses_);
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
      Vdouble* _dLikelihoods_node_i_c = &(*_dLikelihoods_node_i)[c];
      Vdouble* _d2Likelihoods_node_i_c = &(*_d2Likelihoods_node_i)[c];
      _likelihoods_node_i_c->resize(nbStates_);
      _dLikelihoods_node_i_c->resize(nbStates_);
      _d2Likelihoods_node_i_c->resize(nbStates_);
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_likelihoods_node_i_c)[s] = 1; // All likelihoods are initialized to 1.
        (*_dLikelihoods_node_i_c)[s] = 0; // All dLikelihoods are initialized to 0.
        (*_d2Likelihoods_node_i_c)[s] = 0; // All d2Likelihoods are initialized to 0.
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
      throw SequenceNotFoundException("DRASRTreeLikelihoodData::initTreelikelihoods. Leaf name in tree not found in site conainer: ", (node->getName()));
    }
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
      int state = seq->getValue(i);
      for (size_t c = 0; c < nbClasses_; c++)
      {
        Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
        double test = 0.;
        for (size_t s = 0; s < nbStates_; s++)
        {
          // Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
          // otherwise value set to 0:
          (*_likelihoods_node_i_c)[s] = model.getInitValue(s, state);
          test += (*_likelihoods_node_i_c)[s];
        }
        if (test < 0.000001)
          std::cerr << "WARNING!!! Likelihood will be 0 for this site." << std::endl;
      }
    }
  }
  else
  {
    // 'node' is an internal node.
    std::map<int, std::vector<size_t> >* patternLinks__node = &patternLinks_[node->getId()];
    size_t nbSonNodes = node->getNumberOfSons();
    for (size_t l = 0; l < nbSonNodes; l++)
    {
      // For each son node,
      const Node* son = (*node)[static_cast<int>(l)];
      initLikelihoods(son, sequences, model);
      std::vector<size_t>* patternLinks__node_son = &(*patternLinks__node)[son->getId()];

      // Init map:
      patternLinks__node_son->resize(nbDistinctSites_);

      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        (*patternLinks__node_son)[i] = i;
      }
    }
  }
}

/******************************************************************************/

SitePatterns* DRASRTreeLikelihoodData::initLikelihoodsWithPatterns(const Node* node, const SiteContainer& sequences, const SubstitutionModel& model) throw (Exception)
{
  SiteContainer* tmp = PatternTools::getSequenceSubset(sequences, *node);
  SitePatterns* patterns = new SitePatterns(tmp, true);
  SiteContainer* subSequences = patterns->getSites();

  size_t nbSites = subSequences->getNumberOfSites();

  // Initialize likelihood vector:
  DRASRTreeLikelihoodNodeData* nodeData = &nodeData_[node->getId()];
  nodeData->setNode(node);
  VVVdouble* _likelihoods_node = &nodeData->getLikelihoodArray();
  VVVdouble* _dLikelihoods_node = &nodeData->getDLikelihoodArray();
  VVVdouble* _d2Likelihoods_node = &nodeData->getD2LikelihoodArray();
  _likelihoods_node->resize(nbSites);
  _dLikelihoods_node->resize(nbSites);
  _d2Likelihoods_node->resize(nbSites);

  for (size_t i = 0; i < nbSites; i++)
  {
    VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
    VVdouble* _dLikelihoods_node_i = &(*_dLikelihoods_node)[i];
    VVdouble* _d2Likelihoods_node_i = &(*_d2Likelihoods_node)[i];
    _likelihoods_node_i->resize(nbClasses_);
    _dLikelihoods_node_i->resize(nbClasses_);
    _d2Likelihoods_node_i->resize(nbClasses_);
    for (size_t c = 0; c < nbClasses_; c++)
    {
      Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
      Vdouble* _dLikelihoods_node_i_c = &(*_dLikelihoods_node_i)[c];
      Vdouble* _d2Likelihoods_node_i_c = &(*_d2Likelihoods_node_i)[c];
      _likelihoods_node_i_c->resize(nbStates_);
      _dLikelihoods_node_i_c->resize(nbStates_);
      _d2Likelihoods_node_i_c->resize(nbStates_);
      for (size_t s = 0; s < nbStates_; s++)
      {
        (*_likelihoods_node_i_c)[s] = 1; // All likelihoods are initialized to 1.
        (*_dLikelihoods_node_i_c)[s] = 0; // All dLikelihoods are initialized to 0.
        (*_d2Likelihoods_node_i_c)[s] = 0; // All d2Likelihoods are initialized to 0.
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
      throw SequenceNotFoundException("HomogeneousTreeLikelihood::initTreelikelihoodsWithPatterns. Leaf name in tree not found in site conainer: ", (node->getName()));
    }
    for (size_t i = 0; i < nbSites; i++)
    {
      VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
      int state = seq->getValue(i);
      for (size_t c = 0; c < nbClasses_; c++)
      {
        Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
        double test = 0.;
        for (size_t s = 0; s < nbStates_; s++)
        {
          // Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
          // otherwise value set to 0:
          // cout << "i=" << i << "\tc=" << c << "\ts=" << s << endl;
          (*_likelihoods_node_i_c)[s] = model.getInitValue(s, state);
          test += (*_likelihoods_node_i_c)[s];
        }
        if (test < 0.000001)
          std::cerr << "WARNING!!! Likelihood will be 0 for this site." << std::endl;
      }
    }
  }
  else
  {
    // 'node' is an internal node.
    std::map<int, std::vector<size_t> >* patternLinks__node = &patternLinks_[node->getId()];

    // Now initialize pattern links:
    size_t nbSonNodes = node->getNumberOfSons();
    for (int l = 0; l < static_cast<int>(nbSonNodes); l++)
    {
      // For each son node,
      const Node* son = (*node)[l];

      std::vector<size_t>* patternLinks__node_son = &(*patternLinks__node)[son->getId()];

      // Initialize subtree 'l' and retrieves corresponding subSequences:
      SitePatterns* subPatterns = initLikelihoodsWithPatterns(son, *subSequences, model);
      (*patternLinks__node_son) = subPatterns->getIndices();
      delete subPatterns;
    }
  }
  delete subSequences;
  return patterns;
}

/******************************************************************************/

