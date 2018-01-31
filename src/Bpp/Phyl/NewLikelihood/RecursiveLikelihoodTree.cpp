//
// File: RecursiveLikelihoodTree.cpp
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

#include "RecursiveLikelihoodTree.h"
#include "SpeciationComputingNode.h"
#include "../PatternTools.h"

// From bpp-seq:

// From the STL:
#include <memory>

using namespace bpp;
using namespace std;


/******************************************************************************/
/*  RecursiveLikelihoodTree                                         */
/******************************************************************************/

RecursiveLikelihoodTree::RecursiveLikelihoodTree(const SubstitutionProcess& process, bool usepatterns) :
  AbstractLikelihoodTree(process),
  vTree_(),
  patternLinks_(),
  usePatterns_(usepatterns),
  initializedAboveLikelihoods_(false)
{
  for (size_t i = 0; i < nbClasses_; i++)
  {
    std::shared_ptr<LikTree> pTC2(new LikTree(process.getParametrizablePhyloTree()));
    std::vector<std::shared_ptr<RecursiveLikelihoodNode> > vCN=pTC2->getAllNodes();
    
    for (size_t j=0; j<vCN.size(); j++)
      vCN[j]->updateTree(pTC2.get(), pTC2->getNodeIndex(vCN[j]));
    
    vTree_.push_back(pTC2);
  }
}

RecursiveLikelihoodTree::RecursiveLikelihoodTree(const RecursiveLikelihoodTree& data) :
  AbstractLikelihoodTree(data),
  vTree_(),
  patternLinks_(data.patternLinks_),
  usePatterns_(data.usePatterns_),
  initializedAboveLikelihoods_(data.initializedAboveLikelihoods_)
{
  for (size_t i = 0; i < data.vTree_.size(); i++)
  {
    std::shared_ptr<LikTree> pTC2(new LikTree(*data.vTree_[i]));
    std::vector<std::shared_ptr<RecursiveLikelihoodNode> > vCN=pTC2->getAllNodes();
    for (size_t j=0; j<vCN.size(); j++)
      vCN[j]->updateTree(pTC2.get(), pTC2->getNodeIndex(vCN[j]));
    vTree_.push_back(pTC2);
  }
}

RecursiveLikelihoodTree& RecursiveLikelihoodTree::operator=(const RecursiveLikelihoodTree& data)
{
  AbstractLikelihoodTree::operator=(data);

  for (size_t i = 0; i < data.vTree_.size(); i++)
  {
    std::shared_ptr<LikTree> pTC2(new LikTree(*data.vTree_[i]));
    std::vector<std::shared_ptr<RecursiveLikelihoodNode> > vCN=pTC2->getAllNodes();
    for (size_t j=0; j<vCN.size(); j++)
      vCN[j]->updateTree(pTC2.get(), pTC2->getNodeIndex(vCN[j]));

    vTree_.push_back(pTC2);
  }

  patternLinks_      = data.patternLinks_;
  usePatterns_       = data.usePatterns_;
  initializedAboveLikelihoods_ = data.initializedAboveLikelihoods_;

  return *this;
}

RecursiveLikelihoodTree::~RecursiveLikelihoodTree()
{
  vTree_.clear();
}

void RecursiveLikelihoodTree::initLikelihoods(const AlignedValuesContainer& sites, const SubstitutionProcess& process)
throw (Exception)
{
  if (sites.getNumberOfSequences() == 1)
    throw Exception("RecursiveLikelihoodTree::initLikelihoods. Only 1 sequence in data set.");
  if (sites.getNumberOfSequences() == 0)
    throw Exception("RecursiveLikelihoodTree::initLikelihoods. No sequence in data set.");
  if (!process.isCompatibleWith(sites))
    throw Exception("RecursiveLikelihoodTree::initLikelihoods. Data and model are not compatible.");
  alphabet_ = sites.getAlphabet();
  nbStates_ = process.getNumberOfStates();
  nbSites_  = sites.getNumberOfSites();
  unique_ptr<SitePatterns> patterns;

  if (usePatterns_)
  {
    patterns.reset(initLikelihoodsWithPatterns_(vTree_[0]->getRoot().get(), sites, process));
    shrunkData_ = patterns->getSites();
    rootWeights_      = patterns->getWeights();
    rootPatternLinks_ = patterns->getIndices();
    nbDistinctSites_  = shrunkData_->getNumberOfSites();

    setPatterns(patternLinks_);
  }
  else
  {
    patterns.reset(new SitePatterns(&sites));
    shrunkData_       = patterns->getSites();
    rootWeights_      = patterns->getWeights();
    rootPatternLinks_ = patterns->getIndices();
    nbDistinctSites_  = shrunkData_->getNumberOfSites();
    initLikelihoodsWithoutPatterns_(vTree_[0]->getRoot().get(), *shrunkData_, process);
  }
}

/******************************************************************************/

void RecursiveLikelihoodTree::initLikelihoodsWithoutPatterns_(const RecursiveLikelihoodNode* node, const AlignedValuesContainer& sequences, const SubstitutionProcess& process) throw (Exception)
{
  const ParametrizablePhyloTree& tree=process.getParametrizablePhyloTree();
  int nId = node->getId();

  // for Model specific Alphabet State
  const StateMap& statemap = process.getStateMap();

  // Initialize likelihood vector:
  if (!node->hasFather())
  {
    resetAboveLikelihoods(nId, nbDistinctSites_, nbStates_);
    resetLikelihoods(nId, nbDistinctSites_, nbStates_, ComputingNode::D0);

    resetLikelihoods(nId, nbDistinctSites_, nbStates_, ComputingNode::D1);
    resetLikelihoods(nId, nbDistinctSites_, nbStates_, ComputingNode::D2);
  }


  resetBelowLikelihoods(nId, nbDistinctSites_, nbStates_, ComputingNode::D0);
  resetBelowLikelihoods(nId, nbDistinctSites_, nbStates_, ComputingNode::D1);
  resetBelowLikelihoods(nId, nbDistinctSites_, nbStates_, ComputingNode::D2);


  // Now initialize likelihood values and pointers:

  
  if (node->hasNoSon())
  {
    size_t posSeq;
    try
    {
      posSeq=sequences.getSequencePosition(tree.getNode(node->getId())->getName());
    }
    catch (SequenceNotFoundException snfe)
    {
      throw SequenceNotFoundException("RecursiveLikelihoodTree::initTreelikelihoods. Leaf name in tree not found in site container: ", tree.getNode(node->getId())->getName());
    }

    for (size_t c = 0; c < nbClasses_; c++)
    {
      RecursiveLikelihoodNode& lNode = dynamic_cast<RecursiveLikelihoodNode&>(*vTree_[c]->getNode(nId));
      VVdouble& array = lNode.getBelowLikelihoodArray_(ComputingNode::D0);

      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        Vdouble* array_i = &array[i];

        double test = 0.;
        for (size_t s = 0; s < nbStates_; s++)
        {
          double x = sequences.getStateValueAt(i, posSeq, statemap.getAlphabetStateAsInt(s));
          
          if (lNode.usesLog())
          {
            if (x <= 0)
              (*array_i)[s] = -10000;
            else
              (*array_i)[s] = log(x);
          }
          else
            (*array_i)[s] = x;

          test += x;
        }

        if (test < 0.000001)
          std::cerr << "WARNING!!! Likelihood will be 0 for this site " << TextTools::toString(i) << std::endl;
      }
      lNode.updateBelow_(true, ComputingNode::D0);
    }
  }
  else
  {
    // 'node' is an internal node.
    std::map<int, std::vector<size_t> >* patternLinks_node = &patternLinks_[nId];
    size_t nbSonNodes = node->getNumberOfSons();
    for (size_t l = 0; l < nbSonNodes; ++l)
    {
      // For each son node,
      const RecursiveLikelihoodNode* son = dynamic_cast<const RecursiveLikelihoodNode*>((*node)[(int)l]);
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

  if (!node->hasFather())
    setAboveLikelihoods(nId, process.getRootFrequencies());
}

/******************************************************************************/

SitePatterns* RecursiveLikelihoodTree::initLikelihoodsWithPatterns_(const RecursiveLikelihoodNode* node, const AlignedValuesContainer& sequences, const SubstitutionProcess& process) throw (Exception)
{
  const ParametrizablePhyloTree& tree=process.getParametrizablePhyloTree();
  AlignedValuesContainer* tmp = PatternTools::getSequenceSubset(sequences, tree.getNode(node->getId()), tree);

  SitePatterns* patterns = new SitePatterns(tmp, true); // Important: patterns own tmp, otherwise sizes will not be accessible outside this function.
  shared_ptr<AlignedValuesContainer> subSequences(patterns->getSites());

  size_t nbSites = subSequences->getNumberOfSites();
  int nId = node->getId();

  // for Model specific Alphabet State
  const StateMap& statemap = process.getStateMap();

  // Initialize likelihood vectors:

  if (!node->hasFather())
  {
    resetAboveLikelihoods(nId, nbSites, nbStates_);
    resetLikelihoods(nId, nbSites, nbStates_, ComputingNode::D0);
    resetLikelihoods(nId, nbSites, nbStates_, ComputingNode::D1);
    resetLikelihoods(nId, nbSites, nbStates_, ComputingNode::D2);
  }

  resetBelowLikelihoods(nId, nbSites, nbStates_, ComputingNode::D0);

  resetBelowLikelihoods(nId, nbSites, nbStates_, ComputingNode::D1);
  resetBelowLikelihoods(nId, nbSites, nbStates_, ComputingNode::D2);


  // Now initialize likelihood values and pointers:

  if (node->hasNoSon())
  {
    size_t posSeq;
    try
    {
      posSeq=subSequences->getSequencePosition(tree.getNode(node->getId())->getName());
    }
    catch (SequenceNotFoundException snfe)
    {
      throw SequenceNotFoundException("RecursiveLikelihoodTree::initTreelikelihoodsWithPatterns_. Leaf name in tree not found in site container: ", tree.getNode(node->getId())->getName());
    }

    for (size_t c = 0; c < nbClasses_; c++)
    {
      RecursiveLikelihoodNode& lNode = dynamic_cast<RecursiveLikelihoodNode&>(*vTree_[c]->getNode(nId));
      VVdouble& array = lNode.getBelowLikelihoodArray_(ComputingNode::D0);

      for (size_t i = 0; i < nbSites; i++)
      {
        Vdouble* array_i = &array[i];

        double test = 0.;
        for (size_t s = 0; s < nbStates_; s++)
        {
          double x = subSequences->getStateValueAt(i, posSeq, statemap.getAlphabetStateAsInt(s));

          if (lNode.usesLog())
          {
            if (x <= 0)
              (*array_i)[s] = -10000;
            else
              (*array_i)[s] = log(x);
          }
          else
            (*array_i)[s] = x;

          test += x;
        }

        if (test < 0.000001)
          std::cerr << "WARNING!!! Likelihood will be 0 for this site " << TextTools::toString(i) << std::endl;
      }
      lNode.updateBelow_(true, ComputingNode::D0);
    }
  }
  else
  {
    // 'node' is an internal node.
    std::map<int, std::vector<size_t> >* patternLinks_node = &patternLinks_[nId];

    // Now initialize pattern links:
    int nbSonNodes = static_cast<int>(node->getNumberOfSons());
    for (int l = 0; l < nbSonNodes; l++)
    {
      // For each son node,
      const RecursiveLikelihoodNode* son = dynamic_cast<const RecursiveLikelihoodNode*>((*node)[l]);

      std::vector<size_t>* patternLinks_node_son = &(*patternLinks_node)[son->getId()];


// Initialize subtree 'l' and retrieves corresponding subSequences:
      unique_ptr<SitePatterns> subPatterns(initLikelihoodsWithPatterns_(son, *subSequences.get(), process));
      (*patternLinks_node_son) = subPatterns->getIndices();
    }
  }
  
  if (!node->hasFather())
    setAboveLikelihoods(nId, process.getRootFrequencies());

  return patterns;
}


/******************************************************************************/

void RecursiveLikelihoodTree::resetInnerAboveLikelihoods()
{
  // Initialize likelihood vector:
  for (size_t i = 0; i < nbClasses_; i++)
  {
    vector<shared_ptr<RecursiveLikelihoodNode> > vNd = vTree_[i]->getAllNodes();

    for (size_t j = 0; j < vNd.size(); j++)
    {
      if (vNd[j] != vTree_[i]->getRoot())
        vNd[j]->resetAboveLikelihoods(nbDistinctSites_, nbStates_);
    }
  }

  initializedAboveLikelihoods_ = true;
}
