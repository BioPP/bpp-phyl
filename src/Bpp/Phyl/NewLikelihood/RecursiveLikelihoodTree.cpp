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
#include "ComputingNode.h"
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
/*  RecursiveLikelihoodTree                                         */
/******************************************************************************/

RecursiveLikelihoodTree::RecursiveLikelihoodTree(const SubstitutionProcess& process, bool usepatterns) :
  AbstractLikelihoodTree(process),
  vTree_(),
  patternLinks_(),
  usePatterns_(usepatterns),
  initializedAboveLikelihoods_(false)
{
  TreeTemplate<Node> tree=process.getParametrizableTree().getTree();

  RecursiveLikelihoodNode* rCN= TreeTemplateTools::cloneSubtree<RecursiveLikelihoodNode>(*tree.getRootNode());
  
  TreeTemplate<RecursiveLikelihoodNode>* pTC=new TreeTemplate<RecursiveLikelihoodNode>(rCN);

  for (size_t i=0; i<nbClasses_; i++){
    TreeTemplate<RecursiveLikelihoodNode>* pTC2=pTC->clone();
    vTree_.push_back(pTC2);
  }

  delete pTC;
}

RecursiveLikelihoodTree::RecursiveLikelihoodTree(const RecursiveLikelihoodTree& data):
  AbstractLikelihoodTree(data),
  vTree_(),
  patternLinks_(data.patternLinks_),
  usePatterns_(data.usePatterns_),
  initializedAboveLikelihoods_(data.initializedAboveLikelihoods_)
{
  for (size_t i=0; i<data.vTree_.size(); i++)
    vTree_.push_back(data.vTree_[i]->clone());
}

RecursiveLikelihoodTree& RecursiveLikelihoodTree::operator=(const RecursiveLikelihoodTree & data)
{
  AbstractLikelihoodTree::operator=(data);

  for (size_t i=0; i<data.vTree_.size(); i++)
    vTree_.push_back(data.vTree_[i]->clone());
  
  patternLinks_      = data.patternLinks_;
  usePatterns_       = data.usePatterns_;
  initializedAboveLikelihoods_ = data.initializedAboveLikelihoods_;

  return *this;
}

RecursiveLikelihoodTree::~RecursiveLikelihoodTree()
{
  for (size_t i=0;i<vTree_.size();i++)
    TreeTemplateTools::deleteSubtree(vTree_[i]->getRootNode());

  vTree_.clear();
}

void RecursiveLikelihoodTree::initLikelihoods(const SiteContainer& sites, const SubstitutionProcess& process)
throw (Exception)
{
  if (sites.getNumberOfSequences() == 1) throw Exception("RecursiveLikelihoodTree::initLikelihoods. Only 1 sequence in data set.");
  if (sites.getNumberOfSequences() == 0) throw Exception("RecursiveLikelihoodTree::initLikelihoods. No sequence in data set.");
  if (!process.isCompatibleWith(sites))
    throw Exception("RecursiveLikelihoodTree::initLikelihoods. Data and model are not compatible.");
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

    setPatterns(patternLinks_);
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

void RecursiveLikelihoodTree::initLikelihoodsWithoutPatterns_(const Node* node, const SiteContainer& sequences, const SubstitutionProcess& process) throw (Exception)
{
  int nId=node->getId();

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

  if (node->isLeaf())
  {
    const Sequence* seq;
    try
    {
      seq = &sequences.getSequence(node->getName());
    }
    catch (SequenceNotFoundException snfe)
    {
      throw SequenceNotFoundException("RecursiveLikelihoodTree::initTreelikelihoods. Leaf name in tree not found in site conainer: ", (node->getName()));
    }

    for (size_t c = 0; c < nbClasses_; c++)
    {
      RecursiveLikelihoodNode& lNode= *dynamic_cast<RecursiveLikelihoodNode*>(vTree_[c]->getNode(nId));
      VVdouble& array=lNode.getBelowLikelihoodArray_(ComputingNode::D0);
      
      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        Vdouble* array_i=&array[i];
        int state = seq->getValue(i);
        double test = 0.;
        for (size_t s = 0; s < nbStates_; s++)
        {
          double x = process.getInitValue(s, state);
          if (lNode.usesLog())
          {
            if (x<=0)
              (*array_i)[s]=-10000;
            else
              (*array_i)[s] = log(x);
          }
          else
            (*array_i)[s] = x;
            
          test += x;
        }
        if (test < 0.000001) std::cerr << "WARNING!!! Likelihood will be 0 for this site " << TextTools::toString(i) << std::endl;
      }
      lNode.updateBelow_(true, ComputingNode::D0);
      
    }
  }
  else
  {
    // 'node' is an internal node.
    std::map<int, std::vector<size_t> >* patternLinks_node = &patternLinks_[nId];
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

  if (!node->hasFather())
    setAboveLikelihoods(nId, process.getRootFrequencies());
}

/******************************************************************************/

SitePatterns* RecursiveLikelihoodTree::initLikelihoodsWithPatterns_(const Node* node, const SiteContainer& sequences, const SubstitutionProcess& process) throw (Exception)
{
  SiteContainer* tmp = PatternTools::getSequenceSubset(sequences, *node);
  auto_ptr<SitePatterns> patterns(new SitePatterns(tmp, true)); //Important: patterns own tmp, otherwise sizes will not be accessible outside this function.
  auto_ptr<SiteContainer> subSequences(patterns->getSites());

  size_t nbSites = subSequences->getNumberOfSites();
  int nId=node->getId();

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

  if (node->isLeaf())
  {
    const Sequence* seq;
    try
    {
      seq = &subSequences->getSequence(node->getName());
    }
    catch (SequenceNotFoundException snfe)
    {
      throw SequenceNotFoundException("RecursiveLikelihoodTree::initTreelikelihoodsWithPatterns_. Leaf name in tree not found in site conainer: ", (node->getName()));
    }

    for (size_t c = 0; c < nbClasses_; c++)
    {
      RecursiveLikelihoodNode& lNode= *dynamic_cast<RecursiveLikelihoodNode*>(vTree_[c]->getNode(nId));
      VVdouble& array=lNode.getBelowLikelihoodArray_(ComputingNode::D0);

      for (size_t i = 0; i < nbSites; i++)
      {
        Vdouble* array_i=&array[i];
        
        int state = seq->getValue(i);
        double test = 0.;
        for (size_t s = 0; s < nbStates_; s++)
        {
          double x = process.getInitValue(s, state);
          if (lNode.usesLog())
          {
            if (x<=0)
              (*array_i)[s]=-10000;
            else
              (*array_i)[s] = log(x);
          }
          else
            (*array_i)[s] = x;

          test += x;
        }

        if (test < 0.000001) std::cerr << "WARNING!!! Likelihood will be 0 for this site " << TextTools::toString(i) << std::endl;
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

void RecursiveLikelihoodTree::resetInnerAboveLikelihoods() 
{
  // Initialize likelihood vector:
  for (size_t i=0; i<nbClasses_; i++)
  {
    vector<RecursiveLikelihoodNode*> vNd=vTree_[i]->getNodes();

    for (size_t  j=0; j<vNd.size(); j++)
    {
      if (vNd[j]->getId()!=vTree_[i]->getRootId())
        vNd[j]->resetAboveLikelihoods(nbDistinctSites_, nbStates_);
    }
  }
  
  initializedAboveLikelihoods_ = true;
}
  
