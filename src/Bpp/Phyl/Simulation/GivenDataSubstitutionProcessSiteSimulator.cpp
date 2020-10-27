//
// File: GivenDataSubstitutionProcessSiteSimulator.cpp
// Created by: Laurent Guéguen
// Created on: mardi 13 octobre 2020, à 22h 15
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

#include "GivenDataSubstitutionProcessSiteSimulator.h"
#include "../NewLikelihood/DataFlow/ForwardLikelihoodTree.h"
#include <algorithm>

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From SeqLib:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

void GivenDataSubstitutionProcessSiteSimulator::init()
{
  // Initialize sons & fathers of tree_ Nodes    
  // set sequence names

  if (outputInternalSites_) {
    auto vCN= phyloTree_->getAllNodes();
    seqNames_.resize(vCN.size());    
    seqIndexes_.resize(vCN.size());    
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      auto index = phyloTree_->getNodeIndex(vCN[i]);
      seqNames_[i] = (phyloTree_->isLeaf(vCN[i]))?vCN[i]->getName():TextTools::toString(index);
      seqIndexes_[i] = index; 
    }
  }
  else {
    auto vCN= phyloTree_->getAllLeaves();
    seqNames_.resize(vCN.size());    
    seqIndexes_.resize(vCN.size());    
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      seqNames_[i] = vCN[i]->getName();
      seqIndexes_[i] = phyloTree_->getNodeIndex(vCN[i]);
    }
  }

  // Set up cumsum rates

  const auto dRate = process_->getRateDistribution();

  qRates_.clear();
  
  for (size_t i=0; i<nbClasses_; i++)
    qRates_.push_back(calcul_->getSiteLikelihoodsForAClass(i)(pos_, true));

  qRates_ /= VectorTools::sum(qRates_);

  // Initialize root frequencies

  const  auto dRoot = process_->getRootFrequencies();
  qRoots_.resize(nbClasses_);
  Vdouble temp(nbStates_);

  for (size_t c=0;c<nbClasses_;c++)
  {
    const auto& siteForwLik = calcul_->getForwardLikelihoodsAtNodeForClass(tree_.getNodeIndex(tree_.getRoot()), c)->getTargetValue().col(pos_);

    for (size_t x = 0; x < nbStates_; x++)
      temp[x] = siteForwLik(x);
    
    temp /= VectorTools::sum(temp);
    
    qRoots_[c] = VectorTools::cumSum(temp);
  }
  
  // Initialize cumulative pxy for edges that have models

  auto edges = tree_.getAllEdges();

  Vdouble postTrans(nbStates_);

  for (auto& edge : edges)
  {
    const auto model = edge->getModel();
    auto outid = tree_.getEdgeIndex(edge);

    if (edge->useProb())
      continue;

    const auto transmodel = dynamic_cast<const TransitionModel*>(model);
    if (!transmodel)
      throw Exception("SubstitutionProcessSiteSimulator::init : model "  + model->getName() + " on branch " + TextTools::toString(tree_.getEdgeIndex(edge)) + " is not a TransitionModel.");
    
    VVVdouble* cumpxy_node_ = &edge->cumpxy_;
    cumpxy_node_->resize(nbClasses_);
    
    for (size_t c = 0; c < nbClasses_; c++)
    {
      double brlen = dRate->getCategory(c) * phyloTree_->getEdge(edge->getSpeciesIndex())->getLength();

      VVdouble* cumpxy_node_c_ = &(*cumpxy_node_)[c];
    
      cumpxy_node_c_->resize(nbStates_);
    
      // process transition probabilities already consider rates &
      // branch length

      const Matrix<double>* P;

      const auto& vSub(edge->subModelNumbers());

      if (vSub.size()==0)
        P = &transmodel->getPij_t(brlen);
      else
      {
        if (vSub.size()>1)
          throw Exception("SubstitutionProcessSiteSimulator::init : only 1 submodel can be used.");
        
        const auto* mmodel = dynamic_cast<const MixedTransitionModel*>(transmodel);
        
        const auto* model2 = mmodel->getNModel(vSub[0]);
        
        P = &model2->getPij_t(brlen);
      }

      /* Get likelihoods on this node for all states at this position*/

      const auto& siteForwLik = calcul_->getForwardLikelihoodsAtNodeForClass(tree_.getSon(outid), c)->getTargetValue().col(pos_);

      for (size_t x = 0; x < nbStates_; x++)
      {
        for (size_t y = 0; y < nbStates_; y++)
          postTrans[y] = std::max((*P)(x, y), NumConstants::PICO()) * siteForwLik(y); // to avoid null trans prob on short branches, and then null sum of the postTrans

        postTrans /= VectorTools::sum(postTrans);
        
        (*cumpxy_node_c_)[x] = VectorTools::cumSum(postTrans);
      }
    }
  }
  
  // Initialize cumulative prob for mixture nodes
  auto nodes = tree_.getAllNodes();

  for (auto node:nodes)
  {
    if (node->isMixture()) // set probas to chose
    {
      auto outEdges = tree_.getOutgoingEdges(node);
      VVdouble vprob;
      vprob.resize(nbClasses_);

      for (auto edge : outEdges)
      {
        auto outid = tree_.getEdgeIndex(edge);
        
        auto model = dynamic_cast<const MixedTransitionModel*>(edge->getModel());
        if (!model)
          throw Exception("SubstitutionProcessSiteSimulator::init : model in edge " + TextTools::toString(tree_.getEdgeIndex(edge)) + " is not a mixture.");

        // a priori probabilities of the edges
        const auto& vNb(edge->subModelNumbers());

        double x=0.;
        for (auto nb:vNb)
          x += model->getNProbability(nb);

        // forward lik
        for (size_t c = 0; c < nbClasses_; c++)
        {
          double forwLik = calcul_->getForwardLikelihoodTree(c)->getEdge(outid)->getTargetValue().col(pos_).sum();
          vprob[c].push_back(x * forwLik);
        }

        node->sons_.push_back(tree_.getSon(edge));
      }

      for (size_t c = 0; c < nbClasses_; c++)
      {
        vprob[c] /= VectorTools::sum(vprob[c]);
        node->cumProb_[c] = VectorTools::cumSum(vprob[c]);
      }
    }
  }

}

/******************************************************************************/



