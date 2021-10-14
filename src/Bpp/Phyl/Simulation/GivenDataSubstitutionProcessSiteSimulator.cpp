//
// File: GivenDataSubstitutionProcessSiteSimulator.cpp
// Authors:
//   Laurent GuÃÂ©guen
// Created: mardi 13 octobre 2020, ÃÂ  22h 15
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <algorithm>

#include "../Likelihood/DataFlow/ForwardLikelihoodTree.h"
#include "GivenDataSubstitutionProcessSiteSimulator.h"

// From SeqLib:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;
using namespace std;
using namespace numeric;

/******************************************************************************/

void GivenDataSubstitutionProcessSiteSimulator::init()
{
  calcul_->makeLikelihoods();

  // Initialize sons & fathers of tree_ Nodes
  // set sequence names

  outputInternalSites(outputInternalSites_);

  // Set up cumsum rates

  const auto& dRate = process_->getRateDistribution();

  std::vector<DataLik> qR(nbClasses_);
  for (size_t i = 0; i < nbClasses_; i++)
  {
    qR[i] = calcul_->getSiteLikelihoodsForAClass(i)(pos_, true);
  }

  auto sQ = VectorTools::sum(qR);

  qRates_.clear();
  for (size_t i = 0; i < nbClasses_; i++)
  {
    qRates_.push_back(convert(qR[i] / sQ));
  }

  // Initialize root frequencies

  const auto dRoot = process_->getRootFrequencies();
  qRoots_.resize(nbClasses_);
  RowLik temp((int)nbStates_);

  for (size_t c = 0; c < nbClasses_; c++)
  {
    temp = calcul_->getForwardLikelihoodsAtNodeForClass(tree_.getNodeIndex(tree_.getRoot()), c)->getTargetValue().col(pos_);

    temp /= temp.sum();

    Vdouble temp2;
    copyEigenToBpp(temp, temp2);

    qRoots_[c] = VectorTools::cumSum(temp2);
  }

  // Initialize cumulative pxy for edges that have models

  auto edges = tree_.getAllEdges();

  std::vector<DataLik> postTrans(nbStates_);

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

      if (vSub.size() == 0)
        P = &transmodel->getPij_t(brlen);
      else
      {
        if (vSub.size() > 1)
          throw Exception("SubstitutionProcessSiteSimulator::init : only 1 submodel can be used.");

        const auto* mmodel = dynamic_cast<const MixedTransitionModel*>(transmodel);

        const auto* model2 = mmodel->getNModel(vSub[0]);

        P = &model2->getPij_t(brlen);
      }

      /* Get likelihoods on this node for all states at this position*/

      const auto& siteForwLik = calcul_->getForwardLikelihoodsAtNodeForClass(tree_.getSon(outid), c)->getTargetValue().col(pos_);

      for (size_t x = 0; x < size_t(nbStates_); x++)
      {
        for (size_t y = 0; y < size_t(nbStates_); y++)
        {
          postTrans[y] = std::max((*P)(x, y), NumConstants::PICO()) * siteForwLik(Eigen::Index(y)); // to avoid null trans prob on short branches, and then null sum of the postTrans
        }
        postTrans /= VectorTools::sum(postTrans);

        Vdouble pT2(postTrans.size());
        for (size_t i = 0; i < postTrans.size(); i++)
        {
          pT2[i] = convert(postTrans[i]);
        }

        (*cumpxy_node_c_)[x] = VectorTools::cumSum(pT2);
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
      std::vector<std::vector<DataLik> > vprob;
      vprob.resize(nbClasses_);

      for (auto edge : outEdges)
      {
        auto outid = tree_.getEdgeIndex(edge);

        auto model = dynamic_cast<const MixedTransitionModel*>(edge->getModel());
        if (!model)
          throw Exception("SubstitutionProcessSiteSimulator::init : model in edge " + TextTools::toString(tree_.getEdgeIndex(edge)) + " is not a mixture.");

        // a priori probabilities of the edges
        const auto& vNb(edge->subModelNumbers());

        double x = 0.;
        for (auto nb:vNb)
        {
          x += model->getNProbability(nb);
        }

        // forward lik
        for (size_t c = 0; c < nbClasses_; c++)
        {
          auto forwLik = calcul_->getForwardLikelihoodTree(c)->getEdge(outid)->getTargetValue().col(pos_).sum();
          vprob[c].push_back(x * forwLik);
        }

        node->sons_.push_back(tree_.getSon(edge));
      }

      for (size_t c = 0; c < nbClasses_; c++)
      {
        vprob[c] /= VectorTools::sum(vprob[c]);

        Vdouble pT2(vprob.size());

        // convert to double
        for (size_t i = 0; i < vprob.size(); i++)
        {
          pT2[i] = convert(vprob[c][i]);
        }

        node->cumProb_[c] = VectorTools::cumSum(pT2);
      }
    }
  }
}

/******************************************************************************/
