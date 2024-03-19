// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/VectorTools.h>

#include "DRTreeLikelihoodTools.h"

using namespace bpp;

// -----------------------------------------------------------------------------------------

VVVdouble DRTreeLikelihoodTools::getPosteriorProbabilitiesPerStatePerRate(
    const DRTreeLikelihoodInterface& drl,
    int nodeId)
{
  size_t nSites   = drl.likelihoodData().getNumberOfDistinctSites();
  size_t nClasses = drl.getNumberOfClasses();
  size_t nStates  = drl.getNumberOfStates();
  VVVdouble postProb(nSites);

  auto rDist = drl.getRateDistribution();
  Vdouble rcProbs = rDist->getProbabilities();
  if (drl.tree().isLeaf(nodeId))
  {
    VVdouble larray = drl.likelihoodData().getLeafLikelihoods(nodeId);
    for (size_t i = 0; i < nSites; ++i)
    {
      VVdouble* postProb_i = &postProb[i];
      postProb_i->resize(nClasses);
      Vdouble* larray_i = &larray[i];
      // In case of generic character:
      double sumprobs = VectorTools::sum(*larray_i);
      for (size_t c = 0; c < nClasses; c++)
      {
        Vdouble* postProb_i_c = &(*postProb_i)[c];
        postProb_i_c->resize(nStates);
        double* rcProb = &rcProbs[c];
        for (size_t x = 0; x < nStates; x++)
        {
          (*postProb_i_c)[x] = (*larray_i)[x] * (*rcProb) / sumprobs;
        }
      }
    }
  }
  else
  {
    VVVdouble larray;
    drl.computeLikelihoodAtNode(nodeId, larray);

    Vdouble likelihoods(nSites, 0);
    for (size_t i = 0; i < nSites; i++)
    {
      VVdouble* larray_i = &larray[i];
      for (size_t c = 0; c < nClasses; c++)
      {
        Vdouble* larray_i_c = &(*larray_i)[c];
        for (size_t s = 0; s < nStates; s++)
        {
          likelihoods[i] += (*larray_i_c)[s];
        }
      }
    }

    for (size_t i = 0; i < nSites; i++)
    {
      VVdouble* postProb_i = &postProb[i];
      postProb_i->resize(nClasses);
      VVdouble* larray_i = &larray[i];
      double likelihood = likelihoods[i];
      for (size_t c = 0; c < nClasses; c++)
      {
        Vdouble* postProb_i_c = &(*postProb_i)[c];
        postProb_i_c->resize(nStates);
        Vdouble* larray_i_c = &(*larray_i)[c];
        for (size_t x = 0; x < nStates; x++)
        {
          (*postProb_i_c)[x] = (*larray_i_c)[x] / likelihood;
        }
      }
    }
  }
  return postProb;
}

// -----------------------------------------------------------------------------------------

Vdouble DRTreeLikelihoodTools::getPosteriorStateFrequencies(
    const DRTreeLikelihoodInterface& drl,
    int nodeId)
{
  VVVdouble probs = getPosteriorProbabilitiesPerStatePerRate(drl, nodeId);
  Vdouble freqs(drl.getNumberOfStates());
  double sumw = 0, w;
  for (size_t i = 0; i < probs.size(); i++)
  {
    w = drl.likelihoodData().getWeight(i);
    sumw += w;
    for (size_t j = 0; j < drl.getNumberOfClasses(); j++)
    {
      for (size_t k = 0; k < drl.getNumberOfStates(); k++)
      {
        freqs[k] += probs[i][j][k] * w;
      }
    }
  }
  for (size_t k = 0; k < drl.getNumberOfStates(); k++)
  {
    freqs[k] /= sumw;
  }
  return freqs;
}

// -----------------------------------------------------------------------------------------
