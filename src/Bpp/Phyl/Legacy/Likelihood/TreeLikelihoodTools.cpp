// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "TreeLikelihoodTools.h"

using namespace std;
using namespace bpp;

void TreeLikelihoodTools::getAncestralFrequencies(
  const TreeLikelihoodInterface& tl,
  size_t site,
  std::map<int, std::vector<double> >& frequencies,
  bool alsoForLeaves)
{
  int currentId = tl.tree().getRootId();
  vector<double> currentFreqs = tl.getRootFrequencies(tl.getSiteIndex(site));
  getAncestralFrequencies_(tl, tl.getSiteIndex(site), currentId, currentFreqs, frequencies, alsoForLeaves);
}

void TreeLikelihoodTools::getAncestralFrequencies(
  const TreeLikelihoodInterface& tl,
  std::map<int, std::vector<double> >& frequencies,
  bool alsoForLeaves)
{
  size_t n = tl.likelihoodData().getNumberOfDistinctSites();
  size_t ns = tl.getNumberOfStates();
  double sumw = 0, w;
  map<int, vector<double> > siteFrequencies;
  for (size_t i = 0; i < n; ++i)
  {
    w = tl.likelihoodData().getWeight(i);
    sumw += w;
  }
  for (size_t i = 0; i < n; ++i)
  {
    w = tl.likelihoodData().getWeight(i);
    getAncestralFrequencies(tl, i, siteFrequencies, alsoForLeaves);
    // Initialization
    if (i == 0)
    {
      frequencies = siteFrequencies; // Initialize all nodes ids.
      // Now reset to 0:
      for (map<int, vector<double> >::iterator it = frequencies.begin(); it != frequencies.end(); it++)
      {
        VectorTools::fill(it->second, 0.);
      }
    }
    map<int, vector<double> >::iterator it = frequencies.begin();
    map<int, vector<double> >::iterator itSite = siteFrequencies.begin();
    for (size_t j = 0; j < frequencies.size(); ++j)
    {
      for (size_t k = 0; k < ns; ++k)
      {
        it->second[k] += itSite->second[k] * w / sumw;
      }
      it++;
      itSite++;
    }
  }
}

void TreeLikelihoodTools::getAncestralFrequencies_(
  const TreeLikelihoodInterface& tl,
  size_t siteIndex,
  int parentId,
  const std::vector<double>& ancestralFrequencies,
  std::map<int, std::vector<double> >& frequencies,
  bool alsoForLeaves)
{
  if (!tl.tree().isLeaf(parentId) || alsoForLeaves)
    frequencies[parentId] = ancestralFrequencies;
  vector<int> sonsId = tl.tree().getSonsId(parentId);
  for (size_t i = 0; i < sonsId.size(); i++)
  {
    vector<double> sonFrequencies(tl.getNumberOfStates());
    VVdouble pijt = tl.getTransitionProbabilities(sonsId[i], siteIndex);
    for (size_t j = 0; j < tl.getNumberOfStates(); j++)
    {
      double x = 0;
      for (size_t k = 0; k < tl.getNumberOfStates(); k++)
      {
        x += pijt[k][j] * ancestralFrequencies[k];
      }
      sonFrequencies[j] = x;
    }
    getAncestralFrequencies_(tl, siteIndex, sonsId[i], sonFrequencies, frequencies, alsoForLeaves);
  }
}
