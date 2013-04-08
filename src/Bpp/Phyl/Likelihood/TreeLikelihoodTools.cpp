//
// File: TreeLikelihoodTools.cpp
// Created by: Julien Dutheil
// Created on: Tue Jun 30 12:25 2009
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

#include "TreeLikelihoodTools.h"

using namespace std;
using namespace bpp;

void TreeLikelihoodTools::getAncestralFrequencies(
        const TreeLikelihood& tl,
        size_t site,
        std::map<int, std::vector<double> >& frequencies,
        bool alsoForLeaves) throw (Exception)
{
  int currentId = tl.getTree().getRootId();
  vector<double> currentFreqs = tl.getRootFrequencies(tl.getSiteIndex(site));
  getAncestralFrequencies_(tl, tl.getSiteIndex(site), currentId, currentFreqs, frequencies, alsoForLeaves); 
}

void TreeLikelihoodTools::getAncestralFrequencies(
        const TreeLikelihood& tl,
        std::map<int, std::vector<double> >& frequencies,
        bool alsoForLeaves) throw (Exception)
{
  size_t n = tl.getLikelihoodData()->getNumberOfDistinctSites();
  size_t ns = tl.getNumberOfStates();
  double sumw = 0, w;
  map<int, vector<double> > siteFrequencies;
  for (size_t i = 0; i < n; ++i)
  {
    w = tl.getLikelihoodData()->getWeight(i);
    sumw += w;
  }
  for (size_t i = 0; i < n; ++i)
  {
    w = tl.getLikelihoodData()->getWeight(i);
    getAncestralFrequencies(tl, i, siteFrequencies, alsoForLeaves);
    //Initialization
    if (i == 0)
    {
      frequencies = siteFrequencies; //Initialize all nodes ids.
      //Now reset to 0:
      for (map<int, vector<double> >::iterator it = frequencies.begin(); it != frequencies.end(); it++)
        VectorTools::fill(it->second, 0.);
    }
    map<int, vector<double> >::iterator it = frequencies.begin();
    map<int, vector<double> >::iterator itSite = siteFrequencies.begin();
    for (size_t j = 0; j < frequencies.size(); ++j)
    {
      for (size_t k = 0; k < ns; ++k)
        it->second[k] += itSite->second[k] * w / sumw;
      it++;
      itSite++;
    }
  }
}

void TreeLikelihoodTools::getAncestralFrequencies_(
        const TreeLikelihood& tl,
        size_t siteIndex,
        int parentId,
        const std::vector<double>& ancestralFrequencies,
        std::map<int,std::vector<double> >& frequencies,
        bool alsoForLeaves) throw (Exception)
{
  if (!tl.getTree().isLeaf(parentId) || alsoForLeaves)
    frequencies[parentId] = ancestralFrequencies;
  vector<int> sonsId = tl.getTree().getSonsId(parentId);
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

