//
// File: AbstractLikelihoodTreeCalculation.cpp
// Created by: Julien Dutheil, Laurent Guéguen
// Created on: mardi 23 juin 2015, à 14h 06
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

#include "AbstractLikelihoodTreeCalculation.h"
#include "../PatternTools.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

void AbstractLikelihoodTreeCalculation::setData(const AlignedValuesContainer& sites)
{
  data_.reset(PatternTools::getSequenceSubset(sites, process_->getParametrizablePhyloTree().getRoot(), process_->getParametrizablePhyloTree())) ;

  if (verbose_)
    ApplicationTools::displayTask("Initializing data structure");

  getLikelihoodData().initLikelihoods(*data_, *process_);

// We assume here that all models have the same number of states, and
// that they have the same 'init' method, which is a reasonable
// assumption as long as they share the same alphabet.

  if (verbose_)
    ApplicationTools::displayTaskDone();

  nbSites_         = getLikelihoodData().getNumberOfSites();
  nbDistinctSites_ = getLikelihoodData().getNumberOfDistinctSites();
  nbStates_        = getLikelihoodData().getNumberOfStates();

  if (verbose_)
    ApplicationTools::displayResult("Number of distinct sites", TextTools::toString(nbDistinctSites_));

  initialized_ = true;
  vSites_.resize(nbSites_);
}


/******************************************************************************/

double AbstractLikelihoodTreeCalculation::getLogLikelihood()
{
  for (size_t i = 0; i < nbSites_; ++i){
    vSites_[i] = getLogLikelihoodForASite(i);
  }
  
  sort(vSites_.begin(), vSites_.end());
  double ll = 0;
  for (size_t i = nbSites_; i > 0; i--)
  {
    ll += vSites_[i - 1];
  }
  
  return ll;
}

/******************************************************************************/

double AbstractLikelihoodTreeCalculation::getDLogLikelihood()
{
  // Derivative of the sum is the sum of derivatives:

  for (size_t i = 0; i < nbSites_; ++i)
    vSites_[i] = getDLogLikelihoodForASite(i);
  
  sort(vSites_.begin(), vSites_.end());
  double dl = 0;
  for (size_t i = nbSites_; i > 0; --i)
  {
    dl += vSites_[i - 1];
  }

  return dl;
}

/******************************************************************************/

double AbstractLikelihoodTreeCalculation::getD2LogLikelihood() 
{
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t i = 0; i < nbSites_; ++i)
  {
    dl += getD2LogLikelihoodForASite(i);
  }

  return dl;
}


/******************************************************************************/

double AbstractLikelihoodTreeCalculation::getLikelihoodForASiteIndex(size_t siteindex) const
{
  if (usesLogAtRoot(0))
    return exp(getLogLikelihoodForASiteIndex(siteindex));

  double l = 0;
  
  int Rid=getRootId();
  
  for (size_t c = 0; c < nbClasses_; ++c)
  {
    const VVdouble& lla = getLikelihoodData().getLikelihoodArray(Rid, c, ComputingNode::D0);

    if (!usesLogAtRoot(c))
      for (size_t j = 0; j < nbStates_; ++j)
        l += lla[siteindex][j] * process_->getProbabilityForModel(c);
    else
      for (size_t j = 0; j < nbStates_; ++j)
        l += exp(lla[siteindex][j]) * process_->getProbabilityForModel(c);
  }

  if (l < 0) l = 0; //May happen because of numerical errors.

  return l;
}


/******************************************************************************/

double AbstractLikelihoodTreeCalculation::getLogLikelihoodForASiteIndex(size_t siteindex) const
{
  if (!usesLogAtRoot(0))
    return log(getLikelihoodForASiteIndex(siteindex));

  int Rid=getRootId();

  vector<double> v(nbClasses_);
  
  for (size_t c = 0; c < nbClasses_; ++c)
  {
    const VVdouble& lla = getLikelihoodData().getLikelihoodArray(Rid, c, ComputingNode::D0);

    if (!usesLogAtRoot(c))
      v[c]=VectorTools::logSumExp(VectorTools::log(lla[siteindex])) + log(process_->getProbabilityForModel(c));
    else
      v[c]=VectorTools::logSumExp(lla[siteindex]) + log(process_->getProbabilityForModel(c));
  }

  return VectorTools::logSumExp(v);
}

/******************************************************************************/

double AbstractLikelihoodTreeCalculation::getLikelihoodForASiteIndexForAClass(size_t siteindex, size_t classIndex)
{
  const Vdouble& la = getLikelihoodData().getLikelihoodArray(getRootId(),classIndex, ComputingNode::D0)[siteindex];

  if (!usesLogAtRoot(classIndex))
    return VectorTools::sum(la);
  else
    return VectorTools::sumExp(la);
}

/******************************************************************************/

double AbstractLikelihoodTreeCalculation::getLogLikelihoodForASiteIndexForAClass(size_t siteindex, size_t classIndex)
{
  const Vdouble& la = getLikelihoodData().getLikelihoodArray(getRootId(),classIndex, ComputingNode::D0)[siteindex];

  if (!usesLogAtRoot(classIndex))
    return log(VectorTools::sum(la));
  else
    return VectorTools::logSumExp(la);
}

/******************************************************************************/

double AbstractLikelihoodTreeCalculation::getLikelihoodForASiteIndexForAState(size_t siteindex, int state)
{
  int Rid=getRootId();

  double l = 0;
  
  for (size_t c = 0; c < nbClasses_; ++c)
  {
    const VVdouble& lla = getLikelihoodData().getLikelihoodArray(Rid, c, ComputingNode::D0);

    if (!usesLogAtRoot(c))
      l += lla[siteindex][state] * process_->getProbabilityForModel(c);
    else
      l += exp(lla[siteindex][state]) * process_->getProbabilityForModel(c);
  }

  if (l < 0) l = 0; //May happen because of numerical errors.
  return l;
}

/******************************************************************************/

double AbstractLikelihoodTreeCalculation::getLogLikelihoodForASiteIndexForAState(size_t siteindex, int state)
{
  int Rid=getRootId();
  
  vector<double> v(nbClasses_);

  for (size_t c = 0; c < nbClasses_; ++c)
  {
    const VVdouble& lla = getLikelihoodData().getLikelihoodArray(Rid, c, ComputingNode::D0);

    if (!usesLogAtRoot(c))
      v[c] = log(lla[siteindex][state]) + log(process_->getProbabilityForModel(c));
    else
      v[c] = lla[siteindex][state] + log(process_->getProbabilityForModel(c));
  }

  return VectorTools::logSumExp(v);
}

/******************************************************************************/

double AbstractLikelihoodTreeCalculation::getLikelihoodForASiteIndexForAClassForAState(size_t siteindex, size_t classIndex, int state)
{
  if (!usesLogAtRoot(classIndex))
    return getLikelihoodData().getLikelihoodArray(getRootId(),classIndex, ComputingNode::D0)[siteindex][state];
  else
    return exp(getLikelihoodData().getLikelihoodArray(getRootId(),classIndex, ComputingNode::D0)[siteindex][state]);
}

/******************************************************************************/

double AbstractLikelihoodTreeCalculation::getLogLikelihoodForASiteIndexForAClassForAState(size_t siteindex, size_t classIndex, int state)
{
  if (!usesLogAtRoot(classIndex))
    return log(getLikelihoodData().getLikelihoodArray(getRootId(),classIndex, ComputingNode::D0)[siteindex][state]);
  else
    return getLikelihoodData().getLikelihoodArray(getRootId(),classIndex, ComputingNode::D0)[siteindex][state];
}


/******************************************************************************/

double AbstractLikelihoodTreeCalculation::getDLikelihoodForASiteIndex(size_t siteindex) const
{
  if (nullDLogLikelihood_)
    return 0;

  int Rid=process_->getParametrizablePhyloTree().getNodeIndex(process_->getParametrizablePhyloTree().getRoot());
  
  
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t c = 0; c < nbClasses_; c++)
  {
    const VVdouble& ldla = getLikelihoodData().getLikelihoodArray(Rid, c, ComputingNode::D1);
    for (size_t j = 0; j < nbStates_; ++j)
      dl += ldla[siteindex][j] * process_->getProbabilityForModel(c);
  }
  
  return dl;
}

/******************************************************************************/

double AbstractLikelihoodTreeCalculation::getD2LikelihoodForASiteIndex(size_t siteindex) const
{
  if (nullD2LogLikelihood_)
    return 0;

  int Rid=process_->getParametrizablePhyloTree().getNodeIndex(process_->getParametrizablePhyloTree().getRoot());

  // Derivative of the sum is the sum of derivatives:
  double d2l = 0;
  for (size_t c = 0; c < nbClasses_; c++)
  {
    const VVdouble& d2la = getLikelihoodData().getLikelihoodArray(Rid, c, ComputingNode::D2);
    for (size_t j = 0; j < nbStates_; ++j)
      d2l += d2la[siteindex][j] * process_->getProbabilityForModel(c);
  }

  return d2l;
}


/******************************************************************************/

void AbstractLikelihoodTreeCalculation::displayLikelihoodArray(const VVVdouble& likelihoodArray)
{
  size_t nbClasses = likelihoodArray.size();
  size_t nbSites    = likelihoodArray[0].size();
  size_t nbStates  = likelihoodArray[0][0].size();
  for (size_t c = 0; c < nbClasses; ++c)
  {
    cout << "Model class " << c;
    for (size_t i = 0; i < nbSites; ++i)
    {
      cout << "Site " << i << ":" << endl;
      for (size_t s = 0; s < nbStates; ++s)
      {
        cout << "\t" << likelihoodArray[c][i][s];
      }
      cout << endl;
    }
    cout << endl;
  }
}

/******************************************************************************/
/******************************************************************************/


void AbstractLikelihoodTreeCalculation::getAncestralFrequencies(
  std::map<int, std::vector<double> >& frequencies,
  bool alsoForLeaves)
{
  if (!up2date_)
    throw Exception("AbstractLikelihoodTreeCalculation::getAncestralFrequencies not up to date.");
  
  size_t n = getNumberOfDistinctSites();
  size_t ns = getNumberOfStates();
  double sumw = 0, w;
  map<int, vector<double> > siteFrequencies;
  for (size_t i = 0; i < n; ++i)
  {
    w = getLikelihoodData().getWeight(i);
    sumw += w;
  }

  for (size_t i = 0; i < n; ++i)
  {
    w = getLikelihoodData().getWeight(i);
    getAncestralFrequencies(i, siteFrequencies, alsoForLeaves);
    
    //Initialization
    if (i == 0)
    {
      frequencies = siteFrequencies; //Initialize all nodes ids.
      //Now reset to 0:
      for (auto& it : frequencies)
        VectorTools::fill(it.second, 0.);
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

/******************************************************************************/

void AbstractLikelihoodTreeCalculation::getAncestralFrequencies(
  size_t site,
  std::map<int, std::vector<double> >& frequencies,
  bool alsoForLeaves)
{
  const vector<double>& currentFreqs = getSubstitutionProcess()->getRootFrequencies();
  size_t nbClasses = getNumberOfClasses();

  std::vector<std::map<int, std::vector<double> > > vmfreqcl;

  for (size_t nclass=0; nclass<nbClasses; nclass++)
  {
    int currentId = getLikelihoodData().getRootData(nclass).getId();
    
    std::map<int, std::vector<double> > mfreqcl;
    getAncestralFrequencies_(getSiteIndex(site), nclass, currentId, currentFreqs, mfreqcl, alsoForLeaves);
    vmfreqcl.push_back(mfreqcl);
  }
  
  frequencies.clear();
  frequencies.insert(vmfreqcl[0].begin(), vmfreqcl[0].end());
  double prob0=getSubstitutionProcess()->getProbabilityForModel(0);
  for (auto& it:frequencies)
    it.second*=prob0;

  for (size_t nclass=1; nclass<nbClasses; nclass++)
  {
    std::map<int, std::vector<double> >& mfreqcl=vmfreqcl[nclass];
    double prob=getSubstitutionProcess()->getProbabilityForModel(nclass);
    for (const auto& it:mfreqcl)
      frequencies[it.first]+=it.second * prob;
  }
  
}

/******************************************************************************/

void AbstractLikelihoodTreeCalculation::getAncestralFrequencies_(
  size_t siteIndex,
  size_t classIndex,
  int parentId,
  const std::vector<double>& ancestralFrequencies,
  std::map<int,std::vector<double> >& frequencies,
  bool alsoForLeaves)
{
  const AbstractLikelihoodNode& parentNode = getLikelihoodData().getNodeData(parentId, classIndex);

  if (!parentNode.isLeaf() || alsoForLeaves)
    frequencies[parentId] = ancestralFrequencies;

  vector<unsigned int> sonsId = parentNode.getSonsId();
  for (size_t i = 0; i < sonsId.size(); i++)
  {
    vector<double> sonFrequencies(getNumberOfStates());
    const bpp::Matrix<double>& pijt = getSubstitutionProcess()->getTransitionProbabilities(sonsId[i], classIndex);
    for (size_t j = 0; j < getNumberOfStates(); j++)
    {
      double x = 0;
      for (size_t k = 0; k < getNumberOfStates(); k++)
      {
        x += pijt(k,j) * ancestralFrequencies[k];
      }
      sonFrequencies[j] = x;
    }
    getAncestralFrequencies_(siteIndex, classIndex, sonsId[i], sonFrequencies, frequencies, alsoForLeaves);
  }
}



