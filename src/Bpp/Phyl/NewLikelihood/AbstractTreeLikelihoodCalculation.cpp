//
// File: AbstractTreeLikelihoodCalculation.cpp
// Created by: Julien Dutheil
// Created on: Tue July 23 10:50 2013
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

#include "AbstractTreeLikelihoodCalculation.h"
#include "../PatternTools.h"

using namespace bpp;
using namespace std;
using namespace newlik;

/******************************************************************************/

void AbstractTreeLikelihoodCalculation::setData(const SiteContainer& sites)
{
  const TreeTemplate<Node>& tt = dynamic_cast<const TreeTemplate<Node>&>(process_->getTree());

  data_.reset(PatternTools::getSequenceSubset(sites, *tt.getRootNode()));

  if (verbose_)
    ApplicationTools::displayTask("Initializing data structure");
  getLikelihoodData()->initLikelihoods(*data_, *process_);

// We assume here that all models have the same number of states, and
// that they have the same 'init' method, which is a reasonable
// assumption as long as they share the same alphabet.

  if (verbose_)
    ApplicationTools::displayTaskDone();

  nbSites_         = getLikelihoodData()->getNumberOfSites();
  nbDistinctSites_ = getLikelihoodData()->getNumberOfDistinctSites();
  nbStates_        = getLikelihoodData()->getNumberOfStates();

  if (verbose_)
    ApplicationTools::displayResult("Number of distinct sites", TextTools::toString(nbDistinctSites_));

  initialized_ = true;
  resetToCompute();
}


/******************************************************************************/

double AbstractTreeLikelihoodCalculation::getLogLikelihood()
{
//  if (computeLikelihoods_)
    computeTreeLikelihood();
  
  vector<double> la(nbSites_);
  for (size_t i = 0; i < nbSites_; ++i)
  {
    la[i] = log(getLikelihoodForASite(i));
  }

  sort(la.begin(), la.end());
  double ll = 0;
  for (size_t i = nbSites_; i > 0; i--)
  {
    ll += la[i - 1];
  }
  return ll;
}

// /******************************************************************************/

// void DoubleRecursiveTreeLikelihoodCalculation::setData(const SiteContainer& sites) throw (Exception)
// {
//   const TreeTemplate<Node>& tt = dynamic_cast<const TreeTemplate<Node>&>(process_->getTree());

//   data_.reset(PatternTools::getSequenceSubset(sites, *tt.getRootNode()));
//   if (verbose_)
//     ApplicationTools::displayTask("Initializing data structure");
//   likelihoodData_->initLikelihoods(*data_, *process_);
// // We assume here that all models have the same number of states, and
// // that they have the same 'init' method, which is a reasonable
// // assumption as long as they share the same alphabet.

//   if (verbose_)
//     ApplicationTools::displayTaskDone();

//   nbSites_         = likelihoodData_->getNumberOfSites();
//   nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();
//   nbStates_        = likelihoodData_->getNumberOfStates();

//   if (verbose_)
//     ApplicationTools::displayResult("Number of distinct sites", TextTools::toString(nbDistinctSites_));

//   initialized_ = true;
//   compNId_ = -1;
  
//   //Recompute likelihood:
//   computeTreeLikelihood();
// }

/******************************************************************************/

double AbstractTreeLikelihoodCalculation::getDLogLikelihood()
{
  if (computeLikelihoodsD1_){
    throw Exception("getDLogLikelihood() : DLogLikelihood not computed");
    computeLikelihoodsD1_=false;
  }

  // Derivative of the sum is the sum of derivatives:
  vector<double> dla(nbSites_);
  for (size_t i = 0; i < nbSites_; ++i)
  {
    dla[i] = getDLogLikelihoodForASite(i);
  }
  
  sort(dla.begin(), dla.end());
  double dl = 0;
  for (size_t i = nbSites_; i > 0; --i)
  {
    dl += dla[i - 1];
  }

  return dl;
}

/******************************************************************************/

double AbstractTreeLikelihoodCalculation::getD2LogLikelihood() 
{
  if (computeLikelihoodsD2_){
    throw Exception("getD2LogLikelihood() : D2LogLikelihood not computed");
    computeLikelihoodsD2_=false;
  }
  
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  for (size_t i = 0; i < nbSites_; ++i)
  {
    dl += getD2LogLikelihoodForASite(i);
  }

  return dl;
}


/******************************************************************************/

void AbstractTreeLikelihoodCalculation::displayLikelihoodArray(const VVVdouble& likelihoodArray)
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


void AbstractTreeLikelihoodCalculation::getAncestralFrequencies(
  std::map<int, std::vector<double> >& frequencies,
  bool alsoForLeaves)
{
  size_t n = getNumberOfDistinctSites();
  size_t ns = getNumberOfStates();
  double sumw = 0, w;
  map<int, vector<double> > siteFrequencies;
  for (size_t i = 0; i < n; ++i)
  {
    w = getLikelihoodData()->getWeight(i);
    sumw += w;
  }

  for (size_t i = 0; i < n; ++i)
  {
    w = getLikelihoodData()->getWeight(i);
    getAncestralFrequencies(i, siteFrequencies, alsoForLeaves);
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

/******************************************************************************/

void AbstractTreeLikelihoodCalculation::getAncestralFrequencies(
  size_t site,
  std::map<int, std::vector<double> >& frequencies,
  bool alsoForLeaves)
{
  int currentId = getSubstitutionProcess()->getTree().getRootId();
  const vector<double>& currentFreqs = getSubstitutionProcess()->getRootFrequencies();
  size_t nbClasses = getNumberOfClasses();

  std::vector<std::map<int, std::vector<double> > > vmfreqcl;

  for (size_t nclass=0; nclass<nbClasses; nclass++)
  {
    std::map<int, std::vector<double> > mfreqcl;
    getAncestralFrequencies_(getSiteIndex(site), nclass, currentId, currentFreqs, mfreqcl, alsoForLeaves);
    vmfreqcl.push_back(mfreqcl);
  }

  frequencies.clear();
  frequencies.insert(vmfreqcl[0].begin(), vmfreqcl[0].end());

  for (size_t nclass=1; nclass<nbClasses; nclass++)
  {
    std::map<int, std::vector<double> >& mfreqcl=vmfreqcl[nclass];
    double prob=getSubstitutionProcess()->getProbabilityForModel(nclass);
    
    for (std::map<int, std::vector<double> >::const_iterator it=mfreqcl.begin(); it!= mfreqcl.end(); it++)
      frequencies[it->first]+=it->second * prob;
  }
  
}

/******************************************************************************/

void AbstractTreeLikelihoodCalculation::getAncestralFrequencies_(
  size_t siteIndex,
  size_t classIndex,
  int parentId,
  const std::vector<double>& ancestralFrequencies,
  std::map<int,std::vector<double> >& frequencies,
  bool alsoForLeaves)
{
  if (!getSubstitutionProcess()->getTree().isLeaf(parentId) || alsoForLeaves)
    frequencies[parentId] = ancestralFrequencies;

  vector<int> sonsId = getSubstitutionProcess()->getTree().getSonsId(parentId);
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



