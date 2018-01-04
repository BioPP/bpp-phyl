//
// File: OneProcessSequencePhyloLikelihood.cpp
// Created by: Julien Dutheil
// Created on: mardi 28 avril 2015, à 13h 11
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

#include "OneProcessSequencePhyloLikelihood.h"
#include "../RecursiveLikelihoodTreeCalculation.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

OneProcessSequencePhyloLikelihood::OneProcessSequencePhyloLikelihood(
  OneProcessSequenceEvolution& evol,
  size_t nSeqEvol,
  bool verbose,
  bool patterns) :
  AbstractPhyloLikelihood(),
  AbstractAlignedPhyloLikelihood(0),
  AbstractSequencePhyloLikelihood(evol, nSeqEvol),
  mSeqEvol_(evol),
  tlComp_()
{
  SubstitutionProcess& sp = evol.getSubstitutionProcess();

  tlComp_ = std::unique_ptr<LikelihoodTreeCalculation>(new RecursiveLikelihoodTreeCalculation(&sp, true, patterns));
}

/******************************************************************************/

OneProcessSequencePhyloLikelihood::OneProcessSequencePhyloLikelihood(
  const AlignedValuesContainer& data,
  OneProcessSequenceEvolution& evol,
  size_t nSeqEvol,
  size_t nData,
  bool verbose,
  bool patterns) :
  AbstractPhyloLikelihood(),
  AbstractAlignedPhyloLikelihood(data.getNumberOfSites()),
  AbstractSequencePhyloLikelihood(evol, nSeqEvol),
  mSeqEvol_(evol),
  tlComp_()
{
  SubstitutionProcess& sp = evol.getSubstitutionProcess();

  tlComp_ = std::unique_ptr<LikelihoodTreeCalculation>(new RecursiveLikelihoodTreeCalculation(&sp, true, patterns));

  setData(data, nData);
}

/******************************************************************************/

VVdouble OneProcessSequencePhyloLikelihood::getLikelihoodPerSitePerState() const
{
  updateLikelihood();
  computeLikelihood();

  VVdouble l(getNumberOfSites());
  for (size_t i = 0; i < l.size(); ++i)
  {
    Vdouble* l_i = &l[i];
    l_i->resize(getNumberOfStates());
    for (size_t x = 0; x < l_i->size(); ++x)
    {
      (*l_i)[x] = tlComp_->getLikelihoodForASiteForAState(i, static_cast<int>(x));
    }
  }
  return l;
}

/******************************************************************************/

VVdouble OneProcessSequencePhyloLikelihood::getLikelihoodPerSitePerClass() const
{
  updateLikelihood();
  computeLikelihood();

  VVdouble l(getNumberOfSites());
  for (size_t i = 0; i < l.size(); ++i)
  {
    Vdouble* l_i = &l[i];
    l_i->resize(getNumberOfClasses());
    for (size_t c = 0; c < l_i->size(); ++c)
    {
      (*l_i)[c] = tlComp_->getLikelihoodForASiteForAClass(i, c);
    }
  }
  return l;
}

/******************************************************************************/

VVVdouble OneProcessSequencePhyloLikelihood::getLikelihoodPerSitePerClassPerState() const
{
  updateLikelihood();
  computeLikelihood();

  VVVdouble l(getNumberOfSites());
  for (size_t i = 0; i < l.size(); ++i)
  {
    VVdouble* l_i = &l[i];
    l_i->resize(getNumberOfClasses());
    for (size_t c = 0; c < l_i->size(); ++c)
    {
      Vdouble* l_ic = &(*l_i)[c];
      l_ic->resize(getNumberOfStates());
      for (size_t x = 0; x < l_ic->size(); ++x)
      {
        (*l_ic)[x] = tlComp_->getLikelihoodForASiteForAClassForAState(i, c, static_cast<int>(x));
      }
    }
  }
  return l;
}

/******************************************************************************/

VVdouble OneProcessSequencePhyloLikelihood::getPosteriorProbabilitiesPerClass() const
{
  updateLikelihood();
  computeLikelihood();

  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  VVdouble pb = getLikelihoodPerSitePerClass();
  Vdouble l = getLikelihoodPerSite();

  for (size_t i = 0; i < nbSites; ++i)
  {
    for (size_t j = 0; j < nbClasses; ++j)
    {
      pb[i][j] = pb[i][j] * getSubstitutionProcess().getProbabilityForModel(j) / l[i];
    }
  }
  return pb;
}

/******************************************************************************/

vector<size_t> OneProcessSequencePhyloLikelihood::getClassWithMaxPostProbPerSite() const
{
  updateLikelihood();
  computeLikelihood();

  size_t nbSites = getNumberOfSites();
  VVdouble l = getLikelihoodPerSitePerClass();
  vector<size_t> classes(nbSites);
  for (size_t i = 0; i < nbSites; ++i)
  {
    classes[i] = VectorTools::whichMax<double>(l[i]);
  }
  return classes;
}


/******************************************************************************/

Vdouble OneProcessSequencePhyloLikelihood::getPosteriorRatePerSite() const
{
  updateLikelihood();
  computeLikelihood();

  size_t nbSites   = getNumberOfSites();
  size_t nbClasses = getNumberOfClasses();
  VVdouble pb = getLikelihoodPerSitePerClass();
  Vdouble l  = getLikelihoodPerSite();
  Vdouble rates(nbSites, 0.);
  for (size_t i = 0; i < nbSites; i++)
  {
    for (size_t j = 0; j < nbClasses; j++)
    {
      rates[i] += (pb[i][j] / l[i]) * getSubstitutionProcess().getProbabilityForModel(j) *  getSubstitutionProcess().getProbabilityForModel(j);
    }
  }
  return rates;
}

/******************************************************************************/

void OneProcessSequencePhyloLikelihood::computeDLogLikelihood_(const string& variable) const
{
  // check it is a "BrLen" variable

  if (!hasParameter(variable) || (variable.compare(0,5,"BrLen")!=0))
    return;
  
  // Get the node with the branch whose length must be derivated:
  Vuint vbrId;
  
  try {
    vbrId.push_back((unsigned int)atoi(variable.substr(5).c_str()));
  }
  catch (std::exception const& e)
  {
    return;
  }

  tlComp_->computeTreeDLogLikelihood(vbrId);
  dValues_[variable]= std::nan("");
}

/******************************************************************************/

void OneProcessSequencePhyloLikelihood::computeD2LogLikelihood_(const string& variable) const
{
  // check it is a "BrLen" variable

  if (!hasParameter(variable) || (variable.compare(0,5,"BrLen")!=0))
    return;

  // Get the node with the branch whose length must be derivated:
  Vuint vbrId;
  
  try {
    vbrId.push_back((unsigned int)atoi(variable.substr(5).c_str()));
  }
  catch (std::exception const& e)
  {
    return;
  }

  tlComp_->computeTreeD2LogLikelihood(vbrId);
  d2Values_[variable]= std::nan("");
}
