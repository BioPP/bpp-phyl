//
// File: MixtureProcessPhyloLikelihood.cpp
// Created by: Laurent Guéguen
// Created on: vendredi 12 juillet 2013, à 14h 55
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

#include "MixtureProcessPhyloLikelihood.h"

#include <Bpp/Numeric/VectorTools.h>

using namespace std;
using namespace bpp;

/******************************************************************************/

MixtureProcessPhyloLikelihood::MixtureProcessPhyloLikelihood(
  const SiteContainer& data,
  MixtureSequenceEvolution& processSeqEvol,
  size_t nSeqEvol,
  size_t nData,
  bool verbose,
  bool patterns) :
  AbstractPhyloLikelihood(),
  AbstractAlignedPhyloLikelihood(data.getNumberOfSites()),
  MultiProcessSequencePhyloLikelihood(data, processSeqEvol, nSeqEvol, nData, verbose, patterns),
  mSeqEvol_(processSeqEvol)
{
}

/******************************************************************************/

double MixtureProcessPhyloLikelihood::getLogLikelihood() const
{
  updateLikelihood();
  computeLikelihood();

  vector<double> la(nbSites_);
  for (size_t i = 0; i < nbSites_; ++i)
    {
      la[i] = getLogLikelihoodForASite(i);
    }

  sort(la.begin(), la.end());
  double ll = 0;
  for (size_t i = nbSites_; i > 0; i--)
  {
    ll += la[i - 1];
  }
  return ll;
}

/******************************************************************************/

double MixtureProcessPhyloLikelihood::getDLogLikelihoodForASite(size_t site) const
{
  Vdouble vD;

  for (size_t i = 0; i < vpTreelik_.size(); i++)
    vD.push_back(vpTreelik_[i]->getDLogLikelihoodForASite(site));

  return VectorTools::logSumExp(vD,mSeqEvol_.getSubProcessProbabilities());
}

/******************************************************************************/

double MixtureProcessPhyloLikelihood::getD2LogLikelihoodForASite(size_t site) const
{
  Vdouble vD2;

  for (size_t i = 0; i < vpTreelik_.size(); i++)
    vD2.push_back(vpTreelik_[i]->getD2LogLikelihoodForASite(site));

  return VectorTools::logSumExp(vD2,mSeqEvol_.getSubProcessProbabilities());
}

/******************************************************************************/

double MixtureProcessPhyloLikelihood::getDLogLikelihood() const
{
  vector<double> la(nbSites_);
  for (size_t i = 0; i < nbSites_; ++i)
  {
    la[i] = getDLogLikelihoodForASite(i);
  }
  sort(la.begin(), la.end());
  double ll = 0;
  for (size_t i = nbSites_; i > 0; i--)
  {
    ll += la[i - 1];
  }
  return ll;
}

/******************************************************************************/

double MixtureProcessPhyloLikelihood::getD2LogLikelihood() const
{
  // Derivative of the sum is the sum of derivatives:
  vector<double> la(nbSites_);
  for (size_t i = 0; i < nbSites_; ++i)
  {
    la[i] = getD2LogLikelihoodForASite(i);
  }
  sort(la.begin(), la.end());
  double ll = 0;
  for (size_t i = nbSites_; i > 0; i--)
  {
    ll += la[i - 1];
  }
  return ll;
}

/******************************************************************************/

double MixtureProcessPhyloLikelihood::getLikelihoodForASite(size_t site) const
{
  updateLikelihood();
  computeLikelihood();

  double x = 0;
  for (size_t i = 0; i < vpTreelik_.size(); i++)
  {
    x += vpTreelik_[i]->getLikelihoodForASite(site) * mSeqEvol_.getSubProcessProb(i);
  }

  return x;
}

/******************************************************************************/

double MixtureProcessPhyloLikelihood::getLogLikelihoodForASite(size_t site) const
{
  updateLikelihood();
  computeLikelihood();

  vector<double> v(vpTreelik_.size());
  
  for (size_t i = 0; i < vpTreelik_.size(); i++)
    v[i]= vpTreelik_[i]->getLogLikelihoodForASite(site) + log(mSeqEvol_.getSubProcessProb(i));

  return VectorTools::logSumExp(v);
}


/******************************************************************************/

VVdouble MixtureProcessPhyloLikelihood::getPosteriorProbabilitiesForEachSiteForEachProcess() const
{
  updateLikelihood();
  computeLikelihood();

  size_t nbProcess = getNumberOfSubstitutionProcess();

  VVdouble pb = getLikelihoodForEachSiteForEachProcess();
  Vdouble l = getLikelihoodForEachSite();

  for (size_t i = 0; i < nbSites_; ++i)
  {
    for (size_t j = 0; j < nbProcess; ++j)
    {
      pb[i][j] = pb[i][j] * getSubProcessProb(j) / l[i];
    }
  }
  return pb;
}


/******************************************************************************/

void MixtureProcessPhyloLikelihood::computeD2LogLikelihood_(const std::string& variable) const
{
  for (size_t i = 0; i < vpTreelik_.size(); i++)
  {
    vpTreelik_[i]->computeTreeD2LogLikelihood(variable);
  }
}

/******************************************************************************/

void MixtureProcessPhyloLikelihood::computeDLogLikelihood_(const std::string& variable) const
{
  for (size_t i = 0; i < vpTreelik_.size(); i++)
  {
    vpTreelik_[i]->computeTreeDLogLikelihood(variable);
  }
}

/******************************************************************************/
