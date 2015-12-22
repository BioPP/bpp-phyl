//
// File: MixtureOfAlignedPhyloLikelihood.cpp
// Created by: Laurent Guéguen
// Created on: lundi 26 octobre 2015, à 21h 56
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

#include "MixtureOfAlignedPhyloLikelihood.h"

using namespace bpp;
using namespace std;

MixtureOfAlignedPhyloLikelihood::MixtureOfAlignedPhyloLikelihood(PhyloLikelihoodContainer* pC, const std::vector<size_t>& nPhylo) :
  AbstractPhyloLikelihood(),
  AbstractAlignedPhyloLikelihood(0),
  SetOfAlignedPhyloLikelihood(pC, nPhylo, ""),
  simplex_(getNumbersOfPhyloLikelihoods().size(), 1, false, "Mixture.")
{
  addParameters_(simplex_.getParameters());
  update();  
}


MixtureOfAlignedPhyloLikelihood::MixtureOfAlignedPhyloLikelihood(const MixtureOfAlignedPhyloLikelihood& sd) :
  AbstractPhyloLikelihood(sd),
  AbstractAlignedPhyloLikelihood(sd),
  SetOfAlignedPhyloLikelihood(sd),
  simplex_(sd.simplex_)
{
}

MixtureOfAlignedPhyloLikelihood& MixtureOfAlignedPhyloLikelihood::operator=(const MixtureOfAlignedPhyloLikelihood& sd)
{
  SetOfAlignedPhyloLikelihood::operator=(sd);
  simplex_=sd.simplex_;
  
  return *this;
}

void MixtureOfAlignedPhyloLikelihood::setPhyloProb(const Simplex& si)
{
  simplex_.setFrequencies(si.getFrequencies());
  matchParametersValues(simplex_.getParameters());
}


/******************************************************************************/

void MixtureOfAlignedPhyloLikelihood::fireParameterChanged(const ParameterList& parameters)
{
  simplex_.matchParametersValues(parameters);
  SetOfAbstractPhyloLikelihood::fireParameterChanged(parameters);
  
  update();
}

/******************************************************************************/

double MixtureOfAlignedPhyloLikelihood::getLogLikelihood() const
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

double MixtureOfAlignedPhyloLikelihood::getLikelihoodForASite(size_t site) const
{
  updateLikelihood();
  computeLikelihood();

  double x = 0;
  const std::vector<size_t>& nPhylo=getNumbersOfPhyloLikelihoods();
        
  for (size_t i=0; i<nPhylo.size(); i++)
  {
    x += getAbstractPhyloLikelihood(nPhylo[i])->getLikelihoodForASite(site) * simplex_.prob(i);
  }

  return x;
}

/******************************************************************************/

double MixtureOfAlignedPhyloLikelihood::getLogLikelihoodForASite(size_t site) const
{
  updateLikelihood();
  computeLikelihood();

  const std::vector<size_t>& nPhylo=getNumbersOfPhyloLikelihoods();
  vector<double> v(nPhylo.size());
        
  for (size_t i=0; i<nPhylo.size(); i++)
    v[i] = getAbstractPhyloLikelihood(nPhylo[i])->getLogLikelihoodForASite(site);
  
  return VectorTools::logSumExp(v, simplex_.getFrequencies());
}


/******************************************************************************/

double MixtureOfAlignedPhyloLikelihood::getDLogLikelihood() const
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

double MixtureOfAlignedPhyloLikelihood::getDLogLikelihoodForASite(size_t site) const
{
  const std::vector<size_t>& nPhylo=getNumbersOfPhyloLikelihoods();
  vector<double> vD(nPhylo.size());
        
  for (size_t i=0; i<nPhylo.size(); i++)
    vD[i] = getAbstractPhyloLikelihood(nPhylo[i])->getDLogLikelihoodForASite(site);

  return VectorTools::logSumExp(vD, simplex_.getFrequencies());
}

/******************************************************************************/

double MixtureOfAlignedPhyloLikelihood::getD2LogLikelihood() const
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

double MixtureOfAlignedPhyloLikelihood::getD2LogLikelihoodForASite(size_t site) const
{
  const std::vector<size_t>& nPhylo=getNumbersOfPhyloLikelihoods();
  vector<double> vD2(nPhylo.size());
        
  for (size_t i=0; i<nPhylo.size(); i++)
    vD2[i] = getAbstractPhyloLikelihood(nPhylo[i])->getD2LogLikelihoodForASite(site);

  return VectorTools::logSumExp(vD2, simplex_.getFrequencies());
}


