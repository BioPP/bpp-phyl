//
// File: ProductOfAlignedPhyloLikelihood.cpp
// Created by: Laurent Guéguen
// Created on: vendredi 24 octobre 2014, à 17h 33
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

#include "ProductOfAlignedPhyloLikelihood.h"

using namespace bpp;
using namespace std;

ProductOfAlignedPhyloLikelihood::ProductOfAlignedPhyloLikelihood(PhyloLikelihoodContainer* pC) :
  AbstractPhyloLikelihood(),
  AbstractAlignedPhyloLikelihood(0),
  SetOfAlignedPhyloLikelihood(pC)
{
}

ProductOfAlignedPhyloLikelihood::ProductOfAlignedPhyloLikelihood(PhyloLikelihoodContainer* pC, const std::vector<size_t>& nPhylo) :
  AbstractPhyloLikelihood(),
  AbstractAlignedPhyloLikelihood(0),
  SetOfAlignedPhyloLikelihood(pC, nPhylo)
{
}

ProductOfAlignedPhyloLikelihood::ProductOfAlignedPhyloLikelihood(const ProductOfAlignedPhyloLikelihood& sd) :
  AbstractPhyloLikelihood(sd),
  AbstractAlignedPhyloLikelihood(sd),
  SetOfAlignedPhyloLikelihood(sd)
{
}

ProductOfAlignedPhyloLikelihood& ProductOfAlignedPhyloLikelihood::operator=(const ProductOfAlignedPhyloLikelihood& sd)
{
  SetOfAlignedPhyloLikelihood::operator=(sd);
  
  return *this;
}


double ProductOfAlignedPhyloLikelihood::getDLogLikelihoodForASite(size_t site) const
{
  double x=0;
  const std::vector<size_t>& nPhylo=getNumbersOfPhyloLikelihoods();
        
  for (size_t i=0; i<nPhylo.size(); i++)
    x += getAbstractPhyloLikelihood(nPhylo[i])->getDLogLikelihoodForASite(site);
  return x;
}


double ProductOfAlignedPhyloLikelihood::getD2LogLikelihoodForASite(size_t site) const
{
  double x=0;
  const std::vector<size_t>& nPhylo=getNumbersOfPhyloLikelihoods();
        
  for (size_t i=0; i<nPhylo.size(); i++)
    x += getAbstractPhyloLikelihood(nPhylo[i])->getD2LogLikelihoodForASite(site);
  return x;
}



/******************************************************************************/

double ProductOfAlignedPhyloLikelihood::getLogLikelihood() const
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

double ProductOfAlignedPhyloLikelihood::getDLogLikelihood() const
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

double ProductOfAlignedPhyloLikelihood::getD2LogLikelihood() const
{
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
