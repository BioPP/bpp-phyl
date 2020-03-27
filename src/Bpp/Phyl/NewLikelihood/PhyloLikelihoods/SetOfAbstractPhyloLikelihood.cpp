//
// File: SetOfAbstractPhyloLikelihood.cpp
// Created by: Laurent Guéguen
// Created on: jeudi 14 mai 2015, à 17h 07
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

#include "SetOfAbstractPhyloLikelihood.h"

using namespace bpp;
using namespace std;

SetOfAbstractPhyloLikelihood::SetOfAbstractPhyloLikelihood(Context& context, PhyloLikelihoodContainer* pC, const std::string& prefix) :
  AbstractPhyloLikelihood(context),
  AbstractParametrizable(prefix),
  pPhyloCont_(pC),
  nPhylo_()
{
}

 
SetOfAbstractPhyloLikelihood::SetOfAbstractPhyloLikelihood(const SetOfAbstractPhyloLikelihood& sd) :
  AbstractPhyloLikelihood(sd),
  AbstractParametrizable(sd),
  pPhyloCont_(sd.pPhyloCont_),
  nPhylo_(sd.nPhylo_)
{
}

bool SetOfAbstractPhyloLikelihood::addPhyloLikelihood(size_t nPhyl)
{
  const AbstractPhyloLikelihood* aPL=getAbstractPhyloLikelihood(nPhyl);

  if (aPL!=NULL){
    nPhylo_.push_back(nPhyl);
    includeParameters_(aPL->getParameters());
    return true;
  }
  return false;
}

void SetOfAbstractPhyloLikelihood::addAllPhyloLikelihoods()
{
  vector<size_t> vPhyl=getPhyloContainer()->getNumbersOfPhyloLikelihoods();

  for (size_t i=0; i<vPhyl.size(); i++)
    addPhyloLikelihood(vPhyl[i]);
}

ParameterList SetOfAbstractPhyloLikelihood::getBranchLengthParameters() const
{
  ParameterList pl;
  for (size_t i=0; i<nPhylo_.size(); i++)
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getBranchLengthParameters());
  
  return pl;
}
      
ParameterList SetOfAbstractPhyloLikelihood::getSubstitutionModelParameters() const
{
  ParameterList pl;
  for (size_t i=0; i<nPhylo_.size(); i++)
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getSubstitutionModelParameters());

  return pl;
}

ParameterList SetOfAbstractPhyloLikelihood::getRateDistributionParameters() const
{
  ParameterList pl;
  for (size_t i=0; i<nPhylo_.size(); i++)
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getRateDistributionParameters());

  return pl;
}

ParameterList SetOfAbstractPhyloLikelihood::getRootFrequenciesParameters() const
{
  ParameterList pl;
  for (size_t i=0; i<nPhylo_.size(); i++)
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getRootFrequenciesParameters());

  return pl;
}

ParameterList SetOfAbstractPhyloLikelihood::getNonDerivableParameters() const
{
  ParameterList pl;
  for (size_t i=0; i<nPhylo_.size(); i++)
    pl.includeParameters(getPhyloLikelihood(nPhylo_[i])->getNonDerivableParameters());

  return pl;
}

