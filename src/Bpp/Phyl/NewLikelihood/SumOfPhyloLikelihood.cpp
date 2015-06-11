//
// File: SumOfPhyloLikelihood.cpp
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

#include "SumOfPhyloLikelihood.h"

using namespace bpp;
using namespace std;

SumOfPhyloLikelihood::SumOfPhyloLikelihood() :
  AbstractPhyloLikelihood(),
  AbstractParametrizable(""),
  mSDP_()
{
}

SumOfPhyloLikelihood::SumOfPhyloLikelihood(std::map<size_t, PhyloLikelihood*>& mSDP) :
  AbstractPhyloLikelihood(),  
  AbstractParametrizable(""),
  mSDP_()
{
  for (std::map<size_t, PhyloLikelihood*>::const_iterator it=mSDP.begin(); it != mSDP.end(); it++)
  {
    addPhylolikelihood(it->first, it->second);
    includeParameters_(mSDP_[it->first]->getParameters());
  }
  
}

SumOfPhyloLikelihood::~SumOfPhyloLikelihood()
{
  for (std::map<size_t, AbstractPhyloLikelihood*>::iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    delete it->second;
}

SumOfPhyloLikelihood::SumOfPhyloLikelihood(const SumOfPhyloLikelihood& sd) :
  AbstractPhyloLikelihood(sd),
  AbstractParametrizable(sd),
  mSDP_()
{
  for (std::map<size_t, AbstractPhyloLikelihood*>::const_iterator it=sd.mSDP_.begin(); it != sd.mSDP_.end(); it++)
    addPhylolikelihood(it->first, it->second->clone());
}

SumOfPhyloLikelihood& SumOfPhyloLikelihood::operator=(const SumOfPhyloLikelihood& sd)
{
  AbstractPhyloLikelihood::operator=(sd);
  AbstractParametrizable::operator=(sd);
  
  for (std::map<size_t, AbstractPhyloLikelihood*>::iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    delete it->second;
  
  for (std::map<size_t, AbstractPhyloLikelihood*>::const_iterator it=sd.mSDP_.begin(); it != sd.mSDP_.end(); it++)
    addPhylolikelihood(it->first, it->second->clone());
  
  return *this;
}


void SumOfPhyloLikelihood::addPhylolikelihood(size_t nPhyl, PhyloLikelihood* SDP)
{
  if (mSDP_.find(nPhyl)!=mSDP_.end())
    throw Exception("SumOfPhyloLikelihood::addPhylolikelihood: map number already used : " + TextTools::toString(nPhyl));
  
  if (dynamic_cast<AbstractPhyloLikelihood*>(SDP)!=NULL){
    mSDP_[nPhyl]=dynamic_cast<AbstractPhyloLikelihood*>(SDP);
    includeParameters_(SDP->getParameters());
  }

  update();
}


void SumOfPhyloLikelihood::setParameters(const ParameterList& parameters) throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}

double SumOfPhyloLikelihood::getFirstOrderDerivative(const std::string& variable) const throw (Exception)
{
  // patch, to be fixed properly later
  throw Exception("Derivative is not implemented for " + variable + " parameter.");
  
  double x=0;
  for (std::map<size_t, AbstractPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    if (it->second->hasParameter(variable))
      x+=it->second->getFirstOrderDerivative(variable);

  return x;
}

double SumOfPhyloLikelihood::getSecondOrderDerivative(const std::string& variable) const throw (Exception)
{
  // patch, to be fixed properly later
  throw Exception("Derivative is not implemented for " + variable + " parameter.");

  double x=0;
  for (std::map<size_t, AbstractPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    if (it->second->hasParameter(variable))
      x+=it->second->getSecondOrderDerivative(variable);
  
  return x;
}

std::vector<size_t> SumOfPhyloLikelihood::getNumbersOfPhyloLikelihoods() const
{
  vector<size_t> vNum;
  for (std::map<size_t, AbstractPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    vNum.push_back(it->first);

  return vNum;
}

ParameterList SumOfPhyloLikelihood::getBranchLengthParameters() const
{
  ParameterList pl;
  for (std::map<size_t, AbstractPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    pl.includeParameters(it->second->getBranchLengthParameters());

  return pl;
}
      
ParameterList SumOfPhyloLikelihood::getSubstitutionModelParameters() const
{
  ParameterList pl;
  for (std::map<size_t, AbstractPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    pl.includeParameters(it->second->getSubstitutionModelParameters());

  return pl;
}

ParameterList SumOfPhyloLikelihood::getRateDistributionParameters() const
{
  ParameterList pl;
  for (std::map<size_t, AbstractPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    pl.includeParameters(it->second->getRateDistributionParameters());

  return pl;
}

ParameterList SumOfPhyloLikelihood::getRootFrequenciesParameters() const
{
  ParameterList pl;
  for (std::map<size_t, AbstractPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    pl.includeParameters(it->second->getRootFrequenciesParameters());

  return pl;
}

ParameterList SumOfPhyloLikelihood::getDerivableParameters() const
{
  // patch, to be fixed properly later
  return ParameterList();

  ParameterList pl;
  // for (std::map<size_t, AbstractPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
  //   pl.includeParameters(it->second->getDerivableParameters());
  
  return pl;
}

ParameterList SumOfPhyloLikelihood::getNonDerivableParameters() const
{
  // patch, to be fixed properly later
  return getParameters();

  ParameterList pl;
  for (std::map<size_t, AbstractPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    pl.includeParameters(it->second->getNonDerivableParameters());
  
  return pl;
}

void SumOfPhyloLikelihood::enableDerivatives(bool yn)
{
  AbstractPhyloLikelihood::enableDerivatives(yn);
  
  for (std::map<size_t, AbstractPhyloLikelihood*>::iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    it->second->enableDerivatives(yn);
}

void SumOfPhyloLikelihood::enableFirstOrderDerivatives(bool yn)
{
  AbstractPhyloLikelihood::enableFirstOrderDerivatives(yn);
  
  for (std::map<size_t, AbstractPhyloLikelihood*>::iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    it->second->enableFirstOrderDerivatives(yn);
}

void SumOfPhyloLikelihood::enableSecondOrderDerivatives(bool yn)
{
  AbstractPhyloLikelihood::enableSecondOrderDerivatives(yn);

  for (std::map<size_t, AbstractPhyloLikelihood*>::iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    it->second->enableSecondOrderDerivatives(yn);
}

