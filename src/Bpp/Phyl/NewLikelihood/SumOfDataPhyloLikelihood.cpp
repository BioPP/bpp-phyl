//
// File: SumOfDataPhyloLikelihood.cpp
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

#include "SumOfDataPhyloLikelihood.h"

using namespace bpp;
using namespace newlik;
using namespace std;

SumOfDataPhyloLikelihood::SumOfDataPhyloLikelihood() :
  AbstractParametrizable(""),  
  mSDP_(),
  computeFirstOrderDerivatives_(true),
  computeSecondOrderDerivatives_(true),
  minusLogLik_(0)
{
}

SumOfDataPhyloLikelihood::SumOfDataPhyloLikelihood(std::map<size_t, SingleDataPhyloLikelihood*>& mSDP) :
  AbstractParametrizable(""),  
  mSDP_(),
  computeFirstOrderDerivatives_(true),
  computeSecondOrderDerivatives_(true),
  minusLogLik_(0)
{
  for (std::map<size_t, SingleDataPhyloLikelihood*>::const_iterator it=mSDP.begin(); it != mSDP.end(); it++)
  {
    addSingleDataPhylolikelihood(it->first, it->second);
    includeParameters_(mSDP_[it->first]->getParameters());
  }
  
}

SumOfDataPhyloLikelihood::~SumOfDataPhyloLikelihood()
{
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    delete it->second;
}

SumOfDataPhyloLikelihood::SumOfDataPhyloLikelihood(const SumOfDataPhyloLikelihood& sd) :
  AbstractParametrizable(sd),
  mSDP_(),
  computeFirstOrderDerivatives_(sd.computeFirstOrderDerivatives_),
  computeSecondOrderDerivatives_(sd.computeSecondOrderDerivatives_),
  minusLogLik_(sd.minusLogLik_)
{
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=sd.mSDP_.begin(); it != sd.mSDP_.end(); it++)
    addSingleDataPhylolikelihood(it->first, it->second->clone());
}

SumOfDataPhyloLikelihood& SumOfDataPhyloLikelihood::operator=(const SumOfDataPhyloLikelihood& sd)
{
  AbstractParametrizable::operator=(sd);

  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    delete it->second;
  
  computeFirstOrderDerivatives_=sd.computeFirstOrderDerivatives_;
  computeSecondOrderDerivatives_=sd.computeSecondOrderDerivatives_;

  minusLogLik_=sd.minusLogLik_;
  
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=sd.mSDP_.begin(); it != sd.mSDP_.end(); it++)
    addSingleDataPhylolikelihood(it->first, it->second->clone());
  
  return *this;
}


void SumOfDataPhyloLikelihood::addSingleDataPhylolikelihood(size_t nPhyl, SingleDataPhyloLikelihood* SDP)
{
  if (mSDP_.find(nPhyl)!=mSDP_.end())
    throw Exception("SumOfDataPhyloLikelihood::addSingleDataPhylolikelihood: map number already used : " + TextTools::toString(nPhyl));
  
  if (dynamic_cast<AbstractSingleDataPhyloLikelihood*>(SDP)!=NULL){
    mSDP_[nPhyl]=dynamic_cast<AbstractSingleDataPhyloLikelihood*>(SDP);
    includeParameters_(SDP->getParameters());
  } 
}


void SumOfDataPhyloLikelihood::setParameters(const ParameterList& parameters) throw (ParameterNotFoundException, ConstraintException)
{
  AbstractParametrizable::setParametersValues(parameters);
  
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    it->second->setParameters(parameters);
}

double SumOfDataPhyloLikelihood::getFirstOrderDerivative(const std::string& variable) const throw (Exception)
{
  double x=0;
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    if (it->second->hasParameter(variable))
      x+=it->second->getFirstOrderDerivative(variable);
  
  return x;
}

double SumOfDataPhyloLikelihood::getSecondOrderDerivative(const std::string& variable) const throw (Exception)
{
  double x=0;
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    if (it->second->hasParameter(variable))
      x+=it->second->getSecondOrderDerivative(variable);
  
  return x;
}

std::vector<size_t> SumOfDataPhyloLikelihood::getNumbersOfSingleDataPhyloLikelihoods() const
{
  vector<size_t> vNum;
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    vNum.push_back(it->first);

  return vNum;
}


ParameterList SumOfDataPhyloLikelihood::getBranchLengthParameters() const
{
  ParameterList pl;
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    pl.includeParameters(it->second->getBranchLengthParameters());

  return pl;
}
      
ParameterList SumOfDataPhyloLikelihood::getSubstitutionModelParameters() const
{
  ParameterList pl;
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    pl.includeParameters(it->second->getSubstitutionModelParameters());

  return pl;
}

ParameterList SumOfDataPhyloLikelihood::getRateDistributionParameters() const
{
  ParameterList pl;
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    pl.includeParameters(it->second->getRateDistributionParameters());

  return pl;
}

ParameterList SumOfDataPhyloLikelihood::getRootFrequenciesParameters() const
{
  ParameterList pl;
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    pl.includeParameters(it->second->getRootFrequenciesParameters());

  return pl;
}

ParameterList SumOfDataPhyloLikelihood::getDerivableParameters() const
{
  ParameterList pl;
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    pl.includeParameters(it->second->getDerivableParameters());
  
  return pl;
}

ParameterList SumOfDataPhyloLikelihood::getNonDerivableParameters() const
{
  ParameterList pl;
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::const_iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    pl.includeParameters(it->second->getNonDerivableParameters());
  
  return pl;
}

void SumOfDataPhyloLikelihood::enableDerivatives(bool yn)
{
  for (std::map<size_t, AbstractSingleDataPhyloLikelihood*>::iterator it=mSDP_.begin(); it != mSDP_.end(); it++)
    it->second->enableDerivatives(yn);
}

