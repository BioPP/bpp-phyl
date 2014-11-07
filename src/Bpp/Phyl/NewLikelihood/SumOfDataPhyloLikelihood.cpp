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
  vSDP_(),
  numberOfSDP_(0),
  computeFirstOrderDerivatives_(true),
  computeSecondOrderDerivatives_(true),
  minusLogLik_(0)
{
}

SumOfDataPhyloLikelihood::SumOfDataPhyloLikelihood(std::vector<SingleDataPhyloLikelihood*> vSDP) :
  AbstractParametrizable(""),  
  vSDP_(),
  numberOfSDP_(0),
  computeFirstOrderDerivatives_(true),
  computeSecondOrderDerivatives_(true),
  minusLogLik_(0)
{
  for (size_t i = 0; i < vSDP.size(); i++)
  {
    addSingleDataPhylolikelihood(vSDP[i]);
    includeParameters_(vSDP[i]->getParameters());
  }
  
}

SumOfDataPhyloLikelihood::~SumOfDataPhyloLikelihood()
{
  for (size_t i = 0; i < numberOfSDP_; i++)
    delete vSDP_[i];
}

SumOfDataPhyloLikelihood::SumOfDataPhyloLikelihood(const SumOfDataPhyloLikelihood& sd) :
  AbstractParametrizable(sd),
  vSDP_(),
  numberOfSDP_(0),
  computeFirstOrderDerivatives_(sd.computeFirstOrderDerivatives_),
  computeSecondOrderDerivatives_(sd.computeSecondOrderDerivatives_),
  minusLogLik_(sd.minusLogLik_)
{
  for (size_t i = 0; i < sd.numberOfSDP_; i++)
    addSingleDataPhylolikelihood(sd.vSDP_[i]->clone());
}

SumOfDataPhyloLikelihood& SumOfDataPhyloLikelihood::operator=(const SumOfDataPhyloLikelihood& sd)
{
  AbstractParametrizable::operator=(sd);

  for (size_t i = 0; i < numberOfSDP_; i++)
    delete vSDP_[i];
  
  numberOfSDP_=0;
  computeFirstOrderDerivatives_=sd.computeFirstOrderDerivatives_;
  computeSecondOrderDerivatives_=sd.computeSecondOrderDerivatives_;

  minusLogLik_=sd.minusLogLik_;
  
  for (size_t i = 0; i < sd.numberOfSDP_; i++)
    addSingleDataPhylolikelihood(sd.vSDP_[i]->clone());
  
  return *this;
}


void SumOfDataPhyloLikelihood::addSingleDataPhylolikelihood(SingleDataPhyloLikelihood* SDP)
{
  if (dynamic_cast<AbstractSingleDataPhyloLikelihood*>(SDP)!=NULL){
    vSDP_.push_back(dynamic_cast<AbstractSingleDataPhyloLikelihood*>(SDP));
    numberOfSDP_++;
    includeParameters_(SDP->getParameters());
  }
}


void SumOfDataPhyloLikelihood::setParameters(const ParameterList& parameters) throw (ParameterNotFoundException, ConstraintException)
{
  AbstractParametrizable::setParametersValues(parameters);
  
  for (size_t i = 0; i < numberOfSDP_; i++)
    vSDP_[i]->setParameters(parameters);
}

double SumOfDataPhyloLikelihood::getFirstOrderDerivative(const std::string& variable) const throw (Exception)
{
  double x=0;
  for (size_t i = 0; i < numberOfSDP_; i++) 
    if (vSDP_[i]->hasParameter(variable))
      x+=vSDP_[i]->getFirstOrderDerivative(variable);
  
  return x;
}

double SumOfDataPhyloLikelihood::getSecondOrderDerivative(const std::string& variable) const throw (Exception)
{
  double x=0;
  for (size_t i = 0; i < numberOfSDP_; i++)
    if (vSDP_[i]->hasParameter(variable))
      x+=vSDP_[i]->getSecondOrderDerivative(variable);

  return x;
}

ParameterList SumOfDataPhyloLikelihood::getBranchLengthParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < numberOfSDP_; i++)
    pl.includeParameters(vSDP_[i]->getBranchLengthParameters());

  return pl;
}
      
ParameterList SumOfDataPhyloLikelihood::getSubstitutionModelParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < numberOfSDP_; i++)
    pl.includeParameters(vSDP_[i]->getSubstitutionModelParameters());

  return pl;
}

ParameterList SumOfDataPhyloLikelihood::getRateDistributionParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < numberOfSDP_; i++)
    pl.includeParameters(vSDP_[i]->getRateDistributionParameters());

  return pl;
}

ParameterList SumOfDataPhyloLikelihood::getRootFrequenciesParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < numberOfSDP_; i++)
    pl.includeParameters(vSDP_[i]->getRootFrequenciesParameters());

  return pl;
}

ParameterList SumOfDataPhyloLikelihood::getDerivableParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < numberOfSDP_; i++)
    pl.includeParameters(vSDP_[i]->getDerivableParameters());
  
  return pl;
}

ParameterList SumOfDataPhyloLikelihood::getNonDerivableParameters() const
{
  ParameterList pl;
  for (size_t i = 0; i < numberOfSDP_; i++)
    pl.includeParameters(vSDP_[i]->getNonDerivableParameters());
  
  return pl;
}

void SumOfDataPhyloLikelihood::enableDerivatives(bool yn)
{
  for (size_t i = 0; i < numberOfSDP_; i++)
    vSDP_[i]->enableDerivatives(yn);
}

