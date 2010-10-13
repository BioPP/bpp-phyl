//
// File: MixtureOfSubstitutionModels.cpp
// Created by: Laurent Gueguen
// Date: mardi 14 septembre 2010, à 20h 43
//
/*
  Copyright or © or Copr. CNRS, (November 16, 2004)

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

#include "MixtureOfSubstitutionModels.h"

#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Text/TextTools.h>

#include <string>

using namespace bpp;
using namespace std;

MixtureOfSubstitutionModels::MixtureOfSubstitutionModels(const Alphabet* alpha,
                                                         vector<SubstitutionModel*> vpModel_) :
  MixedSubstitutionModel(alpha, "Mixture."),
  modelsContainer_(),
  Vprobas_(),
  Vrates_()
{
  unsigned int i, nbmod = vpModel_.size();

  for (i = 0; i < nbmod; i++){
    if (vpModel_[i]==NULL)
      throw Exception("Empty model number " + TextTools::toString(i) + " in MixtureOfSubstitutionModels constructor");
    for (unsigned int j=i+1;j<nbmod;j++)
      if (vpModel_[i]==vpModel_[j])
        throw Exception("Same model at positions " + TextTools::toString(i) + " and " +
                        TextTools::toString(j) + " in MixtureOfSubstitutionModels constructor");
  }
  
  // Initialization of modelsContainer_.

  for (i = 0; i < nbmod; i++)
    {
      modelsContainer_.push_back(vpModel_[i]);
      Vprobas_.push_back(1.0/nbmod);
      Vrates_.push_back(1.0);
    }

  // Initialization of parameters_.

  // relative rates and probas
  for (i = 0; i < nbmod - 1; i++)
    {
      addParameter_(Parameter("Mixture.relproba" + TextTools::toString(i+1), 1.0 / (nbmod - i ), &Parameter::PROP_CONSTRAINT_EX));
      addParameter_(Parameter("Mixture.relrate" + TextTools::toString(i+1), 1.0 / (nbmod - i), &Parameter::PROP_CONSTRAINT_EX));
    }

  // models parameters

  for (i = 0; i < nbmod; i++)
    {
      modelsContainer_[i]->setNamespace("Mixture." + TextTools::toString(i+1) + "_" + vpModel_[i]->getNamespace());
      addParameters_(vpModel_[i]->getParameters());
    }

  updateMatrices();
}

MixtureOfSubstitutionModels::MixtureOfSubstitutionModels(const Alphabet* alpha,
                                                         vector<SubstitutionModel*> vpModel_,
                                                         Vdouble& vproba,
                                                         Vdouble& vrate) :
  MixedSubstitutionModel(alpha, "Mixture."),
  modelsContainer_(),
  Vprobas_(vproba),
  Vrates_(vrate)
{
  unsigned int i, nbmod = vpModel_.size();

  for (i = 0; i < nbmod; i++){
    if (vpModel_[i]==NULL)
      throw Exception("Empty model number " + TextTools::toString(i) + " in MixtureOfSubstitutionModels constructor");
    for (unsigned int j=i+1;j<nbmod;j++)
      if (vpModel_[i]==vpModel_[j])
        throw Exception("Same model at positions " + TextTools::toString(i) + " and " +
                        TextTools::toString(j) + " in MixtureOfSubstitutionModels constructor");
  }

  double x=0;
  double y=0;
  
  for (i=0; i < nbmod; i++){
    if (vrate[i]<=0)
      throw Exception("Non positive rate: " + TextTools::toString(vrate[i]) + " in MixtureOfSubstitutionModels constructor.");
    if (vproba[i]<=0)
      throw Exception("Non positive probability: " + TextTools::toString(vproba[i]) + " in MixtureOfSubstitutionModels constructor.");
    x+=vproba[i];
    y+=vproba[i]*vrate[i];
  }

  if (fabs(1. - x) > NumConstants::SMALL)
    throw Exception("Probabilities must equal 1 (sum = " + TextTools::toString(x) + ").");
  if (fabs(1. - y) > NumConstants::SMALL)
    throw Exception("Expectation on rates must equal 1 (E =" + TextTools::toString(y) + ").");
  
  
  // Initialization of modelsContainer_.

  for (i = 0; i < nbmod; i++) 
      modelsContainer_.push_back(vpModel_[i]);

  // Initialization of parameters_.

  
  // relative rates and probas
  x=0;y=0;
  
  for (i = 0; i < nbmod - 1; i++)
    {
      addParameter_(Parameter("Mixture.relproba" + TextTools::toString(i+1), vproba[i] / (1 - x), &Parameter::PROP_CONSTRAINT_EX));
      x+=vproba[i];
      addParameter_(Parameter("Mixture.relrate" + TextTools::toString(i+1), vproba[i] * vrate[i] / (1- y), &Parameter::PROP_CONSTRAINT_EX));
      y+=vproba[i]*vrate[i];
    }

  // models parameters

  for (i = 0; i < nbmod; i++)
    {
      modelsContainer_[i]->setNamespace("Mixture." + TextTools::toString(i+1) + "_" + vpModel_[i]->getNamespace());
      addParameters_(vpModel_[i]->getParameters());
    }

  updateMatrices();
}

MixtureOfSubstitutionModels::MixtureOfSubstitutionModels(const MixtureOfSubstitutionModels& msm) :
  MixedSubstitutionModel(msm),
  modelsContainer_(),
  Vprobas_(),
  Vrates_()
{
  for (unsigned int i = 0; i < msm.modelsContainer_.size(); i++)
    {
      modelsContainer_.push_back(msm.modelsContainer_[i]->clone());
      Vprobas_.push_back(msm.Vprobas_[i]);
      Vrates_.push_back(msm.Vrates_[i]);
    }
}

MixtureOfSubstitutionModels& MixtureOfSubstitutionModels::operator=(const MixtureOfSubstitutionModels& msm)
{
  MixedSubstitutionModel::operator=(msm);
  
  //Clear existing containers:
  modelsContainer_.clear();
  Vprobas_.clear();
  Vrates_.clear();
  
  //Now copy new containers:

  for (unsigned int i = 0; i < msm.modelsContainer_.size(); i++)
    {
      modelsContainer_.push_back(msm.modelsContainer_[i]->clone());
      Vprobas_.push_back(msm.Vprobas_[i]);
      Vrates_.push_back(msm.Vrates_[i]);
    }
  return *this;
}


MixtureOfSubstitutionModels::~MixtureOfSubstitutionModels()
{
  for (unsigned int i = 0; i < modelsContainer_.size(); i++)
    {
      delete modelsContainer_[i];
    }
}



void MixtureOfSubstitutionModels::updateMatrices()
{
  unsigned int i, j, nbmod = modelsContainer_.size();
  
  double x,y;
  x = 1.0;

  for (i = 0; i < nbmod-1; i++){
    y =getParameterValue("relproba" + TextTools::toString(i+1));
    Vprobas_[i] = x*y;
    x *= 1 - y;      
  }
  Vprobas_[nbmod-1]=x;

  x = 1.0;
  for (i = 0; i < nbmod-1; i++){
    y =getParameterValue("relrate" + TextTools::toString(i+1));
    Vrates_[i] = x*y/Vprobas_[i];
    x *= 1 - y;      
  }
  Vrates_[nbmod-1]=x/Vprobas_[nbmod-1];

  /// models
  
  for ( i = 0; i < nbmod; i++){
    modelsContainer_[i]->matchParametersValues(getParameters());
    modelsContainer_[i]->setRate(Vrates_[i]);
  }

  /// freq_
  
  for (i = 0; i < getNumberOfStates(); i++){
    freq_[i] = 0;
    for (j = 0; j < modelsContainer_.size(); j++)
      freq_[i] += Vprobas_[i]*modelsContainer_[j]->freq(i);
  }
  
}

unsigned int MixtureOfSubstitutionModels::getNumberOfStates() const
{
  return modelsContainer_[0]->getNumberOfStates();
}


double MixtureOfSubstitutionModels::Pij_t(unsigned int i, unsigned int j, double t) const
{
  double x=0;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
      x+= modelsContainer_[n]->Pij_t(i,j,t)*Vprobas_[n];

  return x;
}


double MixtureOfSubstitutionModels::dPij_dt(unsigned int i, unsigned int j, double t) const
{
  double x=0;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
    x+= modelsContainer_[n]->dPij_dt(i,j,t)*Vprobas_[n];

  return x;
}


double MixtureOfSubstitutionModels::d2Pij_dt2(unsigned int i, unsigned int j, double t) const
{
  double x=0;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
    x+= modelsContainer_[n]->d2Pij_dt2(i,j,t)*Vprobas_[n];

  return x;
}


const Matrix<double>& MixtureOfSubstitutionModels::getPij_t(double t) const
{
  for (unsigned int i=0; i< getNumberOfStates(); i++)
    for (unsigned int j=0; j< getNumberOfStates(); j++)
      pijt_(i,j)=Pij_t(i,j,t);

  return pijt_;
}


const Matrix<double>& MixtureOfSubstitutionModels::getdPij_dt(double t) const
{
  for (unsigned int i=0; i< getNumberOfStates(); i++)
    for (unsigned int j=0; j< getNumberOfStates(); j++)
      dpijt_(i,j)=dPij_dt(i,j,t);

  return dpijt_;
}


const Matrix<double>& MixtureOfSubstitutionModels::getd2Pij_dt2(double t) const
{
  for (unsigned int i=0; i< getNumberOfStates(); i++)
    for (unsigned int j=0; j< getNumberOfStates(); j++)
      dpijt_(i,j)=d2Pij_dt2(i,j,t);

  return d2pijt_;
}


const Vdouble& MixtureOfSubstitutionModels::getFrequencies()
{
  for (unsigned int i=0; i< getNumberOfStates(); i++)
    freq_[i]=freq(i);
  return freq_;
}


double MixtureOfSubstitutionModels::freq(unsigned int i) const
{
  double x=0;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
    x+= modelsContainer_[n]->freq(i)*Vprobas_[n];

  return x;
}



