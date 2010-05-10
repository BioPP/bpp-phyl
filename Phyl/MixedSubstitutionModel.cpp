//
// File: MixedSubstitutionModel.cpp
// Created by: David Fournier, Laurent Gueguen
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

#include "MixedSubstitutionModel.h"
#include <string>

using namespace bpp;
using namespace std;


MixedSubstitutionModel::MixedSubstitutionModel(
    const Alphabet* alpha,
    SubstitutionModel* model,
    std::map<std::string, DiscreteDistribution*> parametersDistributionsList) :
  AbstractSubstitutionModel(alpha, ""),
  distributionMap_(),
  modelsContainer_(),
  probas_()
{
  unsigned int c, i;
  double d;
  string s1, s2, t;
  map<string, DiscreteDistribution*>::iterator it;

  // Initialization of distributionMap_.

  vector<string> parnames = model->getParameters().getParameterNames();

  for (i = 0; i < model->getNumberOfParameters(); i++)
  {
    s1 = parnames[i];
    s2 = model->getParameterNameWithoutNamespace(s1);
   
    if (parametersDistributionsList.find(s2) != parametersDistributionsList.end())
      distributionMap_[s1] = dynamic_cast<DiscreteDistribution*>(parametersDistributionsList.find(s2)->second->clone());
    else
      distributionMap_[ s1] = new ConstantDistribution(model->getParameterValue(s2));


    if (dynamic_cast<ConstantDistribution*>(distributionMap_[s1]) == 0)
      distributionMap_[s1]->setNamespace(s1 + "_" + distributionMap_[s1]->getNamespace());
    else
      distributionMap_[s1]->setNamespace(s1 + "_");
  }


  // Initialization of modelsContainer_.

  c = 1;

  for (it = distributionMap_.begin(); it != distributionMap_.end(); it++)
  {
    c *= it->second->getNumberOfCategories();
  }

  for (i = 0; i < c; i++)
  {
    modelsContainer_.push_back(model->clone());
    modelsContainer_[i]->setNamespace(model->getNamespace());
    probas_.push_back(1.0/c);
  }

  // Initialization of parameters_.

  Constraint* pc;

  for (it = distributionMap_.begin(); it != distributionMap_.end(); it++)
  {
    if (dynamic_cast<ConstantDistribution*>(it->second) == NULL)
    {
      for (i = 0; i != it->second->getNumberOfParameters(); i++)
      {
        t = it->second->getParameters().getParameterNames()[i];
        d = it->second->getParameters().getParameter(t).getValue();
        if (it->second->getParameters().getParameter(t).hasConstraint()){
          pc = it->second->getParameters().getParameter(t).getConstraint()->clone();
          addParameter_(Parameter(t,d,pc,true));
        }
        else
          addParameter_(Parameter(t,d));
      }
    }
    else
    {
      t = it->first;
      pc = model->getParameter(model->getParameterNameWithoutNamespace(getParameterNameWithoutNamespace(t))).getConstraint()->clone();
      d = it->second->getCategory(0);
      addParameter_(Parameter(t,d,pc,true));
    }
  }

  updateMatrices();
}

MixedSubstitutionModel::MixedSubstitutionModel(const MixedSubstitutionModel& msm) :
  AbstractSubstitutionModel(msm),
  distributionMap_(),
  modelsContainer_(),
  probas_()
{
  map<string, DiscreteDistribution*>::const_iterator it;

  for (it = msm.distributionMap_.begin(); it != msm.distributionMap_.end(); it++)
  {
    distributionMap_[it->first] = dynamic_cast<DiscreteDistribution*>(it->second->clone());
  }

  for (unsigned int i = 0; i < msm.modelsContainer_.size(); i++)
  {
    modelsContainer_.push_back(msm.modelsContainer_[i]->clone());
    probas_.push_back(msm.probas_[i]);
  }
}

MixedSubstitutionModel& MixedSubstitutionModel::operator=(const MixedSubstitutionModel& msm)
{
  AbstractSubstitutionModel::operator=(msm);
  
  //Clear existing containers:
  distributionMap_.clear();
  modelsContainer_.clear();
  probas_.clear();
  
  //Now copy new containers:
  map<string, DiscreteDistribution*>::const_iterator it;
  for (it = msm.distributionMap_.begin(); it != msm.distributionMap_.end(); it++)
  {
    distributionMap_[it->first] = dynamic_cast<DiscreteDistribution*>(it->second->clone());
  }

  for (unsigned int i = 0; i < msm.modelsContainer_.size(); i++)
  {
    modelsContainer_.push_back(msm.modelsContainer_[i]->clone());
    probas_.push_back(msm.probas_[i]);
  }
  return *this;
}


MixedSubstitutionModel::~MixedSubstitutionModel()
{
   unsigned int i;
   map<string, DiscreteDistribution*>::iterator it;

  for (it = distributionMap_.begin(); it != distributionMap_.end(); it++)
  {
    delete it->second;
  }
  for (i = 0; i < modelsContainer_.size(); i++)
  {
    delete modelsContainer_[i];
  }
}

void MixedSubstitutionModel::updateMatrices()
{
   string s,t;
   unsigned int i, j, l;
   double d;
   ParameterList pl;
   map<string, DiscreteDistribution*>::iterator it;

  // Update of distribution parameters from the parameters_ member
  // data. (reverse operation compared to what has been done in the
  // constructor).

  for (it = distributionMap_.begin(); it != distributionMap_.end(); it++)
  {
    if (dynamic_cast<ConstantDistribution*>(it->second) == NULL)
    {
      for (i = 0; i < it->second->getNumberOfParameters(); i++)
      {
        t = it->second->getParameters().getParameterNames()[i];
        d = getParameter(getParameterNameWithoutNamespace(t)).getValue();
        it->second->setParameterValue(it->second->getParameterNameWithoutNamespace(t),d);
      }
    }
    else
    {
      t = it->second->getNamespace();
      d = getParameter(getParameterNameWithoutNamespace(t.substr(0,t.length() - 1))).getValue();
      it->second->setParameterValue("value",d);
    }
  }


  for (i = 0; i < modelsContainer_.size(); i++)
  {
    probas_[i]=1;
    j = i;
    for (it = distributionMap_.begin(); it != distributionMap_.end(); it++)
    {
      s = it->first;
      l = j % it->second->getNumberOfCategories();


      d = distributionMap_.find(s)->second->getCategory(l);
      probas_[i]*=distributionMap_.find(s)->second->getProbability(l);
      if (pl.hasParameter(s))
        pl.setParameterValue(s,d);
      else
        pl.addParameter(Parameter(s,d));

      j = j / it->second->getNumberOfCategories();
    }

    modelsContainer_[i]->matchParametersValues(pl);
  }


  for (i = 0; i < getNumberOfStates(); i++)
  {
    freq_[i] = 0;
    for (j = 0; j < modelsContainer_.size(); j++)
    {
      freq_[i] += modelsContainer_[j]->freq(i);
    }
    freq_[i] /= modelsContainer_.size();
  }
}


void MixedSubstitutionModel::setFreq(std::map<int,double>& m)
{
  throw Exception("setFreq method is not available for MixedSubstitutionModel.");
}

unsigned int MixedSubstitutionModel::getNumberOfStates() const
{
  return modelsContainer_[0]->getNumberOfStates();
}

