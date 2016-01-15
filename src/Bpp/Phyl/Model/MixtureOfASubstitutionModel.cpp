//
// File: MixtureOfASubstitutionModel.cpp
// Created by: David Fournier, Laurent Gueguen
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

#include "MixtureOfASubstitutionModel.h"

#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Exceptions.h>

#include <string>

using namespace bpp;
using namespace std;


MixtureOfASubstitutionModel::MixtureOfASubstitutionModel(
  const Alphabet* alpha,
  SubstitutionModel* model,
  std::map<std::string, DiscreteDistribution*> parametersDistributionsList,
  int ffrom,
  int tto) throw (Exception) :
  AbstractParameterAliasable(model->getNamespace()),
  AbstractMixedSubstitutionModel(alpha, model->getStateMap().clone(), model->getNamespace()),
  distributionMap_(),
  from_(ffrom),
  to_(tto)
{
  if (to_ >= int(alpha->getSize()))
    throw BadIntegerException("Bad state in alphabet", to_);
  if (from_ >= int(alpha->getSize()))
    throw BadIntegerException("Bad state in alphabet", from_);

  size_t c, i;
  string s1, s2, t;
  map<string, DiscreteDistribution*>::iterator it;

  // Initialization of distributionMap_.

  vector<string> parnames = model->getParameters().getParameterNames();

  for (i = 0; i < model->getNumberOfParameters(); i++)
  {
    s1 = parnames[i];
    s2 = model->getParameterNameWithoutNamespace(s1);

    if (parametersDistributionsList.find(s2) != parametersDistributionsList.end())
      distributionMap_[s1] = parametersDistributionsList.find(s2)->second->clone();
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
//    modelsContainer_[i]->addRateParameter();
    modelsContainer_[i]->setNamespace(model->getNamespace());

    vProbas_.push_back(1.0 / static_cast<double>(c));
    vRates_.push_back(1.0);
  }

  // Initialization of parameters_.

  Parameter pm;
  DiscreteDistribution* pd;

  for (it = distributionMap_.begin(); it != distributionMap_.end(); it++)
  {
    pm = model->getParameter(model->getParameterNameWithoutNamespace(getParameterNameWithoutNamespace(it->first)));
    pd = it->second;

    if (pm.hasConstraint())
      pd->restrictToConstraint(*pm.getConstraint());

    if (!dynamic_cast<ConstantDistribution*>(it->second))
    {
      const ParameterList pl = pd->getParameters();
      for (i = 0; i != it->second->getNumberOfParameters(); i++)
      {
        addParameter_(pl[i].clone());
      }
    }
    else
      addParameter_(new Parameter(it->first, pd->getCategory(0), (pd->getParameter("value").getConstraint()) ? pd->getParameter("value").getConstraint()->clone() : 0, true));
  }
  updateMatrices();
}

MixtureOfASubstitutionModel::MixtureOfASubstitutionModel(const MixtureOfASubstitutionModel& msm) :
  AbstractParameterAliasable(msm),
  AbstractMixedSubstitutionModel(msm),
  distributionMap_(),
  from_(msm.from_),
  to_(msm.to_)
{
  map<string, DiscreteDistribution*>::const_iterator it;

  for (it = msm.distributionMap_.begin(); it != msm.distributionMap_.end(); it++)
  {
    distributionMap_[it->first] = it->second->clone();
  }
}

MixtureOfASubstitutionModel& MixtureOfASubstitutionModel::operator=(const MixtureOfASubstitutionModel& msm)
{
  AbstractParameterAliasable::operator=(msm);
  AbstractMixedSubstitutionModel::operator=(msm);
  from_ = msm.from_;
  to_ = msm.to_;

  // Clear existing containers:
  distributionMap_.clear();

  // Now copy new containers:
  map<string, DiscreteDistribution*>::const_iterator it;
  for (it = msm.distributionMap_.begin(); it != msm.distributionMap_.end(); it++)
  {
    distributionMap_[it->first] = it->second->clone();
  }
  return *this;
}


MixtureOfASubstitutionModel::~MixtureOfASubstitutionModel()
{
  map<string, DiscreteDistribution*>::iterator it;

  for (it = distributionMap_.begin(); it != distributionMap_.end(); it++)
  {
    delete it->second;
  }
}

const DiscreteDistribution* MixtureOfASubstitutionModel::getDistribution(std::string& parName) const
{
  if (distributionMap_.find(parName) != distributionMap_.end())
    return distributionMap_.find(parName)->second;
  else
    return NULL;
}

void MixtureOfASubstitutionModel::updateMatrices()
{
  string s, t;
  size_t i, j, l;
  double d;
  ParameterList pl;
  map<string, DiscreteDistribution*>::iterator it;

  // Update of distribution parameters from the parameters_ member
  // data. (reverse operation compared to what has been done in the
  // constructor).
  //  vector<string> v=getParameters().getParameterNames();

  for (it = distributionMap_.begin(); it != distributionMap_.end(); it++)
  {
    if (dynamic_cast<ConstantDistribution*>(it->second) == NULL)
    {
      vector<string> vDistnames = it->second->getParameters().getParameterNames();
      for (i = 0; i < it->second->getNumberOfParameters(); i++)
      {
        d = getParameterValue(getParameterNameWithoutNamespace(vDistnames[i]));
        pl.addParameter(Parameter(vDistnames[i], d));
      }
      it->second->matchParametersValues(pl);
      pl.reset();
    }
    else
    {
      t = it->second->getNamespace();
      d = getParameter(getParameterNameWithoutNamespace(t.substr(0, t.length() - 1))).getValue();
      it->second->setParameterValue("value", d);
    }
  }

  for (i = 0; i < modelsContainer_.size(); i++)
  {
    vProbas_[i] = 1;
    j = i;
    for (it = distributionMap_.begin(); it != distributionMap_.end(); it++)
    {
      s = it->first;
      l = j % it->second->getNumberOfCategories();

      d = it->second->getCategory(l);
      vProbas_[i] *= it->second->getProbability(l);
      if (pl.hasParameter(s))
        pl.setParameterValue(s, d);
      else
        pl.addParameter(Parameter(s, d));

      j = j / it->second->getNumberOfCategories();
    }

    modelsContainer_[i]->matchParametersValues(pl);
  }

  //  setting the equilibrium freqs
  for (i = 0; i < getNumberOfStates(); i++)
  {
    freq_[i] = 0;
    for (j = 0; j < modelsContainer_.size(); j++)
    {
      freq_[i] += vProbas_[j] * modelsContainer_[j]->freq(i);
    }
  }

  // setting the rates, if to_ & from_ are different from -1

  if (to_ >= 0 && from_ >= 0)
  {
    Vdouble vd;

    for (j = 0; j < modelsContainer_.size(); j++)
    {
      vd.push_back(1 / getNModel(j)->Qij(static_cast<size_t>(from_), static_cast<size_t>(to_)));
    }

    setVRates(vd);
  }
}

void MixtureOfASubstitutionModel::setFreq(std::map<int, double>& m)
{
  modelsContainer_[0]->setFreq(m);
  matchParametersValues(modelsContainer_[0]->getParameters());
}

const SubstitutionModel* MixtureOfASubstitutionModel::getSubModelWithName(const std::string& name) const
{
  size_t nbmod=getNumberOfModels();
  
  for (size_t i=0; i<nbmod; i++)
    if (getNModel(i)->getName()==name)
      return getNModel(i);

  return NULL;
}

Vint MixtureOfASubstitutionModel::getSubmodelNumbers(string& desc) const
{
  vector<string> parnames = modelsContainer_[0]->getParameters().getParameterNames();
  std::map<std::string, size_t> msubn;
  map<string, DiscreteDistribution*>::const_iterator it;

  StringTokenizer st(desc, ",");
  while (st.hasMoreToken())
  {
    string param = st.nextToken();
    string::size_type index = param.rfind("_");
    if (index == string::npos)
      throw Exception("MixtureOfASubstitutionModel::getSubmodelNumbers parameter description should contain a number" + param);
    msubn[param.substr(0, index)] = TextTools::to<size_t>(param.substr(index + 1, 4)) - 1;
  }

  Vint submodnb;
  size_t i, j, l;
  string s;

  bool nameok = false;

  for (i = 0; i < modelsContainer_.size(); i++)
  {
    j = i;
    for (it = distributionMap_.begin(); it != distributionMap_.end(); it++)
    {
      s = it->first;
      l = j % it->second->getNumberOfCategories();

      if (msubn.find(s) != msubn.end())
      {
        nameok = true;
        if (msubn[s] != l)
          break;
      }

      j = j / it->second->getNumberOfCategories();
    }
    if (nameok && it == distributionMap_.end())
      submodnb.push_back(static_cast<int>(i));
  }

  return submodnb;
}
