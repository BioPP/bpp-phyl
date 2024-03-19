// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/VectorTools.h>
#include <string>

#include "MixtureOfATransitionModel.h"

using namespace bpp;
using namespace std;


MixtureOfATransitionModel::MixtureOfATransitionModel(
    shared_ptr<const Alphabet> alpha,
    unique_ptr<TransitionModelInterface> model,
    std::map<std::string, unique_ptr<DiscreteDistributionInterface>>& parametersDistributionsList,
    int ffrom,
    int tto) :
  AbstractParameterAliasable(model->getNamespace()),
  AbstractTransitionModel(alpha, model->getStateMap(), model->getNamespace()),
  AbstractMixedTransitionModel(alpha, getStateMap(), model->getNamespace()),
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

  // Initialization of distributionMap_.

  vector<string> parnames = model->getParameters().getParameterNames();

  for (i = 0; i < model->getNumberOfParameters(); i++)
  {
    s1 = parnames[i];
    s2 = model->getParameterNameWithoutNamespace(s1);

    if (parametersDistributionsList.find(s2) != parametersDistributionsList.end())
      distributionMap_[s1].reset(parametersDistributionsList.find(s2)->second->clone());
    else
      distributionMap_[s1].reset(new ConstantDistribution(model->getParameterValue(s2)));


    if (dynamic_cast<ConstantDistribution*>(distributionMap_[s1].get()) == nullptr)
      distributionMap_[s1]->setNamespace(s1 + "_" + distributionMap_[s1]->getNamespace());
    else
      distributionMap_[s1]->setNamespace(s1 + "_");

    auto constr = model->parameter(s2).getConstraint();
    if (constr)
      distributionMap_[s1]->restrictToConstraint(*constr);
  }

  // Initialization of modelsContainer_.

  c = 1;

  for (auto& it : distributionMap_)
  {
    c *= it.second->getNumberOfCategories();
  }

  for (i = 0; i < c; i++)
  {
    modelsContainer_.push_back(unique_ptr<TransitionModelInterface>(model->clone()));
    vProbas_.push_back(1.0 / static_cast<double>(c));
    vRates_.push_back(1.0);
  }
  // Initialization of parameters_.


  DiscreteDistributionInterface* pd;

  for (auto& it : distributionMap_)
  {
    const Parameter& pm = model->parameter(model->getParameterNameWithoutNamespace(getParameterNameWithoutNamespace(it.first)));
    pd = it.second.get();
    if (pm.hasConstraint())
    {
      pd->restrictToConstraint(*pm.getConstraint());
    }

    if (!dynamic_cast<ConstantDistribution*>(it.second.get()))
    {
      const ParameterList pl = pd->getParameters();
      for (i = 0; i != it.second->getNumberOfParameters(); ++i)
      {
        addParameter_(pl[i].clone());
      }
    }
    else
    {
      // addParameter_(new Parameter(it.first, pd->getCategory(0), (pd->parameter("value").getConstraint()) ? shared_ptr<ConstraintInterface>(pd->parameter("value").getConstraint()->clone()) : nullptr));
      addParameter_(new Parameter(it.first, pd->getCategory(0), pm.hasConstraint() ? shared_ptr<ConstraintInterface>(pm.getConstraint()->clone()) : nullptr));
    }
  }
  updateMatrices_();
}

MixtureOfATransitionModel::MixtureOfATransitionModel(const MixtureOfATransitionModel& msm) :
  AbstractParameterAliasable(msm),
  AbstractTransitionModel(msm),
  AbstractMixedTransitionModel(msm),
  distributionMap_(),
  from_(msm.from_),
  to_(msm.to_)
{
  for (auto& it : msm.distributionMap_)
  {
    distributionMap_[it.first].reset(it.second->clone());
  }
}

MixtureOfATransitionModel& MixtureOfATransitionModel::operator=(const MixtureOfATransitionModel& msm)
{
  AbstractParameterAliasable::operator=(msm);
  AbstractMixedTransitionModel::operator=(msm);
  from_ = msm.from_;
  to_ = msm.to_;

  // Clear existing containers:
  distributionMap_.clear();

  // Now copy new containers:
  for (auto& it : msm.distributionMap_)
  {
    distributionMap_[it.first].reset(it.second->clone());
  }
  return *this;
}


MixtureOfATransitionModel::~MixtureOfATransitionModel() {}

void MixtureOfATransitionModel::updateMatrices_()
{
  string s, t;
  size_t j, l;
  double d;
  ParameterList pl;

  // Update of distribution parameters from the parameters_ member
  // data. (reverse operation compared to what has been done in the
  // constructor).
  //  vector<string> v=getParameters().getParameterNames();

  for (auto& distrib : distributionMap_)
  {
    if (dynamic_cast<ConstantDistribution*>(distrib.second.get()) == nullptr)
    {
      vector<string> vDistnames = distrib.second->getParameters().getParameterNames();
      for (auto& parname : vDistnames)
      {
        d = getParameterValue(getParameterNameWithoutNamespace(parname));
        pl.addParameter(Parameter(parname, d));
      }
      distrib.second->matchParametersValues(pl);
      pl.reset();
    }
    else
    {
      t = distrib.second->getNamespace();
      d = parameter(getParameterNameWithoutNamespace(t.substr(0, t.length() - 1))).getValue();
      distrib.second->setParameterValue("value", d);
    }
  }

  for (size_t i = 0; i < modelsContainer_.size(); i++)
  {
    vProbas_[i] = 1;
    j = i;
    for (auto& distrib : distributionMap_)
    {
      s = distrib.first;
      l = j % distrib.second->getNumberOfCategories();

      d = distrib.second->getCategory(l);
      vProbas_[i] *= distrib.second->getProbability(l);
      if (pl.hasParameter(s))
        pl.setParameterValue(s, d);
      else
        pl.addParameter(Parameter(s, d));

      j = j / distrib.second->getNumberOfCategories();
    }

    modelsContainer_[i]->matchParametersValues(pl);
  }

  //  setting the equilibrium freqs
  for (size_t i = 0; i < getNumberOfStates(); i++)
  {
    freq_[i] = 0;
    for (j = 0; j < modelsContainer_.size(); j++)
    {
      freq_[i] += vProbas_[j] * modelsContainer_[j]->freq(i);
    }
  }
}

void MixtureOfATransitionModel::setFreq(std::map<int, double>& m)
{
  modelsContainer_[0]->setFreq(m);
  matchParametersValues(modelsContainer_[0]->getParameters());
}

Vuint MixtureOfATransitionModel::getSubmodelNumbers(const string& desc) const
{
  vector<string> parnames = modelsContainer_[0]->getParameters().getParameterNames();
  std::map<std::string, size_t> msubn;

  StringTokenizer st(desc, ",");
  while (st.hasMoreToken())
  {
    string param = st.nextToken();
    string::size_type index = param.rfind("_");
    if (index == string::npos)
      throw Exception("MixtureOfATransitionModel::getSubmodelNumbers parameter description should contain a number " + param);
    msubn[param.substr(0, index)] = TextTools::to<size_t>(param.substr(index + 1, 4)) - 1;
  }

  Vuint submodnb;
  size_t l;
  string s;

  for (size_t i = 0; i < modelsContainer_.size(); ++i)
  {
    size_t stopped = 0;
    size_t j = i;
    for (auto& it : distributionMap_)
    {
      s = it.first;
      l = j % it.second->getNumberOfCategories();

      if (msubn.find(s) != msubn.end())
      {
        if (msubn[s] == l)
        {
          stopped++;
        }
      }

      j = j / it.second->getNumberOfCategories();
    }
    if (stopped == msubn.size())// All requests are fulfilled
      submodnb.push_back(uint(i));
  }

  return submodnb;
}
