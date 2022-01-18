//
// File: CompoundTransitionModel.cpp
// Authors:
//   Anaïs Prud'homme
// Date :
//   Vendredi 3 décembre 2021 à 11h30
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/VectorTools.h>
#include <string>

#include "CompoundTransitionModel.h"

using namespace bpp;
using namespace std;


CompoundTransitionModel::CompoundTransitionModel(
  const Alphabet* alpha,
  TransitionModel* model,
  int ffrom,
  int tto) :
  AbstractParameterAliasable(model->getNamespace()),
  AbstractTransitionModel(alpha, model->shareStateMap(), model->getNamespace()),
  from_(ffrom),
  to_(tto),
  modelsContainer_(),
  vProbas_()
{
  if (to_ >= int(alpha->getSize()))
    throw BadIntegerException("Bad state in alphabet", to_);
  if (from_ >= int(alpha->getSize()))
    throw BadIntegerException("Bad state in alphabet", from_);

  size_t c, i;
  string s1, s2, t;

  // Initialization of modelsContainer_.

  for (i = 0; i < c; i++)
  {
    modelsContainer_.push_back(std::shared_ptr<TransitionModel>(model->clone()));
    vProbas_.push_back(1.0 / static_cast<double>(c));
  }

  updateMatrices();
}

CompoundTransitionModel::CompoundTransitionModel(const CompoundTransitionModel& msm) :
  AbstractParameterAliasable(msm),
  AbstractTransitionModel(msm),
  from_(msm.from_),
  to_(msm.to_),
  modelsContainer_(),
  vProbas_(msm.vProbas_)
{
  for (unsigned int i = 0; i < msm.modelsContainer_.size(); i++)
  {
    modelsContainer_.push_back(std::shared_ptr<TransitionModel>(msm.modelsContainer_[i]->clone()));
  }

}

CompoundTransitionModel& CompoundTransitionModel::operator=(const CompoundTransitionModel& msm)
{
  AbstractParameterAliasable::operator=(msm);
  from_ = msm.from_;
  to_ = msm.to_;

  // Clear existing containers:
  modelsContainer_.clear();
  vProbas_.clear();

  for (unsigned int i = 0; i < msm.modelsContainer_.size(); i++)
  {
    modelsContainer_.push_back(std::shared_ptr<TransitionModel>(msm.modelsContainer_[i]->clone()));
    vProbas_.push_back(msm.vProbas_[i]);
  }
  return *this;
}
void CompoundTransitionModel::updateMatrices()
{
  string s, t;
  size_t j;
  ParameterList pl;

  for (size_t i = 0; i < modelsContainer_.size(); i++)
  {
    vProbas_[i] = 1;
    j = i;

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

void CompoundTransitionModel::setFreq(std::map<int, double>& m)
{
  modelsContainer_[0]->setFreq(m);
  matchParametersValues(modelsContainer_[0]->getParameters());
}

const TransitionModel* CompoundTransitionModel::getModel(const std::string& name) const
{
  size_t nbmod = getNumberOfModels();

  for (size_t i = 0; i < nbmod; i++)
  {
    if (getModel(i)->getName() == name)
      return getModel(i);
  }

  return NULL;
}

/*Vuint CompoundTransitionModel::getSubmodelNumbers(const string& desc) const
{
  vector<string> parnames = modelsContainer_[0]->getParameters().getParameterNames();
  std::map<std::string, size_t> msubn;

  StringTokenizer st(desc, ",");
  while (st.hasMoreToken())
  {
    string param = st.nextToken();
    string::size_type index = param.rfind("_");
    if (index == string::npos)
      throw Exception("CompoundTransitionModel::getSubmodelNumbers parameter description should contain a number " + param);
    msubn[param.substr(0, index)] = TextTools::to<size_t>(param.substr(index + 1, 4)) - 1;
  }

  //Vuint submodnb;
  //size_t i, j, l;
  //string s;

  //return submodnb;
}*/

size_t CompoundTransitionModel::getNumberOfStates() const
{
  return modelsContainer_[0]->getNumberOfStates();
}

const Matrix<double>& CompoundTransitionModel::getPij_t(double t) const
{
  vector<const Matrix<double>* > vM;
  double sP = 0;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
  {
    vM.push_back(&modelsContainer_[n]->getPij_t(t));
    sP += vProbas_[n];
  }

  for (unsigned int i = 0; i < getNumberOfStates(); i++)
  {
    for (unsigned int j = 0; j < getNumberOfStates(); j++)
    {
      double x = 0;
      for (unsigned int n = 0; n < modelsContainer_.size(); n++)
      {
        x += (*vM[n])(i, j) * vProbas_[n];
      }
      pijt_(i, j) = x / sP;
    }
  }
  return pijt_;
}


const Matrix<double>& CompoundTransitionModel::getdPij_dt(double t) const
{
  vector<const Matrix<double>* > vM;
  double sP = 0;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
  {
    vM.push_back(&modelsContainer_[n]->getdPij_dt(t));
    sP += vProbas_[n];
  }

  for (unsigned int i = 0; i < getNumberOfStates(); i++)
  {
    for (unsigned int j = 0; j < getNumberOfStates(); j++)
    {
      double x = 0;
      for (unsigned int n = 0; n < modelsContainer_.size(); n++)
      {
        x += (*vM[n])(i, j) * vProbas_[n];
      }
      dpijt_(i, j) = x / sP;
    }
  }
  return dpijt_;
}


const Matrix<double>& CompoundTransitionModel::getd2Pij_dt2(double t) const
{
  vector<const Matrix<double>* > vM;
  double sP = 0;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
  {
    vM.push_back(&modelsContainer_[n]->getd2Pij_dt2(t));
    sP += vProbas_[n];
  }

  for (unsigned int i = 0; i < getNumberOfStates(); i++)
  {
    for (unsigned int j = 0; j < getNumberOfStates(); j++)
    {
      double x = 0;
      for (unsigned int n = 0; n < modelsContainer_.size(); n++)
      {
        x += (*vM[n])(i, j) * vProbas_[n];
      }
      d2pijt_(i, j) = x / sP;
    }
  }
  return d2pijt_;
}
