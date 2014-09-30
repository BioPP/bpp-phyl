//
// File: AbstractMixedSubstitutionModel
// Created by: Laurent Gueguen
// On: vendredi 19 novembre 2010, à 15h 55
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

#include "AbstractMixedSubstitutionModel.h"

#include <string>

using namespace bpp;
using namespace std;


AbstractMixedSubstitutionModel::AbstractMixedSubstitutionModel(
  const Alphabet* alpha, StateMap* stateMap, const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractSubstitutionModel(alpha, stateMap, prefix),
  modelsContainer_(),
  vProbas_(),
  vRates_()
{
  for (unsigned int i = 0; i < size_; i++)
  {
    for (unsigned int j = 0; j < size_; j++)
    {
      generator_(i, j) = 0;
      exchangeability_(i, j) = 0;
      leftEigenVectors_(i, j) = 0;
      rightEigenVectors_(i, j) = 0;
    }
    eigenValues_[i] = 0;
  }
  eigenDecompose_ = false;
}

AbstractMixedSubstitutionModel::AbstractMixedSubstitutionModel(const AbstractMixedSubstitutionModel& msm) :
  AbstractParameterAliasable(msm),
  AbstractSubstitutionModel(msm),
  modelsContainer_(),
  vProbas_(),
  vRates_()
{
  for (unsigned int i = 0; i < msm.modelsContainer_.size(); i++)
  {
    modelsContainer_.push_back(msm.modelsContainer_[i]->clone());
    vProbas_.push_back(msm.vProbas_[i]);
    vRates_.push_back(msm.vRates_[i]);
  }
}

AbstractMixedSubstitutionModel& AbstractMixedSubstitutionModel::operator=(const AbstractMixedSubstitutionModel& model)
{
  AbstractParameterAliasable::operator=(model);
  AbstractSubstitutionModel::operator=(model);

  // Clear existing containers:
  modelsContainer_.clear();
  vProbas_.clear();
  vRates_.clear();

  for (unsigned int i = 0; i < model.modelsContainer_.size(); i++)
  {
    modelsContainer_.push_back(model.modelsContainer_[i]->clone());
    vProbas_.push_back(model.vProbas_[i]);
    vRates_.push_back(model.vRates_[i]);
  }

  return *this;
}

AbstractMixedSubstitutionModel::~AbstractMixedSubstitutionModel()
{
  for (unsigned int i = 0; i < modelsContainer_.size(); i++)
  {
    delete modelsContainer_[i];
  }
}

size_t AbstractMixedSubstitutionModel::getNumberOfStates() const
{
  return modelsContainer_[0]->getNumberOfStates();
}

const Matrix<double>& AbstractMixedSubstitutionModel::getPij_t(double t) const
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


const Matrix<double>& AbstractMixedSubstitutionModel::getdPij_dt(double t) const
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


const Matrix<double>& AbstractMixedSubstitutionModel::getd2Pij_dt2(double t) const
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


void AbstractMixedSubstitutionModel::setRate(double rate)
{
  AbstractSubstitutionModel::setRate(rate);

  double sum = 0;
  double sP = 0;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
  {
    sum += vRates_[n] * vProbas_[n];
    sP += vProbas_[n];
  }
  sum /= sP;

  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
  {
    vRates_[n] *= rate_ / sum;
    modelsContainer_[n]->setRate(vRates_[n]);
  }
}

void AbstractMixedSubstitutionModel::setVRates(const Vdouble& vd)
{
  if (vd.size() != modelsContainer_.size())
    throw Exception("AbstractMixedSubstitutionModel::setVRates  bad size of Vdouble argument.");

  for (unsigned int i = 0; i < vd.size(); i++)
  {
    vRates_[i] = vd[i];
  }

  normalizeVRates();
}

void AbstractMixedSubstitutionModel::normalizeVRates()
{
  double sum = 0;
  double sP = 0;
  for (unsigned int i = 0; i < vRates_.size(); i++)
  {
    sum += vRates_[i] * vProbas_[i];
    sP += vProbas_[i];
  }
  sum /= sP;

  for (unsigned int i = 0; i < vRates_.size(); i++)
  {
    vRates_[i] *= rate_ / sum;
    modelsContainer_[i]->setRate(vRates_[i]);
  }
}
