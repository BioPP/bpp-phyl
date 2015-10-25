//
// File: MixtureOfSubstitutionModels.cpp
// Created by: Laurent Gueguen
// Date: mardi 14 septembre 2010, à 20h 43
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

#include "MixtureOfSubstitutionModels.h"

#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Text/TextTools.h>

#include <string>

using namespace bpp;
using namespace std;

MixtureOfSubstitutionModels::MixtureOfSubstitutionModels(
  const Alphabet* alpha,
  vector<SubstitutionModel*> vpModel) :
  AbstractParameterAliasable("Mixture."),
  AbstractMixedSubstitutionModel(alpha, vpModel[0]->getStateMap().clone(), "Mixture.")
{
  size_t i, nbmod = vpModel.size();

  for (i = 0; i < nbmod; i++)
  {
    if (!vpModel[i])
      throw Exception("Empty model number " + TextTools::toString(i) + " in MixtureOfSubstitutionModels constructor");
    for (size_t j = i + 1; j < nbmod; j++)
    {
      if (vpModel[i] == vpModel[j])
        throw Exception("Same model at positions " + TextTools::toString(i) + " and " +
                        TextTools::toString(j) + " in MixtureOfSubstitutionModels constructor");
    }
  }

  // Initialization of modelsContainer_.

  for (i = 0; i < nbmod; i++)
  {
    modelsContainer_.push_back(vpModel[i]);
    vProbas_.push_back(1.0 / static_cast<double>(nbmod));
    vRates_.push_back(1.0);
  }

  // Initialization of parameters_.

  // relative rates and probas
  for (i = 0; i < nbmod - 1; i++)
  {
    addParameter_(new Parameter("Mixture.relproba" + TextTools::toString(i + 1), 1.0 / static_cast<double>(nbmod - i), &Parameter::PROP_CONSTRAINT_EX));
    addParameter_(new Parameter("Mixture.relrate" + TextTools::toString(i + 1), 1.0 / static_cast<double>(nbmod - i), &Parameter::PROP_CONSTRAINT_EX));
  }

  // models parameters

  for (i = 0; i < nbmod; i++)
  {
    modelsContainer_[i]->setNamespace("Mixture." + TextTools::toString(i + 1) + "_" + vpModel[i]->getNamespace());
    addParameters_(vpModel[i]->getParameters());
  }

  updateMatrices();
}

MixtureOfSubstitutionModels::MixtureOfSubstitutionModels(
    const Alphabet* alpha,
    vector<SubstitutionModel*> vpModel,
    Vdouble& vproba,
    Vdouble& vrate) :
  AbstractParameterAliasable("Mixture."),
  AbstractMixedSubstitutionModel(alpha, vpModel[0]->getStateMap().clone(), "Mixture.")
{
  size_t i, nbmod = vpModel.size();

  for (i = 0; i < nbmod; i++)
  {
    if (!vpModel[i])
      throw Exception("Empty model number " + TextTools::toString(i) + " in MixtureOfSubstitutionModels constructor");
    for (size_t j = i + 1; j < nbmod; j++)
    {
      if (vpModel[i] == vpModel[j])
        throw Exception("Same model at positions " + TextTools::toString(i) + " and " +
                        TextTools::toString(j) + " in MixtureOfSubstitutionModels constructor");
    }
  }

  double x = 0;
  double y = 0;

  for (i = 0; i < nbmod; i++)
  {
    if (vrate[i] <= 0)
      throw Exception("Non positive rate: " + TextTools::toString(vrate[i]) + " in MixtureOfSubstitutionModels constructor.");
    if (vproba[i] <= 0)
      throw Exception("Non positive probability: " + TextTools::toString(vproba[i]) + " in MixtureOfSubstitutionModels constructor.");
    x += vproba[i];
    y += vproba[i] * vrate[i];
  }

  if (fabs(1. - x) > NumConstants::SMALL())
    throw Exception("Probabilities must equal 1 (sum = " + TextTools::toString(x) + ").");
  if (fabs(1. - y) > NumConstants::SMALL())
    throw Exception("Expectation on rates must equal 1 (E =" + TextTools::toString(y) + ").");


  // Initialization of modelsContainer_.

  for (i = 0; i < nbmod; i++)
  {
    modelsContainer_.push_back(vpModel[i]);
  }

  // rates & probas

  for (i = 0; i < nbmod; i++)
  {
    vProbas_.push_back(1.0 / static_cast<double>(nbmod));
    vRates_.push_back(1.0);
  }

  // Initialization of parameters_.


  // relative rates and probas
  x = 0; y = 0;

  for (i = 0; i < nbmod - 1; i++)
  {
    addParameter_(new Parameter("Mixture.relproba" + TextTools::toString(i + 1), vproba[i] / (1 - x), &Parameter::PROP_CONSTRAINT_EX));
    x += vproba[i];
    addParameter_(new Parameter("Mixture.relrate" + TextTools::toString(i + 1), vproba[i] * vrate[i] / (1 - y), &Parameter::PROP_CONSTRAINT_EX));
    y += vproba[i] * vrate[i];
  }

  // models parameters

  for (i = 0; i < nbmod; i++)
  {
    modelsContainer_[i]->setNamespace("Mixture." + TextTools::toString(i + 1) + "_" + vpModel[i]->getNamespace());
    addParameters_(vpModel[i]->getParameters());
  }

  updateMatrices();
}

MixtureOfSubstitutionModels::MixtureOfSubstitutionModels(const MixtureOfSubstitutionModels& msm) :
  AbstractParameterAliasable(msm),
  AbstractMixedSubstitutionModel(msm)
{}

MixtureOfSubstitutionModels& MixtureOfSubstitutionModels::operator=(const MixtureOfSubstitutionModels& msm)
{
  AbstractMixedSubstitutionModel::operator=(msm);

  return *this;
}


MixtureOfSubstitutionModels::~MixtureOfSubstitutionModels()
{}

const SubstitutionModel* MixtureOfSubstitutionModels::getSubModelWithName(const std::string& name) const
{
  size_t nbmod=getNumberOfModels();

  for (size_t i=0; i<nbmod; i++)
    if (getNModel(i)->getName()==name)
      return getNModel(i);

  return NULL;
}

void MixtureOfSubstitutionModels::updateMatrices()
{
  size_t i, j, nbmod = modelsContainer_.size();

  double x, y;
  x = 1.0;

  for (i = 0; i < nbmod - 1; i++)
  {
    y = getParameterValue("relproba" + TextTools::toString(i + 1));
    vProbas_[i] = x * y;
    x *= 1 - y;
  }
  vProbas_[nbmod - 1] = x;

  x = 1.0;
  bool approx = false; // used when some categories are avoided
  double s = 0;
  for (i = 0; i < nbmod - 1; i++)
  {
    y = getParameterValue("relrate" + TextTools::toString(i + 1));
    if (vProbas_[i] < NumConstants::SMALL())
    {
      vRates_[i] = NumConstants::SMALL();
      approx = true;
    }
    else
    {
      vRates_[i] = x * y / vProbas_[i];
      s += x * y;
    }
    x *= 1 - y;
  }

  if (vProbas_[nbmod - 1] < NumConstants::SMALL())
  {
    vRates_[nbmod - 1] = NumConstants::SMALL();
    approx = true;
  }
  else
  {
    vRates_[nbmod - 1] = x / vProbas_[nbmod - 1];
    s += x;
  }

  if (approx)
    for (i = 0; i < nbmod; i++)
    {
      vRates_[i] /= s;
    }

  // / models

  for (i = 0; i < nbmod; i++)
  {
    modelsContainer_[i]->setRate(rate_ * vRates_[i]);
    modelsContainer_[i]->matchParametersValues(getParameters());
  }

  // / freq_

  for (i = 0; i < getNumberOfStates(); i++)
  {
    freq_[i] = 0;
    for (j = 0; j < modelsContainer_.size(); j++)
    {
      freq_[i] += vProbas_[j] * modelsContainer_[j]->freq(i);
    }
  }
}


void MixtureOfSubstitutionModels::setFreq(std::map<int, double>& m)
{
  ParameterList pl;
  for (unsigned int n = 0; n < modelsContainer_.size(); n++)
  {
    modelsContainer_[n]->setFreq(m);
    pl.addParameters(modelsContainer_[n]->getParameters());
  }
  matchParametersValues(pl);
}

void MixtureOfSubstitutionModels::setVRates(const Vdouble& vd)
{
  AbstractMixedSubstitutionModel::setVRates(vd);

  size_t i, nbmod = modelsContainer_.size();
  double sP = 0;
  for (i = 0; i < nbmod - 1; i++)
  {
    sP += vProbas_[i];
  }

  double y = 0;
  for (i = 0; i < nbmod - 1; i++)
  {
    setParameterValue("relrate" + TextTools::toString(i + 1), vProbas_[i] / sP * vRates_[i] / (1 - y));
    y += vProbas_[i] / sP * vRates_[i];
  }
}

Vint MixtureOfSubstitutionModels::getSubmodelNumbers(string& desc) const
{
  size_t i;
  for (i = 0; i < getNumberOfModels(); i++)
  {
    if (getNModel(i)->getName() == desc)
      break;
  }
  if (i == getNumberOfModels())
    throw Exception("MixtureOfSubstitutionModels::getSubmodelNumbers model description do not match " + desc);

  Vint submodnb;
  submodnb.push_back(static_cast<int>(i));

  return submodnb;
}
