//
// File: CodonAdHocSubstitutionModel.cpp
// Created by:  Laurent Gueguen
// Created on: lundi 30 octobre 2017, à 06h 39
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

#include "CodonAdHocSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod,
    std::vector<CoreCodonSubstitutionModel*>& vpmodel,
    const std::string& name) :
  AbstractParameterAliasable(name+"."),
  AbstractCodonSubstitutionModel(gCode, pmod, name+"."),
  vModel_(),
  name_(name),
  freqSet_()
{
  for (auto model : vpmodel)
    if (model!=NULL)
    {
      model->setNamespace(name+".");
      if (model->getFrequencySet())
      {
        if (freqSet_)
          throw Exception("CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel : two sub models with FrequencySet");
        
        freqSet_=model->getFrequencySet();
      }
      vModel_.push_back(unique_ptr<CoreCodonSubstitutionModel>(model));
      addParameters_(model->getParameters());
    }
  computeFrequencies(true);
  updateMatrices();
}

CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod1,
    NucleotideSubstitutionModel* pmod2,
    NucleotideSubstitutionModel* pmod3,
    std::vector<CoreCodonSubstitutionModel*>& vpmodel,
    const std::string& name) :
  AbstractParameterAliasable(name+"."),
  AbstractCodonSubstitutionModel(gCode, pmod1, pmod2, pmod3, name+"."),
  vModel_(),
  name_(name),
  freqSet_()
{
  for (auto& model : vpmodel)
    if (model!=NULL)
    {
      model->setNamespace(name+".");
      if (model->getFrequencySet())
      {
        if (freqSet_)
          throw Exception("CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel : two sub models with FrequencySet");
        
        freqSet_=model->getFrequencySet();
      }
      
      vModel_.push_back(unique_ptr<CoreCodonSubstitutionModel>(model));
      addParameters_(model->getParameters());
    }

  computeFrequencies(true);
  updateMatrices();
}

CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel(const CodonAdHocSubstitutionModel& model) :
  AbstractParameterAliasable(model),
  AbstractCodonSubstitutionModel(model),
  vModel_(),
  name_(model.name_),
  freqSet_()
{
  for (auto& mod : model.vModel_)
    vModel_.emplace_back(mod->clone());

  for (auto& mod : vModel_)
    if (mod->getFrequencySet())
    {
      if (freqSet_)
        throw Exception("CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel : two sub models with FrequencySet");
        
      freqSet_=mod->getFrequencySet();
    }
}

CodonAdHocSubstitutionModel& CodonAdHocSubstitutionModel::operator=(const CodonAdHocSubstitutionModel& model)
{
  AbstractParameterAliasable::operator=(model);
  AbstractCodonSubstitutionModel::operator=(model);
  name_ = model.name_;
  
  vModel_.clear();
  freqSet_=0;  
  
  for (auto& mod : model.vModel_)
    vModel_.emplace_back(mod->clone());

  for (auto& mod : vModel_)
    if (mod->getFrequencySet())
    {
      if (freqSet_)
        throw Exception("CodonAdHocSubstitutionModel::CodonAdHocSubstitutionModel : two sub models with FrequencySet");
      
      freqSet_=mod->getFrequencySet();
    }

  return *this;
}

void CodonAdHocSubstitutionModel::setFreq(std::map<int,double>& frequencies)
{
  for (auto& model : vModel_)
    model->setFreq(frequencies);
}

void CodonAdHocSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  for (auto& model : vModel_)
    model->matchParametersValues(parameters);

  // Beware: must be called at the end
  AbstractCodonSubstitutionModel::fireParameterChanged(parameters);
}

double CodonAdHocSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  double x(1);

  for (auto& model : vModel_)
    x*= model->getCodonsMulRate(i,j);

  return x;
}

