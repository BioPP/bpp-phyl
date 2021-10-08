//
// File: FromMixtureSubstitutionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: samedi 24 octobre 2015, à 18h 50
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)
  
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


#include "FromMixtureSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

FromMixtureSubstitutionModel::FromMixtureSubstitutionModel(const MixedTransitionModel& mixedModel, const std::string& subModelName, const std::string& mixtDesc) :
  AbstractParameterAliasable(mixedModel.getName() + "_" + subModelName + "."),
  subModel_(),
  mixtName_(mixtDesc)
{
  const TransitionModel* tm = mixedModel.getModel(subModelName);

  if (tm == 0)
    throw ParameterNotFoundException("FromMixtureSubstitutionModel::FromMixtureSubstitutionModel : unknown model name", subModelName);

  const auto sm = dynamic_cast<const SubstitutionModel*>(tm);
  if (sm == 0)
    throw Exception("FromMixtureSubstitutionModel::FromMixtureSubstitutionModel : model " + subModelName + " is not a substitution model.");

  subModel_ = std::unique_ptr<SubstitutionModel>(sm->clone());
  subModel_->setNamespace(getNamespace());

  subModel_->setRate(1);
  addParameters_(subModel_->getParameters());
}

FromMixtureSubstitutionModel::FromMixtureSubstitutionModel(const MixedTransitionModel& mixedModel, size_t subModelNumber, const std::string& mixtDesc) :
  AbstractParameterAliasable(mixedModel.getName() + "_" + TextTools::toString(subModelNumber) + "."),
  subModel_(),
  mixtName_(mixtDesc)
{
  if (subModelNumber >= mixedModel.getNumberOfModels())
    throw ParameterNotFoundException("FromMixtureSubstitutionModel::FromMixtureSubstitutionModel : bad model number", TextTools::toString(subModelNumber));

  const auto tm = mixedModel.getNModel(subModelNumber);
  const auto sm = dynamic_cast<const SubstitutionModel*>(tm);

  if (sm == 0)
    throw Exception("FromMixtureSubstitutionModel::FromMixtureSubstitutionModel : model with number " + TextTools::toString(subModelNumber) + " is not a substitution model.");

  subModel_ = std::unique_ptr<SubstitutionModel>(sm->clone());
  subModel_->setNamespace(getNamespace());

  subModel_->setRate(1);
  addParameters_(subModel_->getParameters());
}


/******************************************************************************/

FromMixtureSubstitutionModel::FromMixtureSubstitutionModel(const FromMixtureSubstitutionModel& fmsm) :
  AbstractParameterAliasable(fmsm),
  subModel_(fmsm.subModel_->clone()),
  mixtName_(fmsm.mixtName_)
{}


/******************************************************************************/

FromMixtureSubstitutionModel& FromMixtureSubstitutionModel::operator=(const FromMixtureSubstitutionModel& fmsm)
{
  AbstractParameterAliasable::operator=(fmsm);

  subModel_ = std::unique_ptr<SubstitutionModel>(fmsm.subModel_->clone());

  mixtName_ = fmsm.mixtName_;

  return *this;
}
