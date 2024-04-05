// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "FromMixtureSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

FromMixtureSubstitutionModel::FromMixtureSubstitutionModel(
    const MixedTransitionModelInterface& mixedModel,
    const std::string& subModelName,
    const std::string& mixtDesc) :
  AbstractParameterAliasable(mixedModel.getName() + "_" + subModelName + "."),
  AbstractWrappedModel(mixedModel.getName() + "_" + subModelName + "."),
  AbstractWrappedTransitionModel(mixedModel.getName() + "_" + subModelName + "."),
  AbstractTotallyWrappedTransitionModel(mixedModel.getName() + "_" + subModelName + "."),
  AbstractWrappedSubstitutionModel(mixedModel.getName() + "_" + subModelName + "."),
  AbstractTotallyWrappedSubstitutionModel(mixedModel.getName() + "_" + subModelName + "."),
  subModel_(),
  mixtName_(mixtDesc)
{
  try
  {
    auto& tm = mixedModel.model(subModelName);
    try
    {
      auto& sm = dynamic_cast<const SubstitutionModelInterface&>(tm);
      subModel_ = std::unique_ptr<SubstitutionModelInterface>(sm.clone());
      subModel_->setNamespace(getNamespace());
      subModel_->setRate(1);
      addParameters_(subModel_->getParameters());
    }
    catch (bad_cast&)
    {
      throw Exception("FromMixtureSubstitutionModel::FromMixtureSubstitutionModel : model " + subModelName + " is not a substitution model.");
    }
  }
  catch (NullPointerException&)
  {
    throw ParameterNotFoundException("FromMixtureSubstitutionModel::FromMixtureSubstitutionModel : unknown model name", subModelName);
  }
}

FromMixtureSubstitutionModel::FromMixtureSubstitutionModel(
    const MixedTransitionModelInterface& mixedModel,
    size_t subModelNumber,
    const std::string& mixtDesc) :
  AbstractParameterAliasable(mixedModel.getName() + "_" + TextTools::toString(subModelNumber) + "."),
  AbstractWrappedModel(mixedModel.getName() + "_" + TextTools::toString(subModelNumber) + "."),
  AbstractWrappedTransitionModel(mixedModel.getName() + "_" + TextTools::toString(subModelNumber) + "."),
  AbstractTotallyWrappedTransitionModel(mixedModel.getName() + "_" + TextTools::toString(subModelNumber) + "."),
  AbstractWrappedSubstitutionModel(mixedModel.getName() + "_" + TextTools::toString(subModelNumber) + "."),
  AbstractTotallyWrappedSubstitutionModel(mixedModel.getName() + "_" + TextTools::toString(subModelNumber) + "."),
  subModel_(),
  mixtName_(mixtDesc)
{
  if (subModelNumber >= mixedModel.getNumberOfModels())
    throw ParameterNotFoundException("FromMixtureSubstitutionModel::FromMixtureSubstitutionModel : bad model number", TextTools::toString(subModelNumber));

  auto& tm = mixedModel.nModel(subModelNumber);
  try
  {
    auto& sm = dynamic_cast<const SubstitutionModelInterface&>(tm);
    subModel_ = std::unique_ptr<SubstitutionModelInterface>(sm.clone());
    subModel_->setNamespace(getNamespace());
    subModel_->setRate(1);
    addParameters_(subModel_->getParameters());
  }
  catch (bad_cast&)
  {
    throw Exception("FromMixtureSubstitutionModel::FromMixtureSubstitutionModel : model with number " + TextTools::toString(subModelNumber) + " is not a substitution model.");
  }
}


/******************************************************************************/

FromMixtureSubstitutionModel::FromMixtureSubstitutionModel(const FromMixtureSubstitutionModel& fmsm) :
  AbstractParameterAliasable(fmsm),
  AbstractWrappedModel(fmsm),
  AbstractWrappedTransitionModel(fmsm),
  AbstractTotallyWrappedTransitionModel(fmsm),
  AbstractWrappedSubstitutionModel(fmsm),
  AbstractTotallyWrappedSubstitutionModel(fmsm),
  subModel_(fmsm.subModel_->clone()),
  mixtName_(fmsm.mixtName_)
{}


/******************************************************************************/

FromMixtureSubstitutionModel& FromMixtureSubstitutionModel::operator=(const FromMixtureSubstitutionModel& fmsm)
{
  AbstractTotallyWrappedSubstitutionModel::operator=(fmsm);

  subModel_.reset(fmsm.subModel_->clone());
  mixtName_ = fmsm.mixtName_;

  return *this;
}
