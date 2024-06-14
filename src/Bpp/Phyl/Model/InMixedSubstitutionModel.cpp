// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "InMixedSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

InMixedSubstitutionModel::InMixedSubstitutionModel(
    unique_ptr<MixedTransitionModelInterface> mixedModel,
    const std::string& subModelName,
    const std::string& mixtDesc) :
  AbstractParameterAliasable(mixedModel->getNamespace()),
  AbstractWrappedModel(mixedModel->getNamespace()),
  AbstractWrappedTransitionModel(mixedModel->getNamespace()),
  AbstractWrappedSubstitutionModel(mixedModel->getNamespace()),
  mixedModelPtr_(std::move(mixedModel)),
  subModelNumber_(0),
  mixtName_(mixtDesc)
{
  Vuint vn = mixedModelPtr_->getSubmodelNumbers(subModelName);
  subModelNumber_ = (size_t)vn[0];
  addParameters_(mixedModelPtr_->getParameters());
}


InMixedSubstitutionModel::InMixedSubstitutionModel(
    unique_ptr<MixedTransitionModelInterface> mixedModel,
    size_t subModelNumber,
    const std::string& mixtDesc) :
  AbstractParameterAliasable(mixedModel->getNamespace()),
  AbstractWrappedModel(mixedModel->getNamespace()),
  AbstractWrappedTransitionModel(mixedModel->getNamespace()),
  AbstractWrappedSubstitutionModel(mixedModel->getNamespace()),
  mixedModelPtr_(std::move(mixedModel)),
  subModelNumber_(subModelNumber),
  mixtName_(mixtDesc)
{
  if (subModelNumber >= mixedModelPtr_->getNumberOfModels())
    throw ParameterNotFoundException("InMixedSubstitutionModel::InMixedSubstitutionModel : bad model number", TextTools::toString(subModelNumber));

  addParameters_(mixedModelPtr_->getParameters());
}


/******************************************************************************/

InMixedSubstitutionModel::InMixedSubstitutionModel(const InMixedSubstitutionModel& fmsm) :
  AbstractParameterAliasable(fmsm),
  AbstractWrappedModel(fmsm),
  AbstractWrappedTransitionModel(fmsm),
  AbstractWrappedSubstitutionModel(fmsm),
  mixedModelPtr_(fmsm.mixedModelPtr_->clone()),
  subModelNumber_(fmsm.subModelNumber_),
  mixtName_(fmsm.mixtName_)
{}


/******************************************************************************/

InMixedSubstitutionModel& InMixedSubstitutionModel::operator=(const InMixedSubstitutionModel& fmsm)
{
  AbstractWrappedSubstitutionModel::operator=(fmsm);

  mixedModelPtr_.reset(fmsm.mixedModelPtr_->clone());
  subModelNumber_ = fmsm.subModelNumber_;
  mixtName_ = fmsm.mixtName_;

  return *this;
}
