// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractFromSubstitutionModelTransitionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

AbstractFromSubstitutionModelTransitionModel::AbstractFromSubstitutionModelTransitionModel(
    unique_ptr<SubstitutionModelInterface> subModel,
    const std::string& prefix) :
  // AbstractParameterAliasable(prefix + subModel->getNamespace()),
  subModel_(std::move(subModel)),
  size_(subModel_->getNumberOfStates()),
  pijt_(size_, size_),
  dpijt_(size_, size_),
  d2pijt_(size_, size_),
  nestedPrefix_(subModel_->getNamespace())
{
  subModel_->setNamespace(getNamespace() + nestedPrefix_);
  addParameters_(subModel_->getParameters());
}


/******************************************************************************/

AbstractFromSubstitutionModelTransitionModel::AbstractFromSubstitutionModelTransitionModel(const AbstractFromSubstitutionModelTransitionModel& fmsm) :
  // AbstractParameterAliasable(fmsm),
  subModel_(fmsm.subModel_->clone()),
  size_(fmsm.size_),
  pijt_(fmsm.pijt_),
  dpijt_(fmsm.dpijt_),
  d2pijt_(fmsm.d2pijt_),
  nestedPrefix_(fmsm.nestedPrefix_)
{}


/******************************************************************************/

AbstractFromSubstitutionModelTransitionModel& AbstractFromSubstitutionModelTransitionModel::operator=(const AbstractFromSubstitutionModelTransitionModel& fmsm)
{
  AbstractParameterAliasable::operator=(fmsm);

  subModel_ = std::unique_ptr<SubstitutionModelInterface>(fmsm.subModel_->clone());
  size_ = fmsm.size_;
  pijt_ = fmsm.pijt_;
  dpijt_ = fmsm.dpijt_;
  d2pijt_ = fmsm.d2pijt_;
  nestedPrefix_ = fmsm.nestedPrefix_;

  return *this;
}
