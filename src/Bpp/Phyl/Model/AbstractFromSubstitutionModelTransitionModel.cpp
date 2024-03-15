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
  //AbstractParameterAliasable(prefix + subModel->getNamespace()),
  subModel_(std::move(subModel)),
  size_(subModel_->getNumberOfStates()),
  pij_t(size_, size_),
  dpij_t(size_, size_),
  d2pij_t(size_, size_),
  nestedPrefix_(subModel_->getNamespace())
{
  subModel_->setNamespace(getNamespace() + nestedPrefix_);
  addParameters_(subModel_->getParameters());
}


/******************************************************************************/

AbstractFromSubstitutionModelTransitionModel::AbstractFromSubstitutionModelTransitionModel(const AbstractFromSubstitutionModelTransitionModel& fmsm) :
  //AbstractParameterAliasable(fmsm),
  subModel_(fmsm.subModel_->clone()),
  size_(fmsm.size_),
  pij_t(fmsm.pij_t),
  dpij_t(fmsm.dpij_t),
  d2pij_t(fmsm.d2pij_t),
  nestedPrefix_(fmsm.nestedPrefix_)
{}


/******************************************************************************/

AbstractFromSubstitutionModelTransitionModel& AbstractFromSubstitutionModelTransitionModel::operator=(const AbstractFromSubstitutionModelTransitionModel& fmsm)
{
  AbstractParameterAliasable::operator=(fmsm);

  subModel_ = std::unique_ptr<SubstitutionModelInterface>(fmsm.subModel_->clone());
  size_ = fmsm.size_;
  pij_t = fmsm.pij_t;
  dpij_t = fmsm.dpij_t;
  d2pij_t = fmsm.d2pij_t;
  nestedPrefix_ = fmsm.nestedPrefix_;

  return *this;
}
