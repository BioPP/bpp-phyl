// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractBiblioMixedTransitionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

AbstractBiblioMixedTransitionModel::AbstractBiblioMixedTransitionModel(const std::string& prefix) :
  AbstractBiblioTransitionModel(prefix),
  mixedModelPtr_()
{}

AbstractBiblioMixedTransitionModel::AbstractBiblioMixedTransitionModel(const AbstractBiblioMixedTransitionModel& mod2) :
  AbstractBiblioTransitionModel(mod2),
  mixedModelPtr_(mod2.mixedModelPtr_->clone())
{}

AbstractBiblioMixedTransitionModel& AbstractBiblioMixedTransitionModel::operator=(const AbstractBiblioMixedTransitionModel& mod2)
{
  AbstractBiblioTransitionModel::operator=(mod2);
  mixedModelPtr_.reset(mod2.mixedModelPtr_->clone());
  return *this;
}

AbstractBiblioMixedTransitionModel::~AbstractBiblioMixedTransitionModel()
{}

Vuint AbstractBiblioMixedTransitionModel::getSubmodelNumbers(const string& desc) const
{
  std::string desc2;

  StringTokenizer st(desc, ",");
  while (st.hasMoreToken())
  {
    string param = st.nextToken();

    desc2 += getParameterNameWithoutNamespace(param);
    if (st.hasMoreToken())
      desc2 += ",";
  }

  return mixedModel().getSubmodelNumbers(desc2);
}
