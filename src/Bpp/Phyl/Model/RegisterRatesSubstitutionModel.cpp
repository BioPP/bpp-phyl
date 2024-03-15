// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "RegisterRatesSubstitutionModel.h"

using namespace bpp;
using namespace std;

RegisterRatesSubstitutionModel::RegisterRatesSubstitutionModel(
    unique_ptr<SubstitutionModelInterface> originalModel,
    const SubstitutionRegisterInterface& reg,
    bool isNormalized) :
  AbstractParameterAliasable("FromRegister."),
  AbstractWrappedModel("FromRegister."),
  AbstractWrappedTransitionModel("FromRegister."),
  AbstractWrappedSubstitutionModel("FromRegister."),
  AbstractSubstitutionModel(originalModel->getAlphabet(), originalModel->getStateMap(), "FromRegister."),
  originalModel_(originalModel->clone()),
  registerName_(reg.getName()),
  vRegStates_(),
  nbTypes_(reg.getNumberOfSubstitutionTypes()),
  vRates_(reg.getNumberOfSubstitutionTypes())
{
//  getSubstitutionModel().enableEigenDecomposition(false);

  // record register types
  isScalable_ = isNormalized;


  vRegStates_.resize(nbTypes_);
  for (auto& vreg : vRegStates_)
  {
    vreg.resize(getNumberOfStates());
  }

  for (size_t i = 0; i < getNumberOfStates(); i++)
  {
    for (size_t j = 0; j < getNumberOfStates(); j++)
    {
      if (i != j)
      {
        size_t nR = reg.getType(i, j);
        if (nR != 0)
          vRegStates_[nR - 1][i].push_back((unsigned int)j);
      }
    }
  }

  /// !!! The order of the inclusion of parameters should not be
  // changed (see updateMatrices for vRates_ update).
  // rates for all register types
  for (size_t i = 1; i <= nbTypes_; i++)
  {
    addParameter_(new Parameter("FromRegister.rho_" + reg.getTypeName(i), 1, Parameter::R_PLUS_STAR));
  }

  model_().setNamespace(getNamespace() + model().getNamespace());
  addParameters_(model().getParameters());
  updateMatrices_();
}

/******************************************************************************/

void RegisterRatesSubstitutionModel::updateMatrices_()
{
  for (size_t i = 0; i < vRates_.size(); ++i)
  {
    vRates_[i] = getParameter_(i).getValue();
  }

  RowMatrix<double>& gen = generator_;

  gen = substitutionModel().generator();

  for (size_t t = 0; t < nbTypes_; ++t)
  {
    double rate = vRates_[t];
    const VVuint& v = vRegStates_[t];

    for (size_t i = 0; i < getNumberOfStates(); ++i)
    {
      const Vuint& v2 = v[i];
      for (const auto& j : v2)
      {
        gen(i, j) *= rate;
      }
    }
  }

  setDiagonal();

  AbstractSubstitutionModel::updateMatrices_();
}
