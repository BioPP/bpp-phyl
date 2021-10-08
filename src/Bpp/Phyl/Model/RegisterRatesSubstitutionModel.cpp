//
// File: RegisterRatesSubstitutionModel.cpp
// Created by: Laurent Gueguen
// Created on: samedi 24 octobre 2015, à 18h 50
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

#include "RegisterRatesSubstitutionModel.h"
#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace bpp;
using namespace std;

RegisterRatesSubstitutionModel::RegisterRatesSubstitutionModel(const SubstitutionModel& originalModel, const SubstitutionRegister& reg, bool isNormalized) :
  AbstractParameterAliasable("FromRegister."),
  AbstractSubstitutionModel(originalModel.getAlphabet(), originalModel.shareStateMap(), "FromRegister."),
  originalModel_(originalModel.clone()),
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
    addParameter_(new Parameter("FromRegister.rho_" + reg.getTypeName(i), 1, Parameter::R_PLUS));
  }

  getModel().setNamespace(getNamespace() + getModel().getNamespace());
  addParameters_(getModel().getParameters());
  updateMatrices();
}

/******************************************************************************/

void RegisterRatesSubstitutionModel::updateMatrices()
{
  for (size_t i = 0; i < vRates_.size(); i++)
  {
    vRates_[i] = getParameter_(i).getValue();
  }

  RowMatrix<double>& gen = generator_;

  gen = getSubstitutionModel().getGenerator();

  for (size_t t = 0; t < nbTypes_; t++)
  {
    double rate = vRates_[t];
    const VVuint& v = vRegStates_[t];

    for (size_t i = 0; i < getNumberOfStates(); i++)
    {
      const Vuint& v2 = v[i];
      for (const auto& j : v2)
      {
        gen(i, j) *= rate;
      }
    }
  }

  setDiagonal();

  AbstractSubstitutionModel::updateMatrices();
}
