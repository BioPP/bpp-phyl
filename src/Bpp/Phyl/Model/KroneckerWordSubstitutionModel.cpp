// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/EigenValue.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>

#include "KroneckerWordSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/WordAlphabet.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/******************************************************************************/

KroneckerWordSubstitutionModel::KroneckerWordSubstitutionModel(
    ModelList& modelList,
    const string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Kron." : prefix),
  AbstractKroneckerWordSubstitutionModel(
      modelList,
      (prefix == "") ? "Kron." : prefix)
{
  updateMatrices_();
}

KroneckerWordSubstitutionModel::KroneckerWordSubstitutionModel(
    std::shared_ptr<const Alphabet> alph,
    std::shared_ptr<const StateMapInterface> stateMap,
    const string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Kron." : prefix),
  AbstractKroneckerWordSubstitutionModel(alph, stateMap, (prefix == "") ? "Kron." : prefix)
{
  enableEigenDecomposition(true);
}

KroneckerWordSubstitutionModel::KroneckerWordSubstitutionModel(
    std::unique_ptr<SubstitutionModelInterface> pmodel,
    unsigned int num,
    const string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Kron." : prefix),
  AbstractKroneckerWordSubstitutionModel(std::move(pmodel),
      num,
      (prefix == "") ? "Kron." : prefix)
{
  enableEigenDecomposition(true);
  updateMatrices_();
}

KroneckerWordSubstitutionModel::KroneckerWordSubstitutionModel(
    ModelList& modelList,
    const vector<set< size_t>>& vPos,
    const string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Kron." : prefix),
  AbstractKroneckerWordSubstitutionModel(
      modelList, vPos,
      (prefix == "") ? "Kron." : prefix)
{
  enableEigenDecomposition(true);
  updateMatrices_();
}

KroneckerWordSubstitutionModel::KroneckerWordSubstitutionModel(
    std::unique_ptr<SubstitutionModelInterface> pmodel,
    unsigned int num,
    const vector<set< size_t>>& vPos,
    const string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Kron." : prefix),
  AbstractKroneckerWordSubstitutionModel(std::move(pmodel),
      num, vPos,
      (prefix == "") ? "Kron." : prefix)
{
  enableEigenDecomposition(true);
  updateMatrices_();
}


string KroneckerWordSubstitutionModel::getName() const
{
  return "Kron";
}
