//
// File: KroneckerWordSubstitutionModel.cpp
// Created by: Laurent Gueguen
// Created on: mercredi 27 juillet 2016, à 00h 12
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

#include "KroneckerWordSubstitutionModel.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/EigenValue.h>

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
  const std::string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Kron." : prefix),
  AbstractKroneckerWordSubstitutionModel(
    modelList,
    (prefix == "") ? "Kron." : prefix)
{
  updateMatrices();
}

KroneckerWordSubstitutionModel::KroneckerWordSubstitutionModel(
  const Alphabet* alph,
  std::shared_ptr<const StateMap> stateMap,
  const std::string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Kron." : prefix),
  AbstractKroneckerWordSubstitutionModel(alph, stateMap, (prefix == "") ? "Kron." : prefix)
{
  enableEigenDecomposition(true);
}

KroneckerWordSubstitutionModel::KroneckerWordSubstitutionModel(
  SubstitutionModel* pmodel,
  unsigned int num,
  const std::string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Kron." : prefix),
  AbstractKroneckerWordSubstitutionModel(pmodel,
                                num,
                                (prefix == "") ? "Kron." : prefix)
{
  enableEigenDecomposition(true);
  updateMatrices();
}

KroneckerWordSubstitutionModel::KroneckerWordSubstitutionModel(
  ModelList& modelList,
  const std::vector<std::set< size_t> >& vPos,
  const std::string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Kron." : prefix),
  AbstractKroneckerWordSubstitutionModel(
    modelList, vPos,
    (prefix == "") ? "Kron." : prefix)
{
  enableEigenDecomposition(true);
  updateMatrices();
}

KroneckerWordSubstitutionModel::KroneckerWordSubstitutionModel(
  SubstitutionModel* pmodel,
  unsigned int num,
  const std::vector<std::set< size_t> >& vPos,
  const std::string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Kron." : prefix),
  AbstractKroneckerWordSubstitutionModel(pmodel,
                                         num, vPos,
                                         (prefix == "") ? "Kron." : prefix)
{
  enableEigenDecomposition(true);
  updateMatrices();
}


string KroneckerWordSubstitutionModel::getName() const
{
  return "Kron";
}


