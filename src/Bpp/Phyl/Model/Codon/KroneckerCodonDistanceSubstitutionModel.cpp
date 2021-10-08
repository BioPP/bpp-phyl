//
// File: KroneckerCodonDistanceSubstitutionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: vendredi 23 septembre 2016, Ã  12h 18
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
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


#include "KroneckerCodonDistanceSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

KroneckerCodonDistanceSubstitutionModel::KroneckerCodonDistanceSubstitutionModel(
  const GeneticCode* gCode,
  NucleotideSubstitutionModel* pmod,
  const AlphabetIndex2* pdist) :
  AbstractParameterAliasable("KronCodonDist."),
  AbstractKroneckerCodonSubstitutionModel(gCode, pmod, "KronCodonDist."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "KronCodonDist.")
{
  computeFrequencies(true);
  updateMatrices();
}

KroneckerCodonDistanceSubstitutionModel::KroneckerCodonDistanceSubstitutionModel(
  const GeneticCode* gCode,
  NucleotideSubstitutionModel* pmod1,
  NucleotideSubstitutionModel* pmod2,
  NucleotideSubstitutionModel* pmod3,
  const AlphabetIndex2* pdist) :
  AbstractParameterAliasable("KronCodonDist."),
  AbstractKroneckerCodonSubstitutionModel(gCode, pmod1, pmod2, pmod3, "KronCodonDist."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "KronCodonDist.")
{
  computeFrequencies(true);
  updateMatrices();
}

KroneckerCodonDistanceSubstitutionModel::KroneckerCodonDistanceSubstitutionModel(
  const GeneticCode* gCode,
  NucleotideSubstitutionModel* pmod,
  const std::vector<std::set< size_t> >& vPos,
  const AlphabetIndex2* pdist) :
  AbstractParameterAliasable("KronCodonDist."),
  AbstractKroneckerCodonSubstitutionModel(gCode, pmod, vPos, "KronCodonDist."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "KronCodonDist.")
{
  computeFrequencies(true);
  updateMatrices();
}

KroneckerCodonDistanceSubstitutionModel::KroneckerCodonDistanceSubstitutionModel(
  const GeneticCode* gCode,
  NucleotideSubstitutionModel* pmod1,
  NucleotideSubstitutionModel* pmod2,
  NucleotideSubstitutionModel* pmod3,
  const std::vector<std::set< size_t> >& vPos,
  const AlphabetIndex2* pdist) :
  AbstractParameterAliasable("KronCodonDist."),
  AbstractKroneckerCodonSubstitutionModel(gCode, pmod1, pmod2, pmod3, vPos, "KronCodonDist."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "KronCodonDist.")
{
  computeFrequencies(true);
  updateMatrices();
}
std::string KroneckerCodonDistanceSubstitutionModel::getName() const
{
  return "KronCodonDist";
}

void KroneckerCodonDistanceSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);

  // Beware: must be call at the end
  AbstractKroneckerCodonSubstitutionModel::fireParameterChanged(parameters);
}

double KroneckerCodonDistanceSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i, j)
         * AbstractKroneckerCodonSubstitutionModel::getCodonsMulRate(i, j);
}

void KroneckerCodonDistanceSubstitutionModel::setNamespace(const std::string& st)
{
  AbstractParameterAliasable::setNamespace(st);
  AbstractKroneckerCodonSubstitutionModel::setNamespace(st);
  AbstractCodonDistanceSubstitutionModel::setNamespace(st);
}
