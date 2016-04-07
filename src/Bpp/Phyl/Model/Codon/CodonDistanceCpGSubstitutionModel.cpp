//
// File: CodonDistanceSubstitutionModel.cpp
// Created by:  Laurent Gueguen
// Created on: Feb 2009
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

#include "CodonDistanceCpGSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

CodonDistanceCpGSubstitutionModel::CodonDistanceCpGSubstitutionModel(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod,
    const AlphabetIndex2* pdist) :
  AbstractParameterAliasable("CodonDistCpG."),
  AbstractCodonSubstitutionModel(gCode, pmod, "CodonDistCpG."),
  AbstractCodonDistanceSubstitutionModel(pdist, "CodonDistCpG."),
  AbstractCodonCpGSubstitutionModel("CodonDistCpG.")
{
  updateMatrices();
}

CodonDistanceCpGSubstitutionModel::CodonDistanceCpGSubstitutionModel(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod1,
    NucleotideSubstitutionModel* pmod2,
    NucleotideSubstitutionModel* pmod3,
    const AlphabetIndex2* pdist) :
  AbstractParameterAliasable("CodonDistCpG."),
  AbstractCodonSubstitutionModel(gCode, pmod1, pmod2, pmod3, "CodonDistCpG."),
  AbstractCodonDistanceSubstitutionModel(pdist, "CodonDistCpG."),
  AbstractCodonCpGSubstitutionModel("CodonDistCpG.")
{
  updateMatrices();
}

std::string CodonDistanceCpGSubstitutionModel::getName() const
{
  return "CodonDistCpG";
}

void CodonDistanceCpGSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);
  AbstractCodonCpGSubstitutionModel::fireParameterChanged(parameters);

  // Beware: must be call at the end
  AbstractCodonSubstitutionModel::fireParameterChanged(parameters);
}

double CodonDistanceCpGSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i,j)
    * AbstractCodonCpGSubstitutionModel::getCodonsMulRate(i,j)
    * AbstractCodonSubstitutionModel::getCodonsMulRate(i,j);
}

