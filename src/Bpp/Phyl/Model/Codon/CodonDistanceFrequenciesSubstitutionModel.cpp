//
// File: CodonDistanceFrequenciesSubstitutionModel.cpp
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

#include "CodonDistanceFrequenciesSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

CodonDistanceFrequenciesSubstitutionModel::CodonDistanceFrequenciesSubstitutionModel(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod,
    FrequenciesSet* pfreq,
    const AlphabetIndex2* pdist,
    bool paramSynRate) :
  AbstractParameterAliasable("CodonDistFreq."),
  AbstractCodonSubstitutionModel(gCode, pmod, "CodonDistFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "CodonDistFreq.", paramSynRate),
  AbstractCodonFrequenciesSubstitutionModel(pfreq, "CodonDistFreq.")
{
  computeFrequencies(true); // for initialization
  updateMatrices();
  computeFrequencies(false);
}

CodonDistanceFrequenciesSubstitutionModel::CodonDistanceFrequenciesSubstitutionModel(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod1,
    NucleotideSubstitutionModel* pmod2,
    NucleotideSubstitutionModel* pmod3,
    FrequenciesSet* pfreq,
    const AlphabetIndex2* pdist,
    bool paramSynRate) :
  AbstractParameterAliasable("CodonDistFreq."),
  AbstractCodonSubstitutionModel(gCode, pmod1, pmod2, pmod3, "CodonDistFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "CodonDistFreq.", paramSynRate),
  AbstractCodonFrequenciesSubstitutionModel(pfreq, "CodonDistFreq.")
{
  computeFrequencies(true); // for initialization
  updateMatrices();
  computeFrequencies(false);
}

std::string CodonDistanceFrequenciesSubstitutionModel::getName() const
{
  return "CodonDistFreq";
}

void CodonDistanceFrequenciesSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);
  AbstractCodonFrequenciesSubstitutionModel::fireParameterChanged(parameters);
  getFrequencies_()=AbstractCodonFrequenciesSubstitutionModel::getFrequenciesSet()->getFrequencies();

  // Beware: must be call at the end
  AbstractCodonSubstitutionModel::fireParameterChanged(parameters);
}

double CodonDistanceFrequenciesSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i,j)
    * AbstractCodonSubstitutionModel::getCodonsMulRate(i,j)
    * AbstractCodonFrequenciesSubstitutionModel::getCodonsMulRate(i,j);
}

void CodonDistanceFrequenciesSubstitutionModel::setNamespace(const std::string& st)
{
  AbstractParameterAliasable::setNamespace(st);
  AbstractCodonSubstitutionModel::setNamespace(st);
  AbstractCodonDistanceSubstitutionModel::setNamespace(st);
  AbstractCodonFrequenciesSubstitutionModel::setNamespace(st);
}

void CodonDistanceFrequenciesSubstitutionModel::setFreq(map<int,double>& frequencies)
{
  AbstractCodonFrequenciesSubstitutionModel::setFreq(frequencies);
  getFrequencies_()=AbstractCodonFrequenciesSubstitutionModel::getFrequenciesSet()->getFrequencies();
}
