//
// File: DFPDistanceFrequenciesSubstitutionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: mercredi 4 novembre 2020, ÃÂ  13h 19
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
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


#include "DFPDistanceFrequenciesSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

DFPDistanceFrequenciesSubstitutionModel::DFPDistanceFrequenciesSubstitutionModel(
    shared_ptr<const GeneticCode> gCode,
    unique_ptr<CodonFrequencySetInterface> pfreq,
    shared_ptr<const AlphabetIndex2> pdist) :
  AbstractParameterAliasable("DFPDistFreq."),
  AbstractDFPSubstitutionModel(gCode, "DFPDistFreq."),
  AbstractCodonDistanceSubstitutionModel(pdist, gCode, "DFPDistFreq."),
  AbstractCodonFrequenciesSubstitutionModel(std::move(pfreq), "DFPDistFreq.")
{
  updateMatrices_();
}

std::string DFPDistanceFrequenciesSubstitutionModel::getName() const
{
  return "DFPDistFreq";
}

void DFPDistanceFrequenciesSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);
  AbstractCodonFrequenciesSubstitutionModel::fireParameterChanged(parameters);

  // Beware: must be called last, since it calls updateMatrices
  AbstractDFPSubstitutionModel::fireParameterChanged(parameters);
}

double DFPDistanceFrequenciesSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i, j)
         * AbstractDFPSubstitutionModel::getCodonsMulRate(i, j)
         * AbstractCodonFrequenciesSubstitutionModel::getCodonsMulRate(i, j);
}

void DFPDistanceFrequenciesSubstitutionModel::setNamespace(const std::string& st)
{
  AbstractParameterAliasable::setNamespace(st);
  AbstractDFPSubstitutionModel::setNamespace(st);
  AbstractCodonDistanceSubstitutionModel::setNamespace(st);
  AbstractCodonFrequenciesSubstitutionModel::setNamespace(st);
}

void DFPDistanceFrequenciesSubstitutionModel::setFreq(map<int, double>& frequencies)
{
  AbstractCodonFrequenciesSubstitutionModel::setFreq(frequencies);
  getFrequencies_() = AbstractCodonFrequenciesSubstitutionModel::codonFrequencySet().getFrequencies();
}
