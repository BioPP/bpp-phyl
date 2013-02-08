//
// File: CodonDistanceFitnessPhaseFrequenciesSubstitutionModel.cpp
// Created by: Fanny Pouyet 
// Created on: February 2012
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
  with loading,  using,  modifying and/or developi_ng or reproducing the
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


#include "CodonDistanceFitnessPhaseFrequenciesSubstitutionModel.h"
using namespace bpp;
using namespace std;

CodonDistanceFitnessPhaseFrequenciesSubstitutionModel::CodonDistanceFitnessPhaseFrequenciesSubstitutionModel(const GeneticCode* palph,
                                                             NucleotideSubstitutionModel* pmod,
                                                             FrequenciesSet* pfit,
                                                             FrequenciesSet* pfreq,
                                                             const AlphabetIndex2* pdist) :
  AbstractParameterAliasable("CodonDistFitPhasFreq."),
  AbstractSubstitutionModel(dynamic_cast<const CodonAlphabet*>(palph->getSourceAlphabet()), "CodonDistFitPhasFreq."),
  AbstractWordSubstitutionModel(dynamic_cast<const CodonAlphabet*>(palph->getSourceAlphabet()), "CodonDistFitPhasFreq."),
  AbstractCodonSubstitutionModel(dynamic_cast<const CodonAlphabet*>(palph->getSourceAlphabet()), pmod, "CodonDistFitPhasFreq."),
  AbstractCodonDistanceSubstitutionModel(palph, pdist, "CodonDistFitPhasFreq."),
  AbstractCodonPhaseFrequenciesSubstitutionModel(pfreq, "CodonDistFitPhasFreq."),
  AbstractCodonFitnessSubstitutionModel(pfit, "CodonDistFitPhasFreq.")
{
  updateMatrices();
}

CodonDistanceFitnessPhaseFrequenciesSubstitutionModel::CodonDistanceFitnessPhaseFrequenciesSubstitutionModel(const GeneticCode* palph,
                                                             NucleotideSubstitutionModel* pmod1,
                                                             NucleotideSubstitutionModel* pmod2,
                                                             NucleotideSubstitutionModel* pmod3,
                                                             FrequenciesSet* pfit,
                                                             FrequenciesSet* pfreq,
                                                             const AlphabetIndex2* pdist) :
  AbstractParameterAliasable("CodonDistFitPhasFreq."),
  AbstractSubstitutionModel(dynamic_cast<const CodonAlphabet*>(palph->getSourceAlphabet()), "CodonDistFitPhasFreq."),
  AbstractWordSubstitutionModel(dynamic_cast<const CodonAlphabet*>(palph->getSourceAlphabet()), "CodonDistFitPhasFreq."),
  AbstractCodonSubstitutionModel(dynamic_cast<const CodonAlphabet*>(palph->getSourceAlphabet()), pmod1, pmod2, pmod3, "CodonDistFitPhasFreq."),
  AbstractCodonDistanceSubstitutionModel(palph, pdist, "CodonDistFitPhasFreq."),
  AbstractCodonPhaseFrequenciesSubstitutionModel(pfreq, "CodonDistFitPhasFreq."),
  AbstractCodonFitnessSubstitutionModel(pfit,"CodonDistFitPhasFreq.")
{
  updateMatrices();
}

string CodonDistanceFitnessPhaseFrequenciesSubstitutionModel::getName() const
{
  return ("CodonDistFitPhasFreq");
}

void CodonDistanceFitnessPhaseFrequenciesSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);
  AbstractCodonPhaseFrequenciesSubstitutionModel::fireParameterChanged(parameters);
  AbstractCodonFitnessSubstitutionModel::fireParameterChanged(parameters);

  AbstractCodonSubstitutionModel::fireParameterChanged(parameters);
}

double CodonDistanceFitnessPhaseFrequenciesSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i,j)
    * AbstractCodonSubstitutionModel::getCodonsMulRate(i,j)
    * AbstractCodonPhaseFrequenciesSubstitutionModel::getCodonsMulRate(i,j)
    * AbstractCodonFitnessSubstitutionModel::getCodonsMulRate(i,j);
}

void CodonDistanceFitnessPhaseFrequenciesSubstitutionModel::setNamespace(const std::string& st)
{
  AbstractParameterAliasable::setNamespace(st);
  AbstractCodonSubstitutionModel::setNamespace(st);
  AbstractCodonDistanceSubstitutionModel::setNamespace(st);
  AbstractCodonPhaseFrequenciesSubstitutionModel::setNamespace(st); 
  AbstractCodonFitnessSubstitutionModel::setNamespace(st);
}

void CodonDistanceFitnessPhaseFrequenciesSubstitutionModel::setFreq(map<int,double>& frequencies)
{
  AbstractCodonPhaseFrequenciesSubstitutionModel::setFreq(frequencies);
}


