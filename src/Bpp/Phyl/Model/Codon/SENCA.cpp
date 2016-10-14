//
// File: SENCA.cpp
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


#include "SENCA.h"
using namespace bpp;
using namespace std;

SENCA::SENCA(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod,
    FrequenciesSet* pfit,
    const AlphabetIndex2* pdist) :
  AbstractParameterAliasable("SENCA."),
  AbstractCodonSubstitutionModel(gCode, pmod, "SENCA."),
  AbstractCodonDistanceSubstitutionModel(pdist, "SENCA."),
  AbstractCodonFitnessSubstitutionModel(pfit, "SENCA.")
{
  updateMatrices();
}

SENCA::SENCA(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod1,
    NucleotideSubstitutionModel* pmod2,
    NucleotideSubstitutionModel* pmod3,
    FrequenciesSet* pfit,
    const AlphabetIndex2* pdist) :
  AbstractParameterAliasable("SENCA."),
  AbstractCodonSubstitutionModel(gCode, pmod1, pmod2, pmod3, "SENCA."),
  AbstractCodonDistanceSubstitutionModel(pdist, "SENCA."),
  AbstractCodonFitnessSubstitutionModel(pfit,"SENCA.")
{
  updateMatrices();
}

string SENCA::getName() const
{
  return ("SENCA");
}

void SENCA::fireParameterChanged(const ParameterList& parameters)
{
  AbstractCodonDistanceSubstitutionModel::fireParameterChanged(parameters);
  AbstractCodonFitnessSubstitutionModel::fireParameterChanged(parameters);

  AbstractCodonSubstitutionModel::fireParameterChanged(parameters);
}

double SENCA::getCodonsMulRate(size_t i, size_t j) const
{
  return AbstractCodonDistanceSubstitutionModel::getCodonsMulRate(i,j)
    * AbstractCodonSubstitutionModel::getCodonsMulRate(i,j)
    * AbstractCodonFitnessSubstitutionModel::getCodonsMulRate(i,j);
}

void SENCA::setNamespace(const std::string& st)
{
  AbstractParameterAliasable::setNamespace(st);
  AbstractCodonSubstitutionModel::setNamespace(st);
  AbstractCodonDistanceSubstitutionModel::setNamespace(st);
  AbstractCodonFitnessSubstitutionModel::setNamespace(st);
}

void SENCA::setFreq(map<int,double>& frequencies)
{ 
   AbstractCodonSubstitutionModel::setFreq(frequencies);
   const Vdouble& freq1 = AbstractCodonSubstitutionModel::getFrequencies();
   const Alphabet* alphabet = getAlphabet();

   map<int, double> freq2;
   double s=0;
   map<int, double>::iterator it;

    for (it=frequencies.begin();it!=frequencies.end();it++)
    {
      freq2[it->first]=(freq1[alphabet->getStateIndex(it->first)-1] != 0 ? it->second/freq1[alphabet->getStateIndex(it->first)-1] : 0);
      s += freq2[it->first];
    }
  
    for (it = freq2.begin(); it != freq2.end(); it++)
      freq2[it->first] /= s;

    AbstractCodonFitnessSubstitutionModel::setFreq(freq2);

  updateMatrices();
}


