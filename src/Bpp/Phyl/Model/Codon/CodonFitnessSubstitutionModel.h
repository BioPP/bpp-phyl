//AbstractCodonBiasUse


//
// File: CodonFitnessSubstitutionModel.h
// Created by: Fanny Pouyet 
// Created on: mars 2012
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
 
#ifndef _CODONFITNESSSUBSTITUTIONMODEL_H_
#define _CODONFITNESSSUBSTITUTIONMODEL_H_
 
#include "AbstractCodonFitnessSubstitutionModel.h"
#include "AbstractCodonSubstitutionModel.h"
#include "AbstractCodonDistanceSubstitutionModel.h"
#include "AbstractCodonPhaseFrequenciesSubstitutionModel.h"

namespace bpp
{

  class CodonFitnessSubstitutionModel :
    public AbstractCodonSubstitutionModel,
    public AbstractCodonDistanceSubstitutionModel,
    public AbstractCodonPhaseFrequenciesSubstitutionModel,
    public AbstractCodonFitnessSubstitutionModel
  {
  public:
    CodonFitnessSubstitutionModel(const GeneticCode* palph,
                                  NucleotideSubstitutionModel* pmod,
                                  FrequenciesSet* pfit,
                                  FrequenciesSet* pfreq,
                                  const AlphabetIndex2<double>* pdist = 0);
    CodonFitnessSubstitutionModel(const GeneticCode* palph,
                                  NucleotideSubstitutionModel* pmod1,
                                  NucleotideSubstitutionModel* pmod2,
                                  NucleotideSubstitutionModel* pmod3,
                                  FrequenciesSet* pfit,
                                  FrequenciesSet* pfreq,
                                  const AlphabetIndex2<double>* pdist = 0);

    CodonFitnessSubstitutionModel(const CodonFitnessSubstitutionModel& model) :
      AbstractParameterAliasable(model),
      AbstractSubstitutionModel(model),
      AbstractWordSubstitutionModel(model),
      AbstractCodonSubstitutionModel(model),
      AbstractCodonDistanceSubstitutionModel(model),
      AbstractCodonPhaseFrequenciesSubstitutionModel(model),
      AbstractCodonFitnessSubstitutionModel(model)
    {}

    CodonFitnessSubstitutionModel& operator=(
                                             const CodonFitnessSubstitutionModel& model)
    {
      AbstractParameterAliasable::operator=(model);
      AbstractSubstitutionModel::operator=(model);
      AbstractWordSubstitutionModel::operator=(model);
      AbstractCodonSubstitutionModel::operator=(model);
      AbstractCodonDistanceSubstitutionModel::operator=(model);
      AbstractCodonPhaseFrequenciesSubstitutionModel::operator=(model);
      AbstractCodonFitnessSubstitutionModel::operator=(model);
      return *this;
    }

    ~CodonFitnessSubstitutionModel() {}

    CodonFitnessSubstitutionModel* clone() const
    {
      return new CodonFitnessSubstitutionModel(*this);
    }

  public:
    void fireParameterChanged(const ParameterList& parameterlist);

    std::string getName() const;

    double getCodonsMulRate(unsigned int, unsigned int) const;

    void setNamespace(const std::string&);

    void setFreq(std::map<int,double>& frequencies);

  };

} // end of namespace bpp.

#endif

