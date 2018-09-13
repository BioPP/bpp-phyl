//
// File: AbstractCodonAAFitnessSubstitutionModel.cpp
// Created by: Laurent Guéguen
// Created on: mercredi 8 novembre 2017, à 21h 10
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

#include "AbstractCodonAAFitnessSubstitutionModel.h"

#include <Bpp/Seq/Alphabet/AlphabetTools.h>

using namespace bpp;
using namespace std;

/****************************************************************************************/

AbstractCodonAAFitnessSubstitutionModel::AbstractCodonAAFitnessSubstitutionModel(FrequenciesSet* pfitset, const GeneticCode* pgencode, const string& prefix):
  AbstractParameterAliasable(prefix), pfitset_(pfitset), pgencode_(pgencode), fitName_(""), stateMap_(new CanonicalStateMap(pgencode->getSourceAlphabet(), false)),
  protStateMap_(&pfitset->getStateMap()), Ns_(1)
{
  if (!AlphabetTools::isProteicAlphabet(pfitset_->getAlphabet()))
    throw Exception("AbstractCodonAAFitnessSubstitutionModel::AbstractCodonAAFitnessSubstitutionModel need Proteic Fitness.");
  
  fitName_ = "fit_" + pfitset_->getNamespace();
  pfitset_->setNamespace(prefix + fitName_);
  
  addParameters_(pfitset_->getParameters());
}

AbstractCodonAAFitnessSubstitutionModel::~AbstractCodonAAFitnessSubstitutionModel()
{
}

void AbstractCodonAAFitnessSubstitutionModel::fireParameterChanged (const ParameterList& parameters)
{
  AbstractParameterAliasable::fireParameterChanged(parameters); 
  if (hasParameter("Ns"))
    Ns_=getParameterValue("Ns");
  
  pfitset_->matchParametersValues(parameters);
}

void AbstractCodonAAFitnessSubstitutionModel::setFreq(map<int, double>& frequencies)
{
  pfitset_->setFrequenciesFromAlphabetStatesFrequencies(frequencies);
  matchParametersValues(pfitset_->getParameters() );
}

double AbstractCodonAAFitnessSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  double mu;

  int aai = pgencode_->translate(stateMap_->getAlphabetStateAsInt(i));
  int aaj = pgencode_->translate(stateMap_->getAlphabetStateAsInt(j));
  
  double phi_i= pfitset_->getFrequencies() [protStateMap_->getModelStates(aai)[0]];
  double phi_j= pfitset_->getFrequencies() [protStateMap_->getModelStates(aaj)[0]];

  if (phi_i == phi_j || Ns_==0)
    mu=1;
  else if (phi_i==0)
    mu=100;
  else if (phi_j==0)
    mu=0;
  else
  {
    if (Ns_==1)
      mu = -(log(phi_i/phi_j)/(1-(phi_i/phi_j)));
    else
    {
      double x=pow(phi_i/phi_j,Ns_);
      mu = -log(x)/(1-x);
    }
  }
  
  return mu;
}

