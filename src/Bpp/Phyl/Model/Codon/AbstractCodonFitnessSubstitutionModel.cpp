//
// File: AbstractCodonFitnessSubstitutionModel.cpp
// Created by:  Fanny Pouyet
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

# include "AbstractCodonFitnessSubstitutionModel.h"
using namespace bpp;
using namespace std;
/****************************************************************************************/
AbstractCodonFitnessSubstitutionModel::AbstractCodonFitnessSubstitutionModel(FrequenciesSet* pfitset, const string& prefix):
  CodonSubstitutionModel(), AbstractParameterAliasable(prefix), pfitset_(pfitset), fitName_("")
{
  if (dynamic_cast<CodonFrequenciesSet*>(pfitset) == NULL)
    throw Exception ("Bad type for fitness parameters"+ pfitset ->getName() );
  fitName_="fit_"+ pfitset_->getNamespace();
  pfitset_->setNamespace(prefix + fitName_);
  addParameters_(pfitset_->getParameters() );
}

AbstractCodonFitnessSubstitutionModel::~AbstractCodonFitnessSubstitutionModel()
{
  if (pfitset_) delete pfitset_;
}

void AbstractCodonFitnessSubstitutionModel::fireParameterChanged (const ParameterList& parameters)
{
  pfitset_->matchParametersValues(parameters);
}

void AbstractCodonFitnessSubstitutionModel::setFreq(map<int, double>& frequencies)
{
  pfitset_->setFrequenciesFromAlphabetStatesFrequencies(frequencies);
  matchParametersValues(pfitset_->getParameters() );
}

double AbstractCodonFitnessSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  double mu;
  double phi_j= pfitset_->getFrequencies() [j];
  double phi_i= pfitset_->getFrequencies() [i];
  if (phi_i == phi_j) mu=1;
  else if (phi_i==0)
    mu=100;
  else if (phi_j==0)
    mu=0;
  else
    mu = -(log(phi_i/phi_j)/(1-(phi_i/phi_j)));
  return mu;
}

