//
// File: AbstractCodonAARateSubstitutionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: lundi 30 octobre 2017, ÃÂ  06h 07
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

#include <Bpp/Numeric/NumConstants.h>

#include "AbstractCodonAARateSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractCodonAARateSubstitutionModel::AbstractCodonAARateSubstitutionModel(
  shared_ptr<ProteinSubstitutionModelInterface> pmodel,
  shared_ptr<const GeneticCode> pgencode,
  const string& prefix,
  bool paramSynRate) :
  AbstractParameterAliasable(prefix),
  pAAmodel_(pmodel),
  pgencode_(pgencode),
  beta_(19),
  gamma_(1),
  stateMap_(new CanonicalStateMap(pgencode->getSourceAlphabet(), false))
{
  if (paramSynRate)
    addParameter_(new Parameter(prefix + "gamma", 1, std::make_shared<IntervalConstraint>(NumConstants::SMALL(), 999, true, true)));

  addParameter_(new Parameter(prefix + "beta", 1, std::make_shared<IntervalConstraint>(NumConstants::SMALL(), 999, true, true)));

  pAAmodel_->enableEigenDecomposition(false);

  pAAmodel_->setNamespace(prefix + pAAmodel_->getNamespace());
  addParameters_(pAAmodel_->getParameters());
}

void AbstractCodonAARateSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  pAAmodel_->matchParametersValues(parameters);

  if (hasParameter("gamma"))
    gamma_ = getParameterValue("gamma");

  beta_ = getParameterValue("beta");
}

double AbstractCodonAARateSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  int si(stateMap_->getAlphabetStateAsInt(i)), sj(stateMap_->getAlphabetStateAsInt(j));

  return pgencode_->areSynonymous(si, sj) ? gamma_ :
         beta_* pAAmodel_->Qij(pAAmodel_->getModelStates(pgencode_->translate(si))[0],
                               pAAmodel_->getModelStates(pgencode_->translate(sj))[0]);
}
