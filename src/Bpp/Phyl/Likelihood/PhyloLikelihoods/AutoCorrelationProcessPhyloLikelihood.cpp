//
// File: AutoCorrelationProcessPhyloLikelihood.cpp
// Authors:
//   Laurent GuÃÂ©guen
// Created: lundi 23 septembre 2013, ÃÂ  22h 56
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

#include <Bpp/Numeric/Hmm/AutoCorrelationTransitionMatrix.h>

#include "AutoCorrelationProcessPhyloLikelihood.h"
#include "HmmLikelihood_DF.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

AutoCorrelationProcessPhyloLikelihood::AutoCorrelationProcessPhyloLikelihood(
  shared_ptr<const AlignmentDataInterface> data,
  shared_ptr<AutoCorrelationSequenceEvolution> processSeqEvol,
  shared_ptr<CollectionNodes> collNodes,
  size_t nSeqEvol,
  size_t nData) :
  AbstractPhyloLikelihood(collNodes->context()),
  AbstractAlignedPhyloLikelihood(collNodes->context(), data->getNumberOfSites()),
  AbstractSingleDataPhyloLikelihood(
      collNodes->context(),
      data->getNumberOfSites(),
      (processSeqEvol->getSubstitutionProcessNumbers().size() != 0)
          ? processSeqEvol->substitutionProcess(processSeqEvol->getSubstitutionProcessNumbers()[0]).getNumberOfStates()
	  : 0,
      nData),
  AbstractSequencePhyloLikelihood(collNodes->context(), processSeqEvol, nData),
  AbstractParametrizable(""),
  MultiProcessSequencePhyloLikelihood(data, processSeqEvol, collNodes, nSeqEvol, nData),
  Hpep_(),
  hmm_()
{
  auto alphyl = make_shared<HmmPhyloAlphabet>(*this);

  Hpep_ = make_shared<HmmPhyloEmissionProbabilities>(alphyl);

  hmm_ = shared_ptr<HmmLikelihood_DF>(new HmmLikelihood_DF(context(), processSeqEvol->getHmmProcessAlphabet(), processSeqEvol->getHmmTransitionMatrix(), Hpep_));

  addParameters_(hmm_->getParameters());
}

void AutoCorrelationProcessPhyloLikelihood::setNamespace(const std::string& nameSpace)
{
  MultiProcessSequencePhyloLikelihood::setNamespace(nameSpace);
  hmm_->setNamespace(nameSpace);
}


void AutoCorrelationProcessPhyloLikelihood::fireParameterChanged(const ParameterList& parameters)
{
  MultiProcessSequencePhyloLikelihood::fireParameterChanged(parameters);

  hmm_->matchParametersValues(parameters);
}
