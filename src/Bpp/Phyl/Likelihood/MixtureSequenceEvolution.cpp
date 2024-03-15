// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "MixtureSequenceEvolution.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

MixtureSequenceEvolution::MixtureSequenceEvolution(
  shared_ptr<SubstitutionProcessCollection> processColl,
  vector<size_t>& nProc) :
  MultiProcessSequenceEvolution(processColl, nProc, ""),
  simplex_(nProc.size(), 1, false, "Mixture.")
{
  // initialize parameters:
  addParameters_(simplex_.getParameters());
}

void MixtureSequenceEvolution::setSubProcessProb(const Simplex& si)
{
  simplex_.setFrequencies(si.getFrequencies());
  matchParametersValues(simplex_.getParameters());
}

void MixtureSequenceEvolution::setNamespace(const std::string& nameSpace)
{
  deleteParameters_(simplex_.getParameters().getParameterNames());

  simplex_.setNamespace(nameSpace);

  addParameters_(simplex_.getParameters());
}


void MixtureSequenceEvolution::fireParameterChanged(const ParameterList& parameters)
{
  MultiProcessSequenceEvolution::fireParameterChanged(parameters);
  simplex_.matchParametersValues(parameters);
}

ParameterList MixtureSequenceEvolution::getNonDerivableParameters() const
{
  ParameterList pl = MultiProcessSequenceEvolution::getNonDerivableParameters();
  pl.includeParameters(simplex_.getParameters());

  return pl;
}
