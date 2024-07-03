// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "OneProcessSequenceEvolution.h"

using namespace bpp;
using namespace std;

OneProcessSequenceEvolution::OneProcessSequenceEvolution(std::shared_ptr<SubstitutionProcessInterface> process, size_t nProc) :
  AbstractParameterAliasable(""),
  subsProc_(process),
  nProc_(nProc),
  vProc_(std::vector<size_t>(1, nProc_))
{
  includeParameters_(process->getSubstitutionModelParameters(true));
  includeParameters_(process->getRateDistributionParameters(true));
  includeParameters_(process->getRootFrequenciesParameters(true));
  includeParameters_(process->getBranchLengthParameters(true));
}

OneProcessSequenceEvolution::OneProcessSequenceEvolution(const OneProcessSequenceEvolution& evol) :
  AbstractParameterAliasable(evol),
  subsProc_(evol.subsProc_),
  nProc_(evol.nProc_),
  vProc_(evol.vProc_)
{}

OneProcessSequenceEvolution& OneProcessSequenceEvolution::operator=(const OneProcessSequenceEvolution& evol)
{
  AbstractParameterAliasable::operator=(evol);
  subsProc_ = evol.subsProc_;
  nProc_ = evol.nProc_;
  vProc_ = evol.vProc_;

  return *this;
}
