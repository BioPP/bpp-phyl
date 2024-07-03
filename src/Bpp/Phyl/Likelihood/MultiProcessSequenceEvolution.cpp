// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "MultiProcessSequenceEvolution.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

MultiProcessSequenceEvolution::MultiProcessSequenceEvolution(
    std::shared_ptr<SubstitutionProcessCollection> processColl,
    vector<size_t> nProc,
    const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  processColl_(processColl),
  nProc_(nProc)
{
  // initialize parameters:

  for (size_t i = 0; i < nProc_.size(); ++i)
  {
    includeParameters_(processColl_->getSubstitutionProcessParameters(nProc_[i], true));
  }
}

/******************************************************************************/

void MultiProcessSequenceEvolution::fireParameterChanged(const ParameterList& parameters)
{
  processColl_->matchParametersValues(parameters);
}

/******************************************************************************/

ParameterList MultiProcessSequenceEvolution::getSubstitutionProcessParameters(bool independent) const
{
  ParameterList pl;

  for (size_t i = 0; i < nProc_.size(); ++i)
  {
    pl.includeParameters(processColl_->getSubstitutionProcessParameters(nProc_[i], independent));
  }

  return pl;
}

ParameterList MultiProcessSequenceEvolution::getSubstitutionModelParameters(bool independent) const
{
  ParameterList pl;

  for (size_t i = 0; i < nProc_.size(); ++i)
  {
    pl.includeParameters(processColl_->substitutionProcess(nProc_[i]).getSubstitutionModelParameters(independent));
  }

  return pl;
}

ParameterList MultiProcessSequenceEvolution::getRateDistributionParameters(bool independent) const
{
  ParameterList pl;

  for (size_t i = 0; i < nProc_.size(); ++i)
  {
    pl.includeParameters(processColl_->substitutionProcess(nProc_[i]).getRateDistributionParameters(independent));
  }

  return pl;
}

ParameterList MultiProcessSequenceEvolution::getRootFrequenciesParameters(bool independent) const
{
  ParameterList pl;

  for (size_t i = 0; i < nProc_.size(); ++i)
  {
    pl.includeParameters(processColl_->substitutionProcess(nProc_[i]).getRootFrequenciesParameters(independent));
  }

  return pl;
}

ParameterList MultiProcessSequenceEvolution::getBranchLengthParameters(bool independent) const
{
  ParameterList pl;

  for (size_t i = 0; i < nProc_.size(); ++i)
  {
    pl.includeParameters(processColl_->substitutionProcess(nProc_[i]).getBranchLengthParameters(independent));
  }

  return pl;
}

ParameterList MultiProcessSequenceEvolution::getNonDerivableParameters() const
{
  ParameterList pl;

  for (size_t i = 0; i < nProc_.size(); ++i)
  {
    pl.includeParameters(processColl_->substitutionProcess(nProc_[i]).getNonDerivableParameters());
  }

  pl.includeParameters(getAliasedParameters(pl));

  return pl;
}

/******************************************************************************/


void MultiProcessSequenceEvolution::setParameters(const ParameterList& parameters)
{
  setParametersValues(parameters);
}


/*************************************************************************/

bool MultiProcessSequenceEvolution::isCompatibleWith(const AlignmentDataInterface& data) const
{
  for (size_t i = 0; i < nProc_.size(); ++i)
  {
    if (!processColl_->substitutionProcess(nProc_[i]).isCompatibleWith(data))
      return false;
  }

  return true;
}
