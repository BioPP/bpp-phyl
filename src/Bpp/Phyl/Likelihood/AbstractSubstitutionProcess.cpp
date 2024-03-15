// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractSubstitutionProcess.h"

using namespace bpp;
using namespace std;

ParameterList AbstractSubstitutionProcess::getNonDerivableParameters() const
{
  ParameterList pl = getSubstitutionModelParameters(true);
  pl.includeParameters(getRootFrequenciesParameters(true));
  pl.includeParameters(getRateDistributionParameters(true));

  pl.includeParameters(getAliasedParameters(pl));
  pl.includeParameters(getFromParameters(pl));

  return pl;
}

