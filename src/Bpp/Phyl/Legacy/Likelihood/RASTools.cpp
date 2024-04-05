// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>

#include "RASTools.h"

using namespace bpp;

// From the STL
#include <map>

using namespace std;

unique_ptr<DiscreteDistributionInterface> RASTools::getPosteriorRateDistribution(
    const DiscreteRatesAcrossSitesTreeLikelihoodInterface& treeLikelihood)
{
  // Get all posterior rate classes for each sites:
  vector<size_t> classes = treeLikelihood.getRateClassWithMaxPostProbPerSite();
  map<size_t, size_t> counts;
  for (size_t i = 0; i < classes.size(); i++)
  {
    counts[classes[i]]++;
  }

  // Now compute the distribution:
  auto rDist = treeLikelihood.getRateDistribution();
  map<double, double> distribution;
  for (auto& i : counts)
  {
    distribution[rDist->getCategory(i.first)] = static_cast<double>(i.second) / static_cast<double>(classes.size());
  }

  // Build a new distribution and return it:
  return make_unique<SimpleDiscreteDistribution>(distribution);
}
