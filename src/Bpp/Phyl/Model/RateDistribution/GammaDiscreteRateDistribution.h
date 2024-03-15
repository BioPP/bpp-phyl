// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_RATEDISTRIBUTION_GAMMADISCRETERATEDISTRIBUTION_H
#define BPP_PHYL_MODEL_RATEDISTRIBUTION_GAMMADISCRETERATEDISTRIBUTION_H



// From bpp-core
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>

namespace bpp
{
class GammaDiscreteRateDistribution :
  public GammaDiscreteDistribution
{
public:
  GammaDiscreteRateDistribution(size_t nbClasses, double alpha = 1.) :
    GammaDiscreteDistribution(nbClasses, alpha, alpha)
  {
    aliasParameters("alpha", "beta");
  }

  GammaDiscreteRateDistribution* clone() const { return new GammaDiscreteRateDistribution(*this); }
};
} // end of namespace bpp;
#endif // BPP_PHYL_MODEL_RATEDISTRIBUTION_GAMMADISCRETERATEDISTRIBUTION_H
