// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_RATEDISTRIBUTION_GAUSSIANDISCRETERATEDISTRIBUTION_H
#define BPP_PHYL_MODEL_RATEDISTRIBUTION_GAUSSIANDISCRETERATEDISTRIBUTION_H


// From bpp-core
#include <Bpp/Numeric/Prob/GaussianDiscreteDistribution.h>

namespace bpp
{
class GaussianDiscreteRateDistribution :
  public GaussianDiscreteDistribution
{
public:
  GaussianDiscreteRateDistribution(size_t nbClasses, double sigma) :
    GaussianDiscreteDistribution(nbClasses, 1., sigma)
  {
    deleteParameter_(0);
  }

  GaussianDiscreteRateDistribution* clone() const { return new GaussianDiscreteRateDistribution(*this); }
};
} // end of namespace bpp;
#endif // BPP_PHYL_MODEL_RATEDISTRIBUTION_GAUSSIANDISCRETERATEDISTRIBUTION_H
