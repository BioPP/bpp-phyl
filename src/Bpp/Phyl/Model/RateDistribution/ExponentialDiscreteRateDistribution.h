// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_RATEDISTRIBUTION_EXPONENTIALDISCRETERATEDISTRIBUTION_H
#define BPP_PHYL_MODEL_RATEDISTRIBUTION_EXPONENTIALDISCRETERATEDISTRIBUTION_H



// From bpp-core
#include <Bpp/Numeric/Prob/ExponentialDiscreteDistribution.h>

namespace bpp
{
class ExponentialDiscreteRateDistribution :
  public ExponentialDiscreteDistribution
{
public:
  ExponentialDiscreteRateDistribution(size_t nbClasses) :
    ExponentialDiscreteDistribution(nbClasses, 1.)
  {
    deleteParameter_(0);
  }

  ExponentialDiscreteRateDistribution* clone() const { return new ExponentialDiscreteRateDistribution(*this); }
};
} // end of namespace bpp;
#endif // BPP_PHYL_MODEL_RATEDISTRIBUTION_EXPONENTIALDISCRETERATEDISTRIBUTION_H
