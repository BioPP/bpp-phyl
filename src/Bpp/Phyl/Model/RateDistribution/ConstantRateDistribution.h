// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_RATEDISTRIBUTION_CONSTANTRATEDISTRIBUTION_H
#define BPP_PHYL_MODEL_RATEDISTRIBUTION_CONSTANTRATEDISTRIBUTION_H


// From bpp-core
#include <Bpp/Numeric/Prob/ConstantDistribution.h>

namespace bpp
{
class ConstantRateDistribution :
  public ConstantDistribution
{
public:
  ConstantRateDistribution(double rate = 1) :
    ConstantDistribution(rate)
  {
    deleteParameter_(0);
  }

  ConstantRateDistribution* clone() const { return new ConstantRateDistribution(*this); }
};
} // end of namespace bpp;
#endif // BPP_PHYL_MODEL_RATEDISTRIBUTION_CONSTANTRATEDISTRIBUTION_H
