//
// File: ExponentialDiscreteRateDistribution.h
// Created by: Julien Dutheil
// Created on: Fri Nov 16 14:37 2012
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2012)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef _EXPONENTIALDISCRETERATEDISTRIBUTION_H_
#define _EXPONENTIALDISCRETERATEDISTRIBUTION_H_

//From bpp-core
#include <Bpp/Numeric/Prob/ExponentialDiscreteDistribution.h>

namespace bpp {

class ExponentialDiscreteRateDistribution:
  public ExponentialDiscreteDistribution
{
  public:
    ExponentialDiscreteRateDistribution(size_t nbClasses):
      AbstractParameterAliasable("Exponential."),
      ExponentialDiscreteDistribution(nbClasses, 1.)
    {
      deleteParameter_(0);
    }

    ExponentialDiscreteRateDistribution* clone() const { return new ExponentialDiscreteRateDistribution(*this); }
    
};

} //end of namespace bpp;

#endif //_EXPONENTIALDISCRETERATEDISTRIBUTION_H_

