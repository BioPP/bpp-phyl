//
// File: ExtendedFloatVectorwiseOp.h
// Authors: Laurent Guéguen (2021)
// Created: samedi 17 avril 2021, à 08h 11
// Last modified: 2018-07-11
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef EXTENDEDFLOATVECTORWISE_OP_H
#define EXTENDEDFLOATVECTORWISE_OP_H

//#include "ExtendedFloatEigen.h"

/*
 * Class of operators to perform columnwise & rowwise operations on ExtendedFloatMatrix objects.
 *
 * This code is strongly influenced from Eigen code.
 */


namespace bpp {
  
  template<typename Derived> class ExtendedFloatEigenBase;

  template<int R, int C, template<int R2, int C2> class EigenType> class ExtendedFloatEigen;

  template<typename DerivedEF, typename MatType, int Direction>
  class ExtendedFloatVectorwiseOp {
  private:
    ExtendedFloatEigenBase<DerivedEF>& efMat_;
    
    Eigen::VectorwiseOp<MatType, Direction> eigenVWiseOp_;

  public:
    ExtendedFloatVectorwiseOp(ExtendedFloatEigenBase<DerivedEF>& der) :
      efMat_(der),
      eigenVWiseOp_(der.derived().float_part()) {}

    template<typename OtherDerived>
    ExtendedFloatEigenBase<DerivedEF>& operator=(const ExtendedFloatEigenBase<OtherDerived>& otherDerived)
    {
      eigenVWiseOp_=otherDerived.derived().float_part();
      efMat_.derived().exponent_part()=otherDerived.derived().exponent_part();
      
      return efMat_;
    }

  };

  
}

#endif // EXTENDEDFLOATVECTORWISE_OP_H
