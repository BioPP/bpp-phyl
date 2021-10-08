//
// File: ExtendedFloatEigenTools.h
// Authors:
//   Laurent GuÃ©guen (2021)
// Created: samedi 17 avril 2021, Ã  08h 11
// Last modified: 2018-07-11 00:00:00
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_NEWLIKELIHOOD_DATAFLOW_EXTENDEDFLOATEIGENTOOLS_H
#define BPP_PHYL_NEWLIKELIHOOD_DATAFLOW_EXTENDEDFLOATEIGENTOOLS_H


#include "ExtendedFloat.h"

/*
 * Class of tool operators  ExtendedFloatMatrix objects, needed from Eigen tools.
 *
 * This code is strongly influenced from Eigen code.
 */


namespace bpp
{
template<typename Derived> class ExtendedFloatEigenBase;

template<int R, int C, template<int R2, int C2> class EigenType> class ExtendedFloatEigen;


/*
 * Class of operators to perform columnwise & rowwise operations on ExtendedFloatMatrix objects.
 *
 * This code is strongly influenced from Eigen code.
 */

template<typename DerivedEF, typename MatType, int Direction>
class ExtendedFloatVectorwiseOp
{
private:
  DerivedEF& efMat_;

  Eigen::VectorwiseOp<MatType, Direction> eigenVWiseOp_;

public:
  ExtendedFloatVectorwiseOp(DerivedEF& der) :
    efMat_(der),
    eigenVWiseOp_(der.float_part()) {}

  template<typename OtherDerived>
  DerivedEF& operator=(const OtherDerived& otherDerived)
  {
    eigenVWiseOp_ = otherDerived.float_part();
    efMat_.exponent_part() = otherDerived.exponent_part();

    return efMat_;
  }

  DerivedEF sum() const
  {
    auto s = eigenVWiseOp_.sum();
    return DerivedEF(s, efMat_.derived().exponent_part());
  }

  DerivedEF mean() const
  {
    auto s = eigenVWiseOp_.mean();
    return DerivedEF(s, efMat_.derived().exponent_part());
  }
};


/*
 * Manage NoAlias for Eigen::Matrix objects in ExtendedFloatMatrix;
 *
 */

template<typename DerivedEF>
class ExtendedFloatNoAlias
{
protected:
  ExtendedFloatEigenBase<DerivedEF>& efMat_;

public:
  explicit ExtendedFloatNoAlias(ExtendedFloatEigenBase<DerivedEF>& der) :
    efMat_(der) {}

  template<typename Otherderived>
  ExtendedFloatEigenBase<DerivedEF>& operator=(const ExtendedFloatEigenBase<Otherderived>& other)
  {
    efMat_.derived().float_part().noalias() = other.derived().float_part();
    efMat_.derived().exponent_part() = other.derived().exponent_part();
    return efMat_;
  }

  template<typename Otherderived>
  ExtendedFloatEigenBase<DerivedEF>& operator+=(const ExtendedFloatEigenBase<Otherderived>& other)
  {
    auto& thisnoal(efMat_.derived().float_part().noalias());
    const auto& rhs(other.derived());

    if (efMat_.derived().exponent_part () >= rhs.exponent_part ())
    {
      thisnoal += rhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, rhs.exponent_part () - efMat_().derived().exponent_part ());
    }
    else
    {
      thisnoal = rhs.float_part () + efMat_.derived().float_part () * constexpr_power<double>(ExtendedFloat::radix, efMat_.derived().exponent_part () - rhs.exponent_part ());
      efMat_.derived().exponent_part() = rhs.exponent_part ();
    }
    efMat_.derived().normalize ();
    return efMat_;
  }

  template<typename Otherderived>
  ExtendedFloatEigenBase<DerivedEF>& operator-=(const ExtendedFloatEigenBase<Otherderived>& other)
  {
    auto& thisnoal(efMat_.derived().float_part().noalias());
    const auto& rhs(other.derived());

    if (efMat_.derived().exponent_part () >= rhs.exponent_part ())
    {
      thisnoal -= rhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, rhs.exponent_part () - efMat_().derived().exponent_part ());
    }
    else
    {
      thisnoal = rhs.float_part () - efMat_.derived().float_part () * constexpr_power<double>(ExtendedFloat::radix, efMat_.derived().exponent_part () - rhs.exponent_part ());
      efMat_.derived().exponent_part() = rhs.exponent_part ();
    }
    efMat_.derived().normalize ();
    return efMat_;
  }
};
}
#endif // BPP_PHYL_NEWLIKELIHOOD_DATAFLOW_EXTENDEDFLOATEIGENTOOLS_H
