// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_EXTENDEDFLOATEIGENTOOLS_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_EXTENDEDFLOATEIGENTOOLS_H


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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++" // Remove EIGEN warning
  template<typename OtherDerived>
  DerivedEF& operator=(const OtherDerived& otherDerived)
  {
    eigenVWiseOp_ = otherDerived.float_part();
    efMat_.exponent_part() = otherDerived.exponent_part();

    return efMat_;
  }
#pragma GCC diagnostic pop

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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++" // Remove EIGEN warning
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
#pragma GCC diagnostic pop
};


/*
 * Manage Row in a ExtendedFloatMatrix;
 *
 */

template< int R,  int C, template< int R2,  int C2> class EigenType>
class ExtendedFloatRow
{
  using VecType = EigenType<1, C>;
  using ExtType =  int;

protected:
  ExtendedFloatEigen<R, C, EigenType>& efMat_;
  Eigen::Index nrow_;

private:
  VecType tmp_;

public:
  ExtendedFloatRow(ExtendedFloatEigen<R, C, EigenType>& der, Eigen::Index nrow) :
    efMat_(der), nrow_(nrow), tmp_(efMat_.float_part().row(nrow_)) {}

  const ExtType& exponent_part () const { return efMat_.exponent_part(); }

  const VecType& float_part () const { return tmp_; }

  ExtType& exponent_part () { return efMat_.exponent_part(); }

  VecType& float_part () { return tmp_; }

  /*
   * @brief Includes a row in an ExtendedFloatMatrix. Exponent parts
   * of both are to be fit, as in denorm_add method.
   *
   */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++" // Remove a EIGEN warning
  ExtendedFloatRow& operator=(const ExtendedFloatEigen<1, C, EigenType>& row)
  {
    efMat_.float_part().row(nrow_) = row.float_part() * constexpr_power<double>(ExtendedFloat::radix, row.exponent_part() - efMat_.exponent_part());
    efMat_.normalize();
    return *this;
  }
#pragma GCC diagnostic pop
};

/*
 * Manage Col in a ExtendedFloatMatrix;
 *
 */

template< int R,  int C, template< int R2,  int C2> class EigenType>
class ExtendedFloatCol : public ExtendedFloatEigen<R, 1, EigenType>
{
  using VecType = EigenType<R, 1>;
  using ExtType =  int;

protected:
  ExtendedFloatEigen<R, C, EigenType>& efMat_;
  Eigen::Index ncol_;

private:
  VecType tmp_;

public:
  ExtendedFloatCol(ExtendedFloatEigen<R, C, EigenType>& der, Eigen::Index ncol) :
    efMat_(der), ncol_(ncol), tmp_(efMat_.float_part().col(ncol_)) {}

  /*
   * @brief Includes a row in an ExtendedFloatMatrix. Exponent parts
   * of both are to be fit, as in denorm_add method.
   *
   */
  const ExtType& exponent_part () const
  {
    return efMat_.exponent_part();
  }

  const VecType& float_part () const
  {
    return tmp_;
  }

  ExtType& exponent_part ()
  {
    return efMat_.exponent_part();
  }

  VecType& float_part ()
  {
    return tmp_;
  }

  ExtendedFloatCol& operator=(const ExtendedFloatEigen<R, 1, EigenType>& col)
  {
    efMat_.float_part().col(ncol_) = col.float_part() * constexpr_power<double>(ExtendedFloat::radix, col.exponent_part() - efMat_.exponent_part());
    efMat_.normalize();
    return *this;
  }

  // std::ostream& operator<<(std::ostream& out, const Self& ef)
  // {
  //   return out << "( " << efMat_.float_part() << " ) *  2^" << efMat_.exponent_part();
  // }

  // const ExtendedFloat& sum() const
  // {
  //   EFtmp_.set_float_part(float_part().sum());
  //   EFtmp_.set_exponent_part(exponent_part());
  //   return EFtmp_;
  // }
};
}

#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_EXTENDEDFLOATEIGENTOOLS_H
