//
// File: ExtendedFloatEigen.h
// Authors:
//   Laurent GuÃÂ©guen (2021)
// Created: samedi 17 avril 2021, ÃÂ  08h 11
// Last modified: 2018-07-11 00:00:00
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_EXTENDEDFLOATEIGEN_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_EXTENDEDFLOATEIGEN_H


#include "ExtendedFloat.h"
#include "ExtendedFloatEigenTools.h"

namespace bpp
{
/*
 * Base Class to allow generic type declaration with no
 * consideration on dimensions. This class knows the
 * Extendedfloatmatrix from which it is inherited (function derived)
 *
 *
 */

struct ExtendedFloatEigenCore {};

template<typename Derived>
class ExtendedFloatEigenBase :
  public ExtendedFloatEigenCore
{
private:
  Derived& der_;

public:
  using Scalar = ExtendedFloat;
  using RealScalar = double;

  ExtendedFloatEigenBase(Derived& der) :
    der_(der) {}

  Derived& derived()
  {
    return der_;
  }

  const Derived& derived() const
  {
    return der_;
  }

  // No Alias

  ExtendedFloatNoAlias<Derived> noalias()
  {
    return ExtendedFloatNoAlias<Derived>(derived());
  }
};


template< int R,  int C>
using EFMatrix = Eigen::Matrix<double, R, C>;

template< int R,  int C>
using EFArray = Eigen::Array<double, R, C>;

/*
 * Class associating Eigen::Matrix & exponant to hanble underflow.
 *
 *
 */

template< int R,  int C, template< int R2,  int C2> class EigenType>
class ExtendedFloatEigen : public ExtendedFloatEigenBase<ExtendedFloatEigen<R, C, EigenType> >
{
  using ExtType =  int;
  using MatType = EigenType<R, C>;

  using RefMatType = Eigen::Ref<EigenType<R, C> >;

  using Self = ExtendedFloatEigen<R, C, EigenType>;

  using Base = ExtendedFloatEigenBase<Self>;

  template< int R2,  int C2>
  using ExtendedFloatMatrix =  ExtendedFloatEigen<R2, C2, EFMatrix>;

  template< int R2,  int C2>
  using ExtendedFloatArray =  ExtendedFloatEigen<R2, C2, EFArray>;

protected:
  MatType mat_;
  ExtType exp_;

  // For specific output without too many creations
  class OwnedExtendedFloat : public ExtendedFloat
  {
private:
    const ExtendedFloatEigen& eigen_;

    void set_float_part(double x)
    {
      f_ = x;
    }
    void set_exponent_part(ExtendedFloat::ExtType x)
    {
      exp_ = x;
    }

public:
    OwnedExtendedFloat(const ExtendedFloatEigen& eigen) :
      ExtendedFloat(),
      eigen_(eigen) {}

    const ExtendedFloat::ExtType& exponent_part() const
    {
      return eigen_.exponent_part();
    }

    friend class ExtendedFloatEigen;
  };

  mutable OwnedExtendedFloat EFtmp_;

public:
  ExtendedFloatEigen(void) :
    ExtendedFloatEigenBase<Self>(*this),
    mat_(MatType()),
    exp_(0),
    EFtmp_(*this)
  {}

  ExtendedFloatEigen(Eigen::DenseBase<MatType>& mat,
                     ExtType exp = 0) :
    ExtendedFloatEigenBase<Self>(*this),
    mat_(mat.derived()),
    exp_(exp),
    EFtmp_(*this)
  {}

  ExtendedFloatEigen(MatType& mat,
                     ExtType exp = 0) :
    ExtendedFloatEigenBase<Self>(*this),
    mat_(mat),
    exp_(exp),
    EFtmp_(*this)
  {}

  ExtendedFloatEigen(const Eigen::DenseBase<MatType>& mat,
                     ExtType exp = 0) :
    ExtendedFloatEigenBase<Self>(*this),
    mat_(mat.derived()),
    exp_(exp),
    EFtmp_(*this)
  {}

  ExtendedFloatEigen(const Eigen::internal::traits<MatType>& mat,
                     ExtType exp = 0) :
    ExtendedFloatEigenBase<Self>(*this),
    mat_(mat),
    exp_(exp),
    EFtmp_(*this)
  {}

  template<class Other>
  ExtendedFloatEigen(const ExtendedFloatEigenBase<Other>& other) :
    ExtendedFloatEigenBase<Self>(*this),
    mat_(other.derived().float_part()),
    exp_(other.derived().exponent_part()),
    EFtmp_(*this) {}

  ExtendedFloatEigen(const ExtendedFloatEigen& other) :
    ExtendedFloatEigenBase<Self>(*this),
    mat_(other.float_part()),
    exp_(other.exponent_part()),
    EFtmp_(*this) {}

  ExtendedFloatEigen(const MatType& mat,
                     ExtType exp = 0) :
    ExtendedFloatEigenBase<Self>(*this),
    mat_(mat),
    exp_(exp),
    EFtmp_(*this)
  {}

  ExtendedFloatEigen( int rows,  int cols) :
    ExtendedFloatEigenBase<Self>(*this),
    mat_(MatType::Zero(rows, cols)),
    exp_(0),
    EFtmp_(*this) {}

  ExtendedFloatEigen( int cols) :
    ExtendedFloatEigenBase<Self>(*this),
    mat_(MatType::Zero(cols)),
    exp_(0),
    EFtmp_(*this) {}

  // Specific constructors

  static Self Zero(Eigen::Index rows, Eigen::Index cols)
  {
    return Self(MatType::Zero(rows, cols), 0);
  }

  static Self Zero(Eigen::Index rows)
  {
    return Self(MatType::Zero(rows), 0);
  }

  static Self Ones(Eigen::Index rows, Eigen::Index cols)
  {
    return Self(MatType::Ones(rows, cols), 0);
  }

  static Self Ones(Eigen::Index rows)
  {
    return Self(MatType::Ones(rows), 0);
  }

  static Self Identity(Eigen::Index rows)
  {
    return Self(MatType::Identity(rows, rows), 0);
  }

  static Self Constant(Eigen::Index rows, Eigen::Index cols, double value)
  {
    return Self(MatType::Constant(rows, cols, value), 0);
  }

  static Self Constant(Eigen::Index rows, Eigen::Index cols, const ExtendedFloat& value)
  {
    return Self(MatType::Constant(rows, cols, value.float_part()), value.exponent_part());
  }

  static Self Constant(Eigen::Index rows, double value)
  {
    return Self(MatType::Constant(rows, value), 0);
  }

  static Self Constant(Eigen::Index rows, const ExtendedFloat& value)
  {
    return Self(MatType::Constant(rows, value.float_part()), value.exponent_part());
  }

  // operators

  ExtendedFloatEigen& operator=(const ExtendedFloatEigen& other)
  {
    mat_ = other.mat_;
    exp_ = other.exp_;
    return *this;
  }

  ExtendedFloatEigen& operator=(const MatType& other)
  {
    mat_ = other;
    exp_ = 0;
    return *this;
  }

  template<typename Derived>
  Self& operator=(const ExtendedFloatEigenBase<Derived>& other)
  {
    mat_ = other.derived().float_part();
    exp_ = other.derived().exponent_part();
    return *this;
  }

  // access members

  const ExtType& exponent_part () const { return exp_; }

  const MatType& float_part () const { return mat_; }

  ExtType& exponent_part () noexcept { return exp_; }

  MatType& float_part () noexcept { return mat_; }

  // Normalization methods

  bool normalize_big () noexcept
  {
    using namespace std;

    auto max = float_part().cwiseAbs().maxCoeff();
    if (isfinite(max))
    {
      bool normalized = false;
      while (max > ExtendedFloat::biggest_normalized_value)
      {
        float_part() *= (double)ExtendedFloat::normalize_big_factor;
        max *= (double)ExtendedFloat::normalize_big_factor;
        exp_ += ExtendedFloat::biggest_normalized_radix_power;
        normalized = true;
      }
      return normalized;
    }
    return false;
  }

  bool normalize_small ()
  {
    const auto& fabs= float_part().cwiseAbs();
    auto max=fabs.maxCoeff();
    if (max > 0){
      // not a vector of zeros
      auto min=fabs.unaryExpr([max](double d){return d>0?d:max;}).minCoeff();
      bool normalized = false;
      while (min< ExtendedFloat::smallest_normalized_value) {
        //to prevent overflow
        if (max * ExtendedFloat::normalize_small_factor >=ExtendedFloat::biggest_normalized_value){
          break;
        } 
        float_part() *= (double)ExtendedFloat::normalize_small_factor;
        min *= (double)ExtendedFloat::normalize_small_factor;
        max *= (double)ExtendedFloat::normalize_small_factor;
        exp_ -= ExtendedFloat::biggest_normalized_radix_power;
        normalized = true;
      }
      return normalized;
    }
    return false;
  }
  
  void normalize () noexcept
  {
    if (!normalize_big())
    {
      normalize_small();
    }
  }

  // Static methods without normalization
  template< int R2,  int C2>
  inline static ExtendedFloatEigen<R, C2, EigenType> denorm_mul (const Self& lhs, const ExtendedFloatEigen<R2, C2, EigenType>& rhs)
  {
    return {lhs.float_part () * rhs.float_part (), lhs.exponent_part () + rhs.exponent_part ()};
  }

  template< typename Derived>
  inline static Self denorm_mul(const Self& lhs, const Eigen::EigenBase<Derived>& rhs)
  {
    return {lhs.float_part () * rhs.derived(), lhs.exponent_part ()};
  }

  template< typename Derived>
  inline static Self denorm_mul(const Eigen::EigenBase<Derived>& lhs, const Self& rhs)
  {
    return {lhs.derived() * rhs.float_part (), rhs.exponent_part ()};
  }

  template<typename T>
  typename std::enable_if< std::is_floating_point<T>::value || std::is_integral<T>::value, Self>::type
  inline static denorm_mul (const Self& lhs, const T& rhs)
  {
    return {lhs.float_part () * rhs, lhs.exponent_part ()};
  }

  inline static Self denorm_mul (const Self& lhs, const double& rhs)
  {
    return {lhs.float_part () * rhs, lhs.exponent_part ()};
  }

  inline static Self denorm_mul (const Self& lhs, const ExtendedFloat rhs)
  {
    return {lhs.float_part () * rhs.float_part (), lhs.exponent_part () + rhs.exponent_part ()};
  }

  template<typename T>
  typename std::enable_if<std::is_same<T, Self>::value || std::is_same<T, ExtendedFloat>::value, Self>::type
  inline static denorm_div (const Self& lhs, const T& rhs)
  {
    return {lhs.float_part () / rhs.float_part (), lhs.exponent_part () - rhs.exponent_part ()};
  }

  template<typename Derived>
  inline Self static denorm_div (const Self& lhs, const ExtendedFloatEigenBase<Derived>& rhs)
  {
    return {lhs.float_part () / rhs.derived().float_part (), lhs.exponent_part () - rhs.derived().exponent_part ()};
  }

  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value || std::is_integral<T>::value, Self>::type
  inline static denorm_add (const Self& lhs, const T& rhs)
  {
    return (lhs.exponent_part () >= 0) ?
           Self(lhs.float_part () + rhs * constexpr_power<double>(ExtendedFloat::radix, -lhs.exponent_part ()), lhs.exponent_part ()) :
           Self(rhs + lhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, lhs.exponent_part ()), 0);
  }

  template<typename T>
  typename std::enable_if<std::is_same<T, Self>::value || std::is_same<T, ExtendedFloat>::value, Self>::type
  inline static denorm_add (const Self& lhs, const T& rhs)
  {
    return (lhs.exponent_part () >= rhs.exponent_part ()) ?
           Self(lhs.float_part () + rhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, rhs.exponent_part () - lhs.exponent_part ()), lhs.exponent_part ()) :
           Self(rhs.float_part () + lhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, lhs.exponent_part () - rhs.exponent_part ()), rhs.exponent_part ());
  }

  template<typename T>
  typename std::enable_if<std::is_same<T, Self>::value || std::is_same<T, ExtendedFloat>::value, Self>::type
  inline static denorm_sub (const Self& lhs, const T& rhs)
  {
    return (lhs.exponent_part () >= rhs.exponent_part ()) ?
           Self(lhs.float_part () - rhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, rhs.exponent_part () - lhs.exponent_part ()), lhs.exponent_part ()) :
           Self(rhs.float_part () - lhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, lhs.exponent_part () - rhs.exponent_part ()), rhs.exponent_part ());
  }


  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value || std::is_integral<T>::value, Self>::type
  inline static denorm_sub (const Self& lhs, const T& rhs)
  {
    return (lhs.exponent_part () >= 0) ?
           Self(lhs.float_part () - rhs * constexpr_power<double>(ExtendedFloat::radix, -lhs.exponent_part ()), lhs.exponent_part ()) :
           Self(rhs - lhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, lhs.exponent_part ()), 0);
  }

  inline static Self denorm_pow (const Self& arr, double exp)
  {
    double b = arr.exponent_part() * exp;
    ExtendedFloat::ExtType e = ExtendedFloat::ExtType(std::lround(b));
    double rs = std::pow(ExtendedFloat::radix, (b - e));
    return Self(arr.float_part().unaryExpr([exp, rs](double x){return std::pow(x, exp) * rs;}), e);
  }

  inline static Self denorm_pow (const Self& arr, int exp)
  {
    if (exp == 0)
      return Self::Ones(arr.rows(), arr.cols);
    if (exp & 1)
      return exp > 0 ? denorm_mul(arr, denorm_pow(arr, exp - 1)) : denorm_div(denorm_pow(arr, exp + 1), arr);
    else
    {
      ExtendedFloatArray<R, C> r2(arr.float_part().square());
      r2.normalize();
      auto r2k = denorm_pow(r2, exp >> 1);
      r2k.normalize();
      return Self(r2k.float_part(), r2k.exponent_part() + arr.exponent_part() * exp);
    }
  }

  /*********************************
  ** Utilities
  *********************************/

  /***
   * Operators
   *
   */

  template<typename T>
  typename std::enable_if<std::is_same<T, Self>::value || std::is_same<T, ExtendedFloat>::value
                          || std::is_floating_point<T>::value || std::is_integral<T>::value, Self>::type
  inline operator+(const T& rhs) const
  {
    auto r = denorm_add (*this, rhs);
    r.normalize ();
    return r;
  }

  template<template<int R2 = R, int C2 = C> class EigenType2>
  inline Self operator+(const ExtendedFloatEigen<R, C, EigenType2>& rhs) const
  {
    auto r = denorm_add (*this, rhs);
    r.normalize ();
    return r;
  }

  template<typename T>
  typename std::enable_if<std::is_same<T, Self>::value || std::is_same<T, ExtendedFloat>::value
                          || std::is_floating_point<T>::value || std::is_integral<T>::value, Self>::type
  inline operator-(const T& rhs) const
  {
    auto r = denorm_sub (*this, rhs);
    r.normalize ();
    return r;
  }

  template<int R2 = C,  int C2, template<int R3 = R2, int C3 = C2> class EigenType2>
  inline ExtendedFloatEigen<R, C2, EigenType> operator*(const ExtendedFloatEigen<R2, C2, EigenType2>& rhs) const
  {
    auto r = denorm_mul (*this, rhs);
    r.normalize ();
    return r;
  }

  // Not sure the product will be fine, in case of dimension mismatch
  template<typename Derived>
  inline Self operator*(const Eigen::EigenBase<Derived>& rhs) const
  {
    auto r = denorm_mul(*this, rhs);
    r.normalize ();
    return r;
  }

  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value || std::is_integral<T>::value || std::is_same<T, ExtendedFloat>::value, Self>::type
  inline operator*(const T& fact) const
  {
    auto r = denorm_mul(*this, fact);
    r.normalize ();
    return r;
  }

  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value || std::is_integral<T>::value || std::is_same<T, ExtendedFloat>::value, Self>::type
  inline operator/(const T& fact) const
  {
    auto r = denorm_div(*this, fact);
    r.normalize ();
    return r;
  }

  template<typename Derived>
  inline Self operator/(const ExtendedFloatEigenBase<Derived>& rhs) const
  {
    auto r = denorm_div(*this, rhs);
    r.normalize ();
    return r;
  }

  template<typename Derived>
  inline ExtendedFloat dot (const ExtendedFloatEigenBase<Derived>& rhs) const
  {
    auto r = ExtendedFloat(float_part().dot(rhs.derived().float_part ()), exponent_part () + rhs.derived().exponent_part ());
    r.normalize();
    return r;
  }

  template<typename Derived>
  inline ExtendedFloat dot (const Eigen::DenseBase<Derived>& rhs) const
  {
    auto r = ExtendedFloat(float_part().dot(rhs.derived ()), exponent_part ());
    r.normalize();
    return r;
  }

  template<typename Obj>
  inline ExtendedFloat dot (const Eigen::Ref<Obj>& rhs) const
  {
    auto r = ExtendedFloat(float_part().dot(rhs.derived ()), exponent_part ());
    r.normalize();
    return r;
  }

  inline Self operator-() const
  {
    return Self(-float_part(), exponent_part());
  }

  inline Self eval() const
  {
    return Self(float_part().eval(), exponent_part());
  }


  /*
   * Modifying operators
   *
   */

  template<typename T>
  typename std::enable_if<std::is_same<T, Self>::value || std::is_same<T, ExtendedFloat>::value, Self&>::type
  inline operator*=(const T& rhs)
  {
    float_part () *= rhs.float_part ();
    exponent_part () += rhs.exponent_part ();
    normalize();
    return *this;
  }

  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value || std::is_integral<T>::value, Self&>::type
  inline operator*=(const T& div)
  {
    float_part () *= div;
    normalize();
    return *this;
  }

  template<typename Derived>
  inline Self& operator*=(const Eigen::DenseBase<Derived>& div)
  {
    float_part () *= div.derived().float_part();
    normalize();
    return *this;
  }

  template<typename Derived>
  inline Self& operator*=(const ExtendedFloatEigenBase<Derived>& div)
  {
    float_part () *= div.derived().float_part();
    exponent_part () += div.derived().exponent_part ();
    normalize();
    return *this;
  }


  template<typename T>
  typename std::enable_if<std::is_same<T, Self>::value || std::is_same<T, ExtendedFloat>::value, Self&>::type
  inline operator/=(const T& rhs)
  {
    float_part () /= rhs.float_part ();
    exponent_part () -= rhs.exponent_part ();
    normalize();
    return *this;
  }

  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value || std::is_integral<T>::value, Self&>::type
  inline operator/=(const T& div)
  {
    float_part () /= div;
    return *this;
  }


  template<typename T>
  typename std::enable_if<std::is_base_of<ExtendedFloatEigenCore, T>::value || std::is_same<T, ExtendedFloat>::value, Self&>::type
  inline operator+=(const T& rhs)
  {
    if (exponent_part () >= rhs.exponent_part ())
    {
      float_part() += rhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, rhs.exponent_part () - exponent_part ());
    }
    else
    {
      float_part() = rhs.float_part () + float_part () * constexpr_power<double>(ExtendedFloat::radix, exponent_part () - rhs.exponent_part ());
      exponent_part() = rhs.exponent_part ();
    }
    normalize ();
    return *this;
  }

  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value || std::is_integral<T>::value, Self&>::type
  inline operator+=(const T& rhs)
  {
    if (exponent_part () >= 0)
    {
      float_part() += rhs * constexpr_power<double>(ExtendedFloat::radix, -exponent_part ());
    }
    else
    {
      float_part() = rhs.float_part () + float_part () * constexpr_power<double>(ExtendedFloat::radix, exponent_part ());
      exponent_part() = 0;
    }
    normalize ();
    return *this;
  }

  template<typename T>
  typename std::enable_if<std::is_same<T, Self>::value || std::is_same<T, ExtendedFloat>::value, Self&>::type
  inline operator-=(const T& rhs)
  {
    if (exponent_part () >= rhs.exponent_part ())
    {
      float_part() -= rhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, rhs.exponent_part () - exponent_part ());
    }
    else
    {
      float_part() = rhs.float_part () - float_part () * constexpr_power<double>(ExtendedFloat::radix, exponent_part () - rhs.exponent_part ());
      exponent_part() = rhs.exponent_part ();
    }
    normalize ();
    return *this;
  }

  template<typename T>
  typename std::enable_if<std::is_floating_point<T>::value || std::is_integral<T>::value, Self&>::type
  inline operator-=(const T& rhs)
  {
    if (exponent_part () >= 0)
    {
      float_part() -= rhs * constexpr_power<double>(ExtendedFloat::radix, -exponent_part ());
    }
    else
    {
      float_part() = float_part () * constexpr_power<double>(ExtendedFloat::radix, exponent_part ()) - rhs.float_part ();
      exponent_part() = 0;
    }
    normalize ();
    return *this;
  }


  inline Self log () const
  {
    return Self(float_part ().log() + static_cast<double>(exponent_part ()) * ExtendedFloat::ln_radix, 0);
  }

  inline Self exp () const
  {
    // look for max, ie most important exp
    Eigen::Index maxRow, maxCol;
    const auto& arrf = float_part();
    auto max = arrf.maxCoeff(&maxRow, &maxCol);

    auto rcoeff = ExtendedFloat(max / ExtendedFloat::ln_radix, exponent_part()).lround();

    MatType c(arrf / ExtendedFloat::ln_radix);
    c.unaryExpr([rcoeff](double x){return std::pow(ExtendedFloat::radix, x - (double)std::get<0>(rcoeff));});

    // more precision for the max
    c(maxRow, maxCol) = std::pow(ExtendedFloat::radix, std::get<1>(rcoeff));

    Self expM(c, std::get<0>(rcoeff));
    expM.normalize();
    return expM;
  }

  /*
   * Tests
   *
   */
  inline bool operator==(const Self& rhs) const
  {
    return float_part() == rhs.float_part() && exponent_part() == rhs.exponent_part();
  }

  inline bool operator!=(const Self& rhs) const
  {
    return float_part() != rhs.float_part() || exponent_part() != rhs.exponent_part();
  }

  /*
   * Eigen like operators
   *
   */
  void fill(double val)
  {
    float_part().fill(val);
    exponent_part() = 0;
  }

  Eigen::Index cols() const
  {
    return float_part().cols();
  }

  Eigen::Index rows() const
  {
    return float_part().rows();
  }

  template<typename CustomNullaryOp>
  static Self NullaryExpr(Eigen::Index rows, Eigen::Index cols, const CustomNullaryOp& func)
  {
    return Self(MatType::NullaryExpr(rows, cols, func), func.exponent_part());
  }

  ExtendedFloatVectorwiseOp<Self, MatType, Eigen::Vertical> colwise()
  {
    return ExtendedFloatVectorwiseOp<Self, MatType, Eigen::Vertical>(*this);
  }

  ExtendedFloatVectorwiseOp<const Self, const MatType, Eigen::Vertical> colwise() const
  {
    return ExtendedFloatVectorwiseOp<const Self, const MatType, Eigen::Vertical>(*this);
  }

  template<typename M = MatType>
  typename std::enable_if<std::is_same<M, EFMatrix<R, C> >::value, ExtendedFloatRow<R, C, EigenType > >::type
  row(Eigen::Index pos)
  {
    return ExtendedFloatRow<R, C, EigenType>(*this, pos);
  }

  ExtendedFloatVectorwiseOp<Self, MatType, Eigen::Horizontal> rowwise()
  {
    return ExtendedFloatVectorwiseOp<Self, MatType, Eigen::Horizontal>(*this);
  }

  ExtendedFloatVectorwiseOp<const Self, const MatType, Eigen::Horizontal> rowwise() const
  {
    return ExtendedFloatVectorwiseOp<const Self, const MatType, Eigen::Horizontal>(*this);
  }


  /****
   * @brief Access to elements of the matrix
   *
   */
  const ExtendedFloat& operator()(Eigen::Index c) const
  {
    EFtmp_.set_float_part(float_part()(c));
    EFtmp_.set_exponent_part(exponent_part());
    return EFtmp_;
  }

  const ExtendedFloat& operator()(Eigen::Index r, Eigen::Index c) const
  {
    EFtmp_.set_float_part(float_part()(r, c));
    EFtmp_.set_exponent_part(exponent_part());
    return EFtmp_;
  }

  template<typename M = MatType>
  typename std::enable_if<std::is_same<M, EFArray<R, C> >::value, const ExtendedFloat&>::type
  operator[](Eigen::Index pos) const
  {
    EFtmp_.set_exponent_part(exponent_part());
    EFtmp_.set_float_part(float_part()[pos]);
    return EFtmp_;
  }

  const ExtendedFloat& sum() const
  {
    EFtmp_.set_float_part(float_part().sum());
    EFtmp_.set_exponent_part(exponent_part());
    return EFtmp_;
  }

  const ExtendedFloat& mean() const
  {
    EFtmp_.set_exponent_part(exponent_part());
    EFtmp_.set_float_part(float_part().mean());
    return EFtmp_;
  }

  const ExtendedFloat& maxCoeff(size_t* pos = 0) const
  {
    EFtmp_.set_exponent_part(exponent_part());
    if (pos)
      EFtmp_.set_float_part(float_part().maxCoeff(pos));
    else
      EFtmp_.set_float_part(float_part().maxCoeff());
    return EFtmp_;
  }

  ExtendedFloatEigen<R, 1, EigenType> col(Eigen::Index col) const
  {
    return ExtendedFloatEigen<R, 1, EigenType>(float_part().col(col), exponent_part());
  }

  /*********************************************/
  /*** Modifications  ******/

  ExtendedFloatMatrix<C, R> transpose() const
  {
    return ExtendedFloatMatrix<C, R>(float_part().transpose(), exponent_part());
  }

  void resize(Eigen::Index rows, Eigen::Index cols)
  {
    float_part().resize(rows, cols);
  }

  void resize(Eigen::Index rows)
  {
    float_part().resize(rows);
  }

  Eigen::Index size() const
  {
    return float_part().size();
  }

  /**************************************/
  /* Output */

  friend std::ostream& operator<<(std::ostream& out, const Self& ef)
  {
    return out << "( " << ef.float_part() << " ) *  2^" << ef.exponent_part();
  }
};


/*
 * Convenient shortnames.
 *
 */

template< int R,  int C>
using ExtendedFloatMatrix =  ExtendedFloatEigen<R, C, EFMatrix>;

template< int R,  int C>
using ExtendedFloatArray =  ExtendedFloatEigen<R, C, EFArray>;

template< int R>
using ExtendedFloatVector =  ExtendedFloatMatrix<R, 1>;

template< int C>
using ExtendedFloatRowVector = ExtendedFloatMatrix<1, C>;

typedef ExtendedFloatMatrix<Eigen::Dynamic, Eigen::Dynamic> ExtendedFloatMatrixXd;

typedef ExtendedFloatRowVector<Eigen::Dynamic> ExtendedFloatRowVectorXd;

typedef ExtendedFloatVector<Eigen::Dynamic> ExtendedFloatVectorXd;

typedef ExtendedFloatArray<Eigen::Dynamic, 1> ExtendedFloatArrayXd;


/*  Extern Methods */

template< int R,  int C>
inline ExtendedFloatArray<R, C> log (const ExtendedFloatArray<R, C>& arr)
{
  return arr.log();
}


template< int R,  int C, template< int R2 = R,  int C2 = C> class MatType>
inline ExtendedFloatEigen<R, C, MatType> exp (const ExtendedFloatEigen<R, C, MatType>& arr)
{
  return arr.exp();
}

template< int R,  int C>
inline ExtendedFloatArray<R, C> pow (const ExtendedFloatArray<R, C>& obj, double exp)
{
  auto r = ExtendedFloatArray<R, C>::denorm_pow(obj, exp);
  r.normalize ();
  return r;
}

template< int R,  int C>
inline ExtendedFloatArray<R, C> pow (const ExtendedFloatArray<R, C>& obj, int exp)
{
  auto r = ExtendedFloatArray<R, C>::denorm_pow(obj, exp);
  r.normalize ();
  return r;
}


template<int R, int C, template< int R2 = R,  int C2 = C> class MatType, typename T>
typename std::enable_if<std::is_same<T, ExtendedFloat>::value
  || std::is_floating_point<T>::value || std::is_integral<T>::value,
                        ExtendedFloatEigen<R, C, MatType> >::type
inline operator+(const T& d, const ExtendedFloatEigen<R, C, MatType> rhs)
{
  return rhs + d;
}

template<int R, int C, template< int R2 = R, int C2 = C> class MatType, typename T>
typename std::enable_if<std::is_same<T, ExtendedFloat>::value
                        || std::is_floating_point<T>::value || std::is_integral<T>::value,
                        ExtendedFloatEigen<R, C, MatType> >::type
inline operator-(const T& d, const ExtendedFloatEigen<R, C, MatType> rhs)
{
  return -(rhs - d);
}

template<int R, int C, template< int R2 = R,  int C2 = C> class MatType, typename T>
typename std::enable_if<std::is_same<T, ExtendedFloat>::value
                        || std::is_floating_point<T>::value || std::is_integral<T>::value,
                        ExtendedFloatEigen<R, C, MatType> >::type
inline operator*(const T& fact, const ExtendedFloatEigen<R, C, MatType> mat)
{
  auto r = ExtendedFloatEigen<R, C, MatType>::denorm_mul (mat, fact);
  r.normalize ();
  return r;
}

template<typename Derived, typename EFType>
inline EFType operator*(const Eigen::EigenBase<Derived>& lhs,
                        const ExtendedFloatEigenBase<EFType>& rhs)
{
  auto r = EFType::denorm_mul(lhs.derived(), rhs);
  r.normalize ();
  return r;
}


template< int R,  int C>
class ExtendedFloatArrayWrapper
{
  using Self = ExtendedFloatArrayWrapper<R, C>;

  using MatType = EFArray<R, C>;

  using ExtType =  int;

  using Array = ExtendedFloatArray<R, C>;

private:
  ExtendedFloatMatrix<R, C>* const efm_;

public:
  ExtendedFloatArrayWrapper(ExtendedFloatMatrix<R, C>& other) :
    efm_(&other) {}

  Self operator=(const Array& rhs)
  {
    efm_->float_part().array() = rhs.float_part();
    efm_->exponent_part() = rhs.exponent_part();
    return *this;
  }

  Self operator+=(const Array& rhs)
  {
    if (efm_->float_part().isZero())
      return operator=(rhs);

    if (efm_->exponent_part () >= rhs.exponent_part ())
    {
      efm_->float_part().array() += rhs.float_part() * constexpr_power<double>(ExtendedFloat::radix, rhs.exponent_part () - efm_->exponent_part ());
    }
    else
    {
      efm_->float_part().array() = rhs.float_part () + efm_->float_part ().array() * constexpr_power<double>(ExtendedFloat::radix, efm_->exponent_part () - rhs.exponent_part ());
      efm_->exponent_part() = rhs.exponent_part ();
    }
    efm_->normalize ();

    return *this;
  }

  Self operator*=(const Array& rhs)
  {
    efm_->float_part().array() *= rhs.float_part ();
    efm_->exponent_part() += rhs.exponent_part ();
    efm_->normalize ();
    return *this;
  }

  Array operator*(const Array& rhs)  const
  {
    Array r(efm_->float_part().array() * rhs.float_part (), efm_->exponent_part () + rhs.exponent_part ());
    r.normalize();
    return r;
  }
};

template<int R, int C>
inline ExtendedFloatArray<R, C> operator*(const ExtendedFloatArray<R, C>& lhs,
                                          const ExtendedFloatArrayWrapper<R, C>& rhs)
{
  return rhs * lhs; // cwise * is commutative
}
}
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_EXTENDEDFLOATEIGEN_H
