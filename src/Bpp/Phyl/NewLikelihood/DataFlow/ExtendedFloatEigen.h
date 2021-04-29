//
// File: Extendedfloatmatrix.h
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

#ifndef EXTENDED_FLOAT_EIGEN_H
#define EXTENDED_FLOAT_EIGEN_H

#include "ExtendedFloat.h"
#include "ExtendedFloatVectorwiseOp.h"

namespace bpp {

  /*
   * Base Class to allow generic type declaration with no
   * consideration on dimensions. This class knows the
   * Extendedfloatmatrix from which it is inherited (function derived)
   *
   *
   */
  
  template<typename Derived>
  class ExtendedFloatEigenBase {
  private:
    Derived& der_;

  public:
    using Scalar = double;
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
  };

  /*
   * Class associating Eigen::Matrix & exponant to hable underflow.
   *
   * BEWARE: Most operators (such as operator (), = ) handle only Eigen object.
   *  
   *
   *
   */

  template <int R, int C>
  using EFMatrix = Eigen::Matrix<double, R, C>;
  
  template <int R, int C>
  using EFArray = Eigen::Array<double, R, C>;

  template<int R, int C, template<int R2, int C2> class EigenType>
  class ExtendedFloatEigen : public ExtendedFloatEigenBase<ExtendedFloatEigen<R, C, EigenType>>  {
    using ExtType = int;
    using MatType = EigenType<R,C>;

    using Self = ExtendedFloatEigen<R, C, EigenType>;

    using Base = ExtendedFloatEigenBase<Self>;//<double, R, C>>;

    template <int R2, int C2>
    using ExtendedFloatMatrix =  ExtendedFloatEigen<R2, C2, EFMatrix>;

    template <int R2, int C2>
    using ExtendedFloatArray =  ExtendedFloatEigen<R2, C2, EFArray>;


  private:
    MatType mat_;
    ExtType exp_;

  public:
    ExtendedFloatEigen(void):
      ExtendedFloatEigenBase<Self>(*this),
      mat_(MatType()),
      exp_(0)
    {}

    ExtendedFloatEigen(Eigen::DenseBase<MatType>& mat,
                       int exp = 0) :
      ExtendedFloatEigenBase<Self>(*this),
      mat_(mat.derived()),
      exp_(exp) 
    {  
    }

    ExtendedFloatEigen(const Eigen::DenseBase<MatType>& mat,
                       int exp = 0) :
      ExtendedFloatEigenBase<Self>(*this),
      mat_(mat.derived()),
      exp_(exp) 
    {  
    }

    ExtendedFloatEigen(const ExtendedFloatEigenBase<Self>& other) :
      ExtendedFloatEigenBase<Self>(*this),
      mat_(other.derived().mat_),
      exp_(other.derived().exp_) {}

    ExtendedFloatEigen(const ExtendedFloatEigen& other) :
      ExtendedFloatEigenBase<Self>(*this),
      mat_(other.mat_),
      exp_(other.exp_) {}

    ExtendedFloatEigen(const MatType& mat,
                       int exp = 0) :
      ExtendedFloatEigenBase<Self>(*this),
      mat_(mat),
      exp_(exp) 
    {  
    }

    ExtendedFloatEigen(int rows, int cols) :
      ExtendedFloatEigenBase<Self>(*this),
      mat_(MatType::Zero(rows, cols)),
      exp_(0) {}
    
    ExtendedFloatEigen(int cols):
      ExtendedFloatEigenBase<Self>(*this),
      mat_(MatType::Zero(cols)),
      exp_(0) {
    }

    ExtendedFloatEigen& operator=(const ExtendedFloatEigen& other)
    {
      mat_ = other.mat_;
      exp_ = other.exp_;
      return *this;
    }

    // access members
    
    const ExtType & exponent_part () const { return exp_; }
    
    const MatType & float_part () const { return mat_; }

    ExtType & exponent_part () noexcept { return exp_; }

    MatType & float_part () noexcept { return mat_; }


    // Normalization methods
    
    void normalize_big () noexcept {
      // if (std::isfinite (f_)) {
      //   while (std::abs(f_) > biggest_normalized_value) {
      //     f_ *= normalize_big_factor;
      //     exp_ += biggest_normalized_radix_power;
      //   }
      // }
    }
    
    void normalize_small () {
      // if (f_!=0) {
      //   while (std::abs(f_) < smallest_normalized_value) {
      //     f_ *= normalize_small_factor;
      //     exp_ += smallest_normalized_radix_power;
      //   }
      // }
    }
    
    void normalize () noexcept {
      normalize_big ();
      normalize_small ();
    }


    // Static methods without normalization
    inline static Self denorm_mul (const Self & lhs, const Self & rhs) {
      return {lhs.float_part () * rhs.float_part (), lhs.exponent_part () + rhs.exponent_part ()};
    }

    inline static Self denorm_dot (const Self & lhs, const Self & rhs) {
      return {lhs.float_part ().dot(rhs.float_part ()), lhs.exponent_part () + rhs.exponent_part ()};
    }

    inline static Self denorm_div (const Self & lhs, const Self & rhs) {
      return {lhs.float_part () / rhs.float_part (), lhs.exponent_part () - rhs.exponent_part ()};
    }

    inline static Self denorm_add (const Self & lhs, const Self & rhs) {
      return (lhs.exponent_part ()>=rhs.exponent_part ())?
        Self(lhs.float_part () + rhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, rhs.exponent_part () - lhs.exponent_part ()), lhs.exponent_part ()):
        Self(rhs.float_part () + lhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, lhs.exponent_part () - rhs.exponent_part ()), rhs.exponent_part ());
    }

    inline static Self denorm_sub (const Self & lhs, const Self & rhs) {
      return (lhs.exponent_part ()>=rhs.exponent_part ())?
        Self(lhs.float_part () - rhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, rhs.exponent_part () - lhs.exponent_part ()), lhs.exponent_part ()):
        Self(rhs.float_part () - lhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, lhs.exponent_part () - rhs.exponent_part ()), rhs.exponent_part ());
    }

    inline static Self denorm_pow (const Self& arr, double exp)
    {
      double b=arr.exponent_part()*exp;
      ExtendedFloat::ExtType e=ExtendedFloat::ExtType(std::lround(b));
      double rs=std::pow(ExtendedFloat::radix,(b-e));
      return Self(arr.float_part().unaryExpr([exp,rs](double x){return(std::pow(x,exp) * rs);}),e);
    }

    inline static Self denorm_pow (const Self& arr, int exp)
    {
      if (exp==0)
        return Self::Ones(arr.rows(),arr.cols);
      if (exp & 1)
        return (exp>0? denorm_mul(arr,denorm_pow(arr, exp-1)): denorm_div(denorm_pow(arr, exp+1),arr));
      else
      {
        ExtendedFloatArray<R,C> r2(arr.float_part().square());
        r2.normalize();
        auto r2k=denorm_pow(r2,exp>>1);
        r2k.normalize();
        return Self(r2k.float_part(), r2k.exponent_part() + arr.exponent_part()*exp);
      }
    }

    /*********************************
     ** Utilities
     *********************************/

    inline Self operator+ (const Self & rhs) const {
      auto r = denorm_add (*this, rhs);
      r.normalize ();
      return r;
    }

    inline Self operator- (const Self & rhs) const {
      auto r = denorm_sub (*this, rhs);
      r.normalize ();
      return r;
    }

    inline Self operator* (const Self & rhs) const {
      auto r = denorm_mul (*this, rhs);
      r.normalize ();
      return r;
    }

    inline Self dot (const Self & rhs) const {
      auto r = denorm_dot (*this, rhs);
      r.normalize ();
      return r;
    }

    inline Self operator/ (const Self & rhs) const {
      auto r = denorm_div (*this, rhs);
      r.normalize ();
      return r;
    }

    inline Self& operator*= (const Self & rhs) {
      float_part () *= rhs.float_part ();
      exponent_part () += rhs.exponent_part ();
      normalize();
      return *this;
    }
    
    inline Self& operator/= (const Self & rhs) {
      float_part () /= rhs.float_part ();
      exponent_part () -= rhs.exponent_part ();
      normalize();
      return *this;
    }

    inline Self& operator+= (const Self & rhs) {
      if (exponent_part ()>=rhs.exponent_part ())
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
    
    inline Self& operator-= (const Self & rhs) {
      if (exponent_part ()>=rhs.exponent_part ())
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
    
    inline Self operator- () const {
      return Self(-float_part(), exponent_part());
    }

    /*
     * Tests
     *
     */
    
    inline bool operator== (const Self & rhs) const {
      return (float_part()==rhs.float_part() && exponent_part()==rhs.exponent_part());
    }

    inline bool operator!= (const Self & rhs) const {
      return (float_part()!=rhs.float_part() || exponent_part()!=rhs.exponent_part());
    }

    template<int R2, int C2, template<int R3=R2, int C3=C2> class ET2> 
    ExtendedFloatEigen& operator=(const ExtendedFloatEigen<R2, C2, ET2>& other)
    {
      mat_ = other.mat_;
      exp_ = other.exp_;
      return *this;
    }

    /*
     * Eigen like operators
     *
     */
     
    void fill(double val)
    {
      mat_.fill(val);
    }

    Eigen::Index cols() const
    {
      return mat_.cols();
    }

    Eigen::Index rows() const
    {
      return mat_.rows();
    }

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

    static Self Identity(Eigen::Index rows, Eigen::Index cols)
    {
      return Self(MatType::Identity(rows, cols), 0);
    }

    static Self Constant(Eigen::Index rows, Eigen::Index cols, double value)
    {
      return Self(MatType::Constant(rows, cols, value), 0);
    }

    static Self Constant(Eigen::Index rows, double value)
    {
      return Self(MatType::Constant(rows, value), 0);
    }

    template<typename CustomNullaryOp>
    static Self NullaryExpr(Eigen::Index rows, Eigen::Index cols, const CustomNullaryOp& func)
    {
      return Self(MatType::NullaryExpr(rows, cols, func), 0);
    }

    // template<typename CustomNullaryOp>
    Self unaryExpr(const double& func) const
    {
      return Self(float_part().rows(), float_part().cols());//, func), 0);
    }

    const double& operator()(Eigen::Index row, Eigen::Index col) const
    {
      return mat_(row,col);
    }

    ExtendedFloatVectorwiseOp<Self, MatType, Eigen::Vertical> colwise()
    {
      return ExtendedFloatVectorwiseOp<Self, MatType, Eigen::Vertical>(*this);
    }

    ExtendedFloatVectorwiseOp<Self, MatType, Eigen::Horizontal> rowwise()
    {
      return ExtendedFloatVectorwiseOp<Self, MatType, Eigen::Horizontal>(*this);
    }

    double& operator()(Eigen::Index row, Eigen::Index col)
    {
      return float_part()(row,col);
    }

    double& operator()(Eigen::Index row)
    {
      return float_part()(row);
    }

    const double& operator()(Eigen::Index row) const
    {
      return float_part()(row);
    }

    ExtendedFloatMatrix<C, R> transpose() const
    {
      return ExtendedFloatMatrix<C, R>(mat_.transpose(), exp_);
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

    inline double maxCoeff() const
    {
      return convert(ExtendedFloat(float_part().Maxcoeff(),exponent_part()));
    }

  };

  template<int R, int C, template<int R2, int C2> class EigenType>
  ExtendedFloatEigen<R,C,EigenType> operator*(const ExtendedFloatEigen<R, C, EigenType>& mat, const double fact)
  {
    return ExtendedFloatEigenBase<ExtendedFloatEigen<R, C, EigenType>>  (mat.float_part()*fact,mat.exponent_part());
  }

  template<int R, int C, template<int R2, int C2> class EigenType>
  ExtendedFloatEigen<R,C,EigenType> operator*(const double fact, const ExtendedFloatEigen<R, C, EigenType>& mat)
  {
    return mat * fact;
  }

  /*
   * Convenient shortnames.
   *
   */

  template <int R, int C>
  using ExtendedFloatMatrix =  ExtendedFloatEigen<R, C, EFMatrix>;

  template <int R, int C>
  using ExtendedFloatArray =  ExtendedFloatEigen<R, C, EFArray>;

  template <int R>
  using ExtendedFloatVector =  ExtendedFloatMatrix<R, 1>;

  template <int C>
  using ExtendedFloatRowVector = ExtendedFloatMatrix<1, C>;

  typedef ExtendedFloatMatrix<Eigen::Dynamic, Eigen::Dynamic> ExtendedFloatMatrixXd;

  typedef ExtendedFloatRowVector<Eigen::Dynamic> ExtendedFloatRowVectorXd;

  typedef ExtendedFloatVector<Eigen::Dynamic> ExtendedFloatVectorXd;

  /* Specific Methods */

  template<int R, int C>
  static inline ExtendedFloatArray<R, C> log (const ExtendedFloatArray<R, C>& arr)
  {
    return ExtendedFloatArray<R,C>(arr.float_part ().log() + static_cast<double> (arr.exponent_part ()) * ExtendedFloat::ln_radix, 0);
  }

  
  template<int R, int C>
  static inline ExtendedFloatArray<R, C> exp (const ExtendedFloatArray<R, C>& arr)
  {
    // look for max, ie most important exp
    Eigen::Index maxRow, maxCol;
    const auto& arrf = arr.float_part();
    auto max = arrf.maxCoeff(&maxRow, &maxCol);

    auto rcoeff = ExtendedFloat(max/ExtendedFloat::ln_radix, arr.exponent_part()).lround();

    auto c = arrf/ExtendedFloat::ln_radix;
    c.unaryExpr([rcoeff](double x){return(std::pow(ExtendedFloat::radix, x-(double)std::get<0>(rcoeff)));});
    
    // more precision for the max
    c(maxRow, maxCol)=std::pow(ExtendedFloat::radix, std::get<1>(rcoeff));
    
    ExtendedFloatArray<R,C>  expM(c, std::get<0>(rcoeff));
    expM.normalize();
    return(expM);
  }

  template<int R, int C>
  static inline ExtendedFloatArray<R, C> pow (const ExtendedFloatArray<R, C>& obj, double exp)
  {
    auto r = denorm_pow(obj, exp);
    r.normalize ();
    return r;
  }

  template<int R, int C>
  static inline ExtendedFloatArray<R, C> pow (const ExtendedFloatArray<R, C>& obj, int exp)
  {
    auto  r = denorm_pow(obj, exp);
    r.normalize ();
    return r;
  }


  
}

#endif // EXTENDEDFLOATEIGEN_H
