//
// File: ExtendedFloat.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-07-04 00:00:00
// Last modified: 2017-07-06
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

#ifndef BPP_NEWPHYL_EXTENDEDFLOAT_H
#define BPP_NEWPHYL_EXTENDEDFLOAT_H

#include <cmath>
#include <limits>
#include <string>
#include <iostream>

namespace bpp {

  template <typename T> constexpr T constexpr_power (T d, int n) {
    return n == 0 ? 1.0 : (n > 0 ? constexpr_power (d, n - 1) * d : constexpr_power (d, n + 1) / d);
  }

  
  class ExtendedFloat {
    // Assumes positive integer
  public:
    using FloatType = double;
    using ExtType = int;

    // Parameter: decide how much product we can do safely before having to normalize (smaller -> less normalizations)
    static constexpr int allowed_product_without_normalization = 50;

    // Radix is the float exponent base
    static constexpr int radix = std::numeric_limits<FloatType>::radix;

    // biggest_repr_radix_power = max { n ; radix^n is representable }
    static constexpr int biggest_repr_radix_power = std::numeric_limits<FloatType>::max_exponent - 1;
    
    // biggest_normalized_radix_power = max { n ; (radix^n)^allowed_product_without_normalization is representable }
    static constexpr int biggest_normalized_radix_power =
      biggest_repr_radix_power / allowed_product_without_normalization;

    // biggest_normalized_value = max { f ; f^allowed_product_without_normalization is representable }
    static constexpr FloatType biggest_normalized_value =
      constexpr_power (FloatType (radix), biggest_normalized_radix_power);

    // smallest_repr_radix_power = min { n ; radix^n is representable }
    static constexpr int smallest_repr_radix_power = std::numeric_limits<FloatType>::min_exponent - 1;

    // smallest_normalized_radix_power = min { n ; (radix^n)^allowed_product_denorm is representable }
    static constexpr int smallest_normalized_radix_power =
      -((-smallest_repr_radix_power) / allowed_product_without_normalization);

    // smallest_normalized_value = min { f ; f^allowed_product_denorm is representable }
    static constexpr FloatType smallest_normalized_value =
      constexpr_power (FloatType (radix), smallest_normalized_radix_power);

    // factors to scale f_ to renormalize.
    static constexpr FloatType normalize_big_factor = 1. / biggest_normalized_value;
    static constexpr FloatType normalize_small_factor = 1. / smallest_normalized_value;

    // TODO add denorm info for sum

    constexpr ExtendedFloat (FloatType f = 0.0, ExtType e = 0) noexcept : f_ (f), exp_ (e) {}

    const FloatType & float_part () const noexcept { return f_; }
    const ExtType & exponent_part () const noexcept { return exp_; }

    void normalize_big () noexcept {
      if (std::isfinite (f_)) {
        while (f_ > biggest_normalized_value) {
          f_ *= normalize_big_factor;
          exp_ += biggest_normalized_radix_power;
        }
      }
    }
    void normalize_small () {
      if (f_ > 0.) {
        while (f_ < smallest_normalized_value) {
          f_ *= normalize_small_factor;
          exp_ += smallest_normalized_radix_power;
        }
      }
    }
    void normalize () noexcept {
      normalize_big ();
      normalize_small ();
    }

  private:
    FloatType f_;
    ExtType exp_;
  };

  inline std::string to_string (const ExtendedFloat & ef) {
    using std::to_string;
    return "double(" + to_string (ef.float_part ()) + " * 2^" + to_string (ef.exponent_part ()) + ")";
  }

  inline ExtendedFloat denorm_mul (const ExtendedFloat & lhs, const ExtendedFloat & rhs) {
    return {lhs.float_part () * rhs.float_part (), lhs.exponent_part () + rhs.exponent_part ()};
  }

  inline ExtendedFloat denorm_div (const ExtendedFloat & lhs, const ExtendedFloat & rhs) {
    return {lhs.float_part () / rhs.float_part (), lhs.exponent_part () - rhs.exponent_part ()};
  }

  inline ExtendedFloat denorm_add (const ExtendedFloat & lhs, const ExtendedFloat & rhs) {
    return (lhs.exponent_part ()>=rhs.exponent_part ())?
      ExtendedFloat(lhs.float_part () + rhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, rhs.exponent_part () - lhs.exponent_part ()), lhs.exponent_part ()):
      ExtendedFloat(rhs.float_part () + lhs.float_part () * constexpr_power<double>(ExtendedFloat::radix, lhs.exponent_part () - rhs.exponent_part ()), rhs.exponent_part ());
  }

  inline ExtendedFloat denorm_pow (const ExtendedFloat & lhs, double exp) {
    double b=lhs.exponent_part()*exp;
    ExtendedFloat::ExtType e=ExtendedFloat::ExtType(lround(b));
    ExtendedFloat r(std::pow(lhs.float_part(),exp)*std::pow(ExtendedFloat::radix,(b-e)), e);
    return r;
  }

  inline ExtendedFloat denorm_pow (const ExtendedFloat & lhs, int exp) {
    if (exp==0)
      return ExtendedFloat(1.0);
    if (exp & 1)
      return (exp>0? denorm_mul(lhs,denorm_pow(lhs, exp-1)): denorm_div(denorm_pow(lhs, exp+1),lhs));
    else
    {
      ExtendedFloat r2(pow(lhs.float_part(),2));
      auto r2k=denorm_pow(r2,exp>>1);
      r2k.normalize_small();
      ExtendedFloat r(r2k.float_part(), r2k.exponent_part() + lhs.exponent_part()*exp);
      return r;
    }
  }

  inline ExtendedFloat operator* (const ExtendedFloat & lhs, const ExtendedFloat & rhs) {
    auto r = denorm_mul (lhs, rhs);
    r.normalize ();
    return r;
  }

  inline ExtendedFloat operator/ (const ExtendedFloat & lhs, const ExtendedFloat & rhs) {
    auto r = denorm_div (lhs, rhs);
    r.normalize ();
    return r;
  }

  inline ExtendedFloat pow (const ExtendedFloat & lhs, double exp) {
    auto  r = denorm_pow(lhs, exp);
    r.normalize ();
    return r;
  }

  inline ExtendedFloat pow (const ExtendedFloat & lhs, int exp) {
    auto  r = denorm_pow(lhs, exp);
    r.normalize ();
    return r;
  }

  inline ExtendedFloat operator+ (const ExtendedFloat & lhs, const ExtendedFloat & rhs) {
    auto r = denorm_add (lhs, rhs);
    r.normalize ();
    return r;
  }

  inline ExtendedFloat operator- (const ExtendedFloat & lhs) {
    return ExtendedFloat(-lhs.float_part(),lhs.exponent_part());
  }

  inline bool operator== (const ExtendedFloat & lhs, const ExtendedFloat & rhs) {
    return (lhs.float_part()==rhs.float_part() && lhs.exponent_part()==rhs.exponent_part());
  }

  inline double log (const ExtendedFloat & ef) {
    static const auto ln_radix = std::log (static_cast<double> (ExtendedFloat::radix));
    return std::log (ef.float_part ()) + static_cast<double> (ef.exponent_part ()) * ln_radix;
  }

  //!!! no check on the validation of the conversion
  inline double convert(const ExtendedFloat & ef) {
    return ef.float_part () * constexpr_power<double>(ExtendedFloat::radix, ef.exponent_part ());
  }

// TODO add Vector<EF> = Vector<double> + one exp (for lik vectors for one site, big tree case)
// TODO add Vector<EF> = Vector<double> + Vector<exps> (for lik vec by site, eigen, delayed_norm)
} // namespace bpp

#endif // BPP_NEWPHYL_EXTENDEDFLOAT_H
