//
// File: ForRange.h
// Authors:
//   Francois Gindraud (2017)
// Created on: mardi 21 février 2017, à 16h55
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _FORRANGE_H_
#define _FORRANGE_H_

/**
 * @file ForRange.h
 * Defines range-related classes and function.
 * Instances of the ForRange<T> class can be used with range-for constructs:
 * @code{.cpp}
 * for (auto i : make_range (size_t (42)))
 * @endcode
 *
 * FIXME try to make it compatible with bpp::Range ?
 * bpp::Range from bpp-core defines begin() and end() for direct access to bounds.
 * Thus it cannot be used as is for range-for constructs.
 */

#if !(__cplusplus >= 201103L)
#error "Bpp/Phyl/ForRange.h requires C++11 support"
#endif

#include <iterator>
#include <type_traits>

namespace bpp
{
  /**
   * @brief ForRange class.
   *
   * Can be used in for-range constructs.
   * Represent the ordered set of values [begin, end[.
   * As constructors have no parameter type deduction, the type must be explicit.
   * @see make_range
   */
  template <typename T>
  class ForRange
  {
    static_assert(std::is_copy_constructible<T>::value, "bpp::ForRange<T>: T must be copyable");

  private:
    T begin_;
    T end_;

  public:
    ForRange(const T& begin_elem, const T& end_elem)
      : begin_(begin_elem)
      , end_(end_elem)
    {
    }

    ForRange(const ForRange&) = default;
    ForRange& operator=(const ForRange&) = default;

    /**
	  * @brief Access range bounds.
	  */
    T& range_begin(void) { return begin_; }
    const T& range_begin(void) const { return begin_; }
    T& range_end(void) { return end_; }
    const T& range_end(void) const { return end_; }

    /**
     * @brief Internal iterator class.
     *
     * input iterator (can only read values and increment).
     */
    class iterator : public std::iterator<std::input_iterator_tag, T>
    {
    private:
      T v_;

    public:
      explicit iterator(const T& v = T())
        : v_(v)
      {
      }
      iterator& operator++(void)
      {
        ++v_;
        return *this;
      }
      iterator operator++(int)
      {
        iterator c = *this;
        operator++();
        return c;
      }
      bool operator==(iterator other) { return v_ == other.v_; }
      bool operator!=(iterator other) { return !operator==(other); }
      T operator*(void) { return v_; }
    };

    iterator begin(void) const { return iterator(begin_); }
    iterator end(void) const { return iterator(end_); }
  };

  /**
	* @brief Build a range with type deduction.
	*
	* Returns a range [begin_elem, end_elem[.
	*/
  template <typename T>
  ForRange<T> make_range(const T& begin_elem, const T& end_elem)
  {
    return {begin_elem, end_elem};
  }

  /**
	* @brief Build a range with type deduction and default begin.
	*
	* Returns a range [T(), end_elem[.
	* For integers T() == 0.
	*/
  template <typename T>
  ForRange<T> make_range(const T& end_elem)
  {
    return {T(), end_elem};
  }
} // end namespace bpp
#endif // _RANGE_H_
