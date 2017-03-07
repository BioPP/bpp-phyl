// File: ForRange.h
// Authors:
//   Francois Gindraud (2017)
// Created: 21/02/2017

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

/** @file ForRange.h
 * Defines range-related classes and function.
 * Instances of the ForRange<T> class can be used with range-for constructs:
 * @code{.cpp}
 * for (auto i : makeRange (size_t (42)))
 * @endcode
 *
 * This class represents the same concept as bpp::Range.
 * However it cannot use bpp::Range as we need begin() and end() to return iterators.
 * bpp::Range defines begin() and end() to return the bounds...
 */

#include <iterator>
#include <type_traits>

namespace bpp
{
  /** Class that represents an interval of values that can be iterated over using the for-range syntax.
   * Represent the ordered set of values [begin, end[.
   * As constructors have no parameter type deduction, the type must be explicit.
   *
   * This class is intended to be used with numeric-like T types.
   * The type T must at least be copyable.
   *
   * @see makeRange
   */
  template <typename T>
  class ForRange
  {
    static_assert(std::is_copy_constructible<T>::value, "bpp::ForRange<T>: T must be copyable");

  private:
    T begin_;
    T end_;

  public:
    ForRange(const T& beginElem, const T& endElem)
      : begin_(beginElem)
      , end_(endElem)
    {
    }

    ForRange(const ForRange&) = default;
    ForRange& operator=(const ForRange&) = default;

    /// Get lower bound.
    T& rangeBegin(void) { return begin_; }
    /// Get lower bound.
    const T& rangeBegin(void) const { return begin_; }
    /// Get upper bound.
    T& rangeEnd(void) { return end_; }
    /// Get upper bound.
    const T& rangeEnd(void) const { return end_; }

    /** Internal iterator class.
     * It is an input iterator (can only read values and increment).
     * It only stores a value of T, returns it on dereference, or increment it.
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

  /** Build a range with type deduction.
   * Returns a range [beginElem, endElem[.
   * Use this function instead of calling the constructor (which requires an explicit type).
   */
  template <typename T>
  ForRange<T> makeRange(const T& beginElem, const T& endElem)
  {
    return {beginElem, endElem};
  }

  /** Build a range with type deduction and default begin.
   * Returns a range [T(), endElem[.
   * For integers T() == 0.
   * Use this function instead of calling the constructor (which requires an explicit type).
   */
  template <typename T>
  ForRange<T> makeRange(const T& endElem)
  {
    return {T(), endElem};
  }
} // end namespace bpp
#endif // _FORRANGE_H_
