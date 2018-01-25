//
// File: IntegerRange.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-05 00:00:00
// Last modified: 2017-12-15
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

#ifndef BPP_NEWPHYL_INTEGERRANGE_H
#define BPP_NEWPHYL_INTEGERRANGE_H

#include <iterator>
#include <type_traits>

namespace bpp {
/** @brief Iterator representing an integer value in the integer space.
 *
 * Used to define integer ranges.
 */
template <typename Int> class IntegerIterator {
	static_assert (std::is_integral<Int>::value, "IntegerIterator<Int>: Int must be an integer type");

public:
	using iterator_category = std::random_access_iterator_tag;
	using value_type = Int;
	using difference_type = std::ptrdiff_t;
	using pointer = const value_type *;
	using reference = value_type; // force copies in templated code (algorithm)

	constexpr IntegerIterator () noexcept = default;
	constexpr IntegerIterator (Int n) noexcept : n_ (n) {}

	// Input / output
	IntegerIterator & operator++ () noexcept { return ++n_, *this; }
	constexpr reference operator* () const noexcept { return n_; }
	constexpr pointer operator-> () const noexcept { return &n_; }
	constexpr bool operator== (const IntegerIterator & o) const noexcept { return n_ == o.n_; }
	constexpr bool operator!= (const IntegerIterator & o) const noexcept { return n_ != o.n_; }

	// Forward
	IntegerIterator operator++ (int) noexcept {
		IntegerIterator tmp (*this);
		++*this;
		return tmp;
	}

	// Bidir
	IntegerIterator & operator-- () noexcept { return --n_, *this; }
	IntegerIterator operator-- (int) noexcept {
		IntegerIterator tmp (*this);
		--*this;
		return tmp;
	}

	// Random access
	IntegerIterator & operator+= (difference_type n) noexcept { return n_ += n, *this; }
	constexpr IntegerIterator operator+ (difference_type n) const noexcept {
		return IntegerIterator (n_ + n);
	}
	friend constexpr IntegerIterator operator+ (difference_type n,
	                                            const IntegerIterator & it) noexcept {
		return it + n;
	}
	IntegerIterator & operator-= (difference_type n) noexcept { return n_ -= n, *this; }
	constexpr IntegerIterator operator- (difference_type n) const noexcept {
		return IntegerIterator (n_ - n);
	}
	constexpr difference_type operator- (const IntegerIterator & o) const noexcept {
		return n_ - o.n_;
	}
	constexpr reference operator[] (difference_type n) const noexcept { return n_ + n; }
	constexpr bool operator< (const IntegerIterator & o) const noexcept { return n_ < o.n_; }
	constexpr bool operator> (const IntegerIterator & o) const noexcept { return n_ > o.n_; }
	constexpr bool operator<= (const IntegerIterator & o) const noexcept { return n_ <= o.n_; }
	constexpr bool operator>= (const IntegerIterator & o) const noexcept { return n_ >= o.n_; }

private:
	Int n_{};
};

/** @brief Stores a pair of iterator representing a range.
 *
 * Instead of std::pair, iterators are accessed using begin/end methods.
 */
template <typename Iterator> class IteratorPair {
public:
	constexpr IteratorPair (Iterator beginIt, Iterator endIt) noexcept
	    : begin_ (beginIt), end_ (endIt) {}

	constexpr Iterator begin () const noexcept { return begin_; }
	constexpr Iterator end () const noexcept { return end_; }

private:
	Iterator begin_;
	Iterator end_;
};

/// Builds the integer range [from, to[.
template <typename Int> IteratorPair<IntegerIterator<Int>> range (Int from, Int to) {
	return {from, to};
}

/// Builds the integer range [0, to[.
template <typename Int> IteratorPair<IntegerIterator<Int>> range (Int to) {
	return {Int{0}, to};
}
} // namespace bpp
#endif // BPP_NEWPHYL_INTEGERRANGE_H
