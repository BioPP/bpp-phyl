//
// File: Range.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-05
// Last modified: 2017-05-05
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

#pragma once
#ifndef BPP_NEWPHYL_RANGE_H
#define BPP_NEWPHYL_RANGE_H

#include <iterator>
#include <type_traits>

namespace bpp {
namespace Iterator {
	template <typename Int> class Integer {
	public:
		using iterator_category = std::random_access_iterator_tag;
		using value_type = Int;
		using difference_type = std::ptrdiff_t;
		using pointer = const value_type *;
		using reference = value_type; // Not standard compliant

		constexpr Integer () noexcept = default;
		constexpr Integer (Int n) noexcept : n_ (n) {}

		// Input / output
		Integer & operator++ () noexcept { return ++n_, *this; }
		constexpr reference operator* () const noexcept { return n_; }
		constexpr pointer operator-> () const noexcept { return &n_; }
		constexpr bool operator== (const Integer & o) const noexcept { return n_ == o.n_; }
		constexpr bool operator!= (const Integer & o) const noexcept { return n_ != o.n_; }

		// Forward
		Integer operator++ (int) noexcept {
			Integer tmp (*this);
			++*this;
			return tmp;
		}

		// Bidir
		Integer & operator-- () noexcept { return --n_, *this; }
		Integer operator-- (int) noexcept {
			Integer tmp (*this);
			--*this;
			return tmp;
		}

		// Random access
		Integer & operator+= (difference_type n) noexcept { return n_ += n, *this; }
		constexpr Integer operator+ (difference_type n) const noexcept { return Integer (n_ + n); }
		friend constexpr Integer operator+ (difference_type n, const Integer & it) noexcept {
			return it + n;
		}
		Integer & operator-= (difference_type n) noexcept { return n_ -= n, *this; }
		constexpr Integer operator- (difference_type n) const noexcept { return Integer (n_ - n); }
		constexpr difference_type operator- (const Integer & o) const noexcept { return n_ - o.n_; }
		constexpr reference operator[] (difference_type n) const noexcept { return n_ + n; }
		constexpr bool operator< (const Integer & o) const noexcept { return n_ < o.n_; }
		constexpr bool operator> (const Integer & o) const noexcept { return n_ > o.n_; }
		constexpr bool operator<= (const Integer & o) const noexcept { return n_ <= o.n_; }
		constexpr bool operator>= (const Integer & o) const noexcept { return n_ >= o.n_; }

	private:
		Int n_{};
	};
}
namespace Range {
	template <typename It> class Base {
		/* Stores a pair of iterators (no iterator requirements except copy).
		 * This can be specialized for iterators that share some data to reduce storage size.
		 * For simplicity, has value semantics and should be passed by value.
		 */
	public:
		constexpr Base (It begin, It end) noexcept : begin_ (begin), end_ (end) {}

		constexpr It begin () const noexcept { return begin_; }
		constexpr It end () const noexcept { return end_; }

	private:
		It begin_;
		It end_;
	};

	template <typename It> class Range : public Base<It> {
		/* Represents a matching pair of iterators.
		 *
		 * This class has constant value semantics.
		 * It can not be modified, but new ranges can be built from it.
		 * It can be passed by value.
		 *
		 * This class has the same properties has the underlying iterator type.
		 * This applies for invalidation, and which operations can be used.
		 * Some operations may work inefficiently for low capability iterators (like size ()).
		 * Beware with less than forward iterators (risk of iterating on them in size()).
		 *
		 * TODO more slices (to_end, from_start).
		 * TODO checked operations ?
		 */
	public:
		using ValueType = typename std::iterator_traits<It>::value_type;
		using ReferenceType = typename std::iterator_traits<It>::reference;
		using DifferenceType = typename std::iterator_traits<It>::difference_type;
		using IteratorCategory = typename std::iterator_traits<It>::iterator_category;

		// Basic constructors
		constexpr Range (Base<It> base) noexcept : Base<It> (base) {}

		using Base<It>::begin;
		using Base<It>::end;

		// Input/forward iterator
		constexpr bool empty () const { return begin () == end (); }
		constexpr ReferenceType front () const { return *begin (); }
		Range pop_front () const { return Base<It>{std::next (begin ()), end ()}; }

		// Bidirectional iterator
		ReferenceType back () const { return *std::prev (end ()); }
		Range pop_back () const { return Base<It>{begin (), std::prev (end ())}; }

		// Random access iterator
		DifferenceType size () const { return std::distance (begin (), end ()); }
		ReferenceType operator[] (DifferenceType n) const { return begin ()[n]; }
		Range pop_front (DifferenceType n) const { return Base<It>{std::next (begin (), n), end ()}; }
		Range pop_back (DifferenceType n) const { return Base<It>{begin (), std::prev (end (), n)}; }

		// Interval-like API
		constexpr bool contains (It it) const { return begin () <= it && it < end (); }
		DifferenceType offset_of (It it) const { return std::distance (begin (), it); }

		// "nicer" api (python like slice ; but at(size ()) return end ())
		// TODO improve...
		It at (DifferenceType n) const {
			auto index = n < 0 ? n + size () : n;
			return std::next (begin (), index);
		}
		Range slice (DifferenceType from, DifferenceType to) const {
			return Base<It>{at (from), at (to)};
		}
		Range slice_to (DifferenceType to) const { return Base<It>{begin (), at (to)}; }
		Range slice_from (DifferenceType from) const { return Base<It>{at (from), end ()}; }

		// Build a container from this range
		template <typename Container> Container to_container () const {
			return Container (begin (), end ());
		}
	};

	// Factory functions

	// From iterator pair (if it is considered an iterator by STL).
	template <typename It, typename = typename std::iterator_traits<It>::iterator_category>
	Range<It> range (It begin, It end) {
		return Base<It>{begin, end};
	}

	// From container (enabled if it supports a begin() function and is lvalue).
	namespace Detail {
		// Utility function to extract the type returned by begin(), in the right context.
		using std::begin;
		template <typename T> auto call_begin (const T & t) -> decltype (begin (t));
		template <typename T> auto call_begin (T & t) -> decltype (begin (t));
	}
	template <typename Container,
	          typename It = decltype (Detail::call_begin (std::declval<Container> ())),
	          typename = typename std::enable_if<std::is_lvalue_reference<Container>::value>::type>
	Range<It> range (Container && container) {
		using std::begin;
		using std::end;
		return range (begin (container), end (container));
	}

	// From integers.
	template <typename Int, typename = typename std::enable_if<std::is_integral<Int>::value>::type>
	Range<Iterator::Integer<Int>> range (Int from, Int to) {
		return range (Iterator::Integer<Int>{from}, Iterator::Integer<Int>{to});
	}
	template <typename Int, typename = typename std::enable_if<std::is_integral<Int>::value>::type>
	Range<Iterator::Integer<Int>> range (Int to) {
		return range (Int{0}, to);
	}

	// Other factory functions
	template <typename Container, typename SizeType = decltype (std::declval<Container> ().size ())>
	Range<Iterator::Integer<SizeType>> index_range (const Container & container) {
		return range (container.size ());
	}
}
using Range::range;
}

#endif // BPP_NEWPHYL_RANGE_H
