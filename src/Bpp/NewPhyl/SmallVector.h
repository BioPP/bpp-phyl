//
// File: SmallVector.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-03
// Last modified: 2017-05-03
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
#ifndef BPP_NEWPHYL_SMALLVECTOR_H
#define BPP_NEWPHYL_SMALLVECTOR_H

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <limits>
#include <type_traits>
#include <utility>

namespace bpp {

namespace Detail {
	template <typename T> void call_destructors (T * a, std::size_t n) {
		for (std::size_t i = 0; i < n; ++i)
			a[i].~T ();
	}
}

template <typename T, typename Allocator = std::allocator<T>>
class SmallVectorBase : private Allocator {
  /* Split SmallVector code:
   * - SmallVectorBase has most of the API, and can grow but not shrink
   *   If created as a SmallVector<T, N>, data_ starts as ptr to T[N] buf.
   *   capacity_ can grow (replacing the buf with allocated data)
   *   it cannot shrink and reuse the local buf (no way to retrieve it from Base).
   *
   * - SmallVector<T, N> has the inline buf only and gives it to Base at init
   *
   * Copy and moves can be optimal between SmallVector.
   *
   * TODO use allocator_traits
   *
   * FIXME need a bit flag to tell if allocated
   * LLVM relies on object layout to do its trick...
   * Use bitfield ?
   * Store size and cap as uint16 ?
   */
private:
	std::size_t size_;
	std::size_t capacity_;
	T * data_;
};

template <typename T, std::size_t threshold_ = sizeof (T *) / sizeof (T),
          typename Allocator = std::allocator<T>>
class SmallVector : private Allocator {
	// capacity_ is never 0
	// TODO will need to review exception behavior
private:
	static_assert (threshold_ > 0, "inline data threshold must be positive");
	using internal_size_type = std::uint32_t;

public:
	using value_type = T;
	using allocator_type = Allocator;
	using size_type = std::size_t;
	using difference_type = std::ptrdiff_t;
	using reference = value_type &;
	using const_reference = const value_type &;
	using pointer = value_type *;
	using const_pointer = const value_type *;
	using iterator = pointer;
	using const_iterator = const_pointer;
	using reverse_iterator = std::reverse_iterator<iterator>;
	using reverse_const_iterator = std::reverse_iterator<const_iterator>;

	static_assert (sizeof (size_type) >= sizeof (internal_size_type),
	               "size_type is smaller than internal_size_type");

	// Basic TODO
	~SmallVector () {
		Detail::call_destructors (data (), size_);
		if (is_allocated ())
			Allocator::deallocate (allocatedData_, capacity_);
	}
	allocator_type get_allocator () const { return *this; }

	// Element access
	reference at (size_type pos) {
		return pos < size_ ? operator[] (pos) : throw std::out_of_range{"at()"};
	}
	constexpr const_reference at (size_type pos) const {
		return pos < size_ ? operator[] (pos) : throw std::out_of_range{"at()"};
	}
	reference operator[] (size_type pos) noexcept { return *(data () + pos); }
	constexpr const_reference operator[] (size_type pos) const noexcept { return *(data () + pos); }
	reference front () noexcept { return *data (); }
	constexpr const_reference front () const noexcept { return *data (); }
	reference back () noexcept { return *(data () + (size_ - 1)); }
	constexpr const_reference back () const noexcept { return *(data () + (size_ - 1)); }
	T * data () noexcept {
		return is_allocated () ? allocatedData_ : reinterpret_cast<T *> (&inlineData_);
	}
	constexpr const T * data () const noexcept {
		return is_allocated () ? allocatedData_ : reinterpret_cast<const T *> (&inlineData_);
	}

	// Iterators
	iterator begin () noexcept { return data (); }
	constexpr const_iterator begin () const noexcept { return data (); }
	constexpr const_iterator cbegin () const noexcept { return data (); }
	iterator end () noexcept { return data () + size_; }
	constexpr const_iterator end () const noexcept { return data () + size_; }
	constexpr const_iterator cend () const noexcept { return data () + size_; }
	reverse_iterator rbegin () noexcept { return end (); }
	constexpr reverse_const_iterator rbegin () const noexcept { return end (); }
	constexpr reverse_const_iterator crbegin () const noexcept { return end (); }
	reverse_iterator rend () noexcept { return begin (); }
	constexpr reverse_const_iterator rend () const noexcept { return begin (); }
	constexpr reverse_const_iterator crend () const noexcept { return begin (); }

	// Capacity
	constexpr bool empty () const noexcept { return size_ == 0; }
	constexpr size_type size () const noexcept { return size_; }
	constexpr size_type max_size () const noexcept {
		return std::numeric_limits<internal_size_type>::max ();
	}
	void reserve (size_type new_cap) {
		if (new_cap > capacity_) {
			// Reserve just moves data to the exact required capacity
			move_to_new_allocated_storage (new_cap);
		}
	}
	constexpr size_type capacity () const noexcept { return capacity_; }
	void shrink_to_fit (); // TODO weird guy, might revert to inline storage

	// Modifiers TODO push_back with x2 growing scheme
	void clear () noexcept {
		call_destructors (data (), size_);
		size_ = 0;
	}
	void push_back (const T & value) {
		auto * storage = size_ < capacity_ ? data () : move_to_new_allocated_storage (capacity_ * 2);
		new (storage + size_) T (value);
		++size_;
	}
	void push_back (T && value) {
		auto * storage = size_ < capacity_ ? data () : move_to_new_allocated_storage (capacity_ * 2);
		new (storage + size_) T (std::move (value));
		++size_;
	}
	template <typename... Args> reference emplace_back (Args &&... args) {
		auto * storage = size_ < capacity_ ? data () : move_to_new_allocated_storage (capacity_ * 2);
		auto * object = new (storage + size_) T (std::forward<Args> (args)...);
		++size_;
		return *object;
	}

	// Small vector specific information
	static constexpr size_type threshold () noexcept { return threshold_; }
	constexpr bool is_allocated () const noexcept { return capacity_ > threshold_; }

private:
	// This function creates a new allocated storage and relocates the current data to it.
	T * move_to_new_allocated_storage (size_type new_cap) {
		// Create new storage
		T * new_storage = Allocator::allocate (new_cap);
		// Move from old, and call the destructors on moved-from T
		T * d = data ();
		std::move (d, d + size_, new_storage);
		Detail::call_destructors (d, size_);
		// Swap and delete old if was allocated
		if (is_allocated ())
			Allocator::deallocate (allocatedData_, capacity_);
		capacity_ = new_cap;
		allocatedData_ = new_storage;
		return new_storage;
	}

	internal_size_type size_{0};
	internal_size_type capacity_{threshold_};
	union {
		typename std::aligned_storage<threshold_ * sizeof (T), alignof (T)>::type inlineData_;
		T * allocatedData_;
	};
};
// Non member functions TODO
}

#endif // BPP_NEWPHYL_SMALLVECTOR_H
