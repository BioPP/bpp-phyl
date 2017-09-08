//
// File: Optional.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-22
// Last modified: 2017-05-22
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
#ifndef BPP_NEWPHYL_OPTIONAL_H
#define BPP_NEWPHYL_OPTIONAL_H

#include <cassert>
#include <type_traits>
#include <utility>

namespace bpp {

// Type tag and global value for in_place construction
// TODO put in its own header when SmallVector arrives to NewPhyl
template <typename T> struct InPlace { constexpr InPlace () = default; };
constexpr InPlace<void> in_place{};

// Type tag and global value for "No value"
struct NullOpt {
	constexpr NullOpt () = default;
};
constexpr NullOpt nullopt{};

template <typename T> class Optional {
	/* Optional<T> is similar to the c++17 std::optional<T>.
	 * It conditionally stores in place a T value.
	 *
	 * It behaves like a pointer to data, also for constness:
	 * A const Optional cannot change its state ; however the stored value is mutable.
	 * An Optional<const T> can change its state but not the stored value after creation.
	 *
	 * Moves are considered to change state, so move overloads require a non const && *this.
	 * A moved-from Optional still contains a value (which is moved-from).
	 *
	 * Types that do not support copy/move assignments use destruction+copy/move construction instead.
	 * Non copyable and non movable object should only use emplace().
	 */
	static_assert (!std::is_reference<T>::value, "base Optional<T> does not support references");

public:
	using value_type = T;

	// Constructors
	constexpr Optional () : has_value_ (false) {}
	constexpr Optional (NullOpt) noexcept : Optional () {}
	Optional (const Optional & other) : Optional () {
		if (other)
			create (*other);
	}
	Optional (Optional && other) noexcept : Optional () {
		if (other)
			create (*std::move (other));
	}
	Optional (const T & t) : Optional () { create (t); }
	Optional (T && t) : Optional () { create (std::move (t)); }
	template <typename... Args> Optional (InPlace<void>, Args &&... args) : Optional () {
		create (std::forward<Args> (args)...);
	}
	template <typename U, typename... Args>
	Optional (InPlace<void>, std::initializer_list<U> ilist, Args &&... args) : Optional () {
		create (std::move (ilist), std::forward<Args> (args)...);
	}

	~Optional () { reset (); }

	// Assignment
	Optional & operator= (NullOpt) noexcept {
		reset ();
		return *this;
	}
	Optional & operator= (const Optional & other) {
		if (has_value () && other)
			replace_value_with (*other);
		else if (!has_value () && other)
			create (*other);
		else if (has_value () && !other)
			destroy ();
		return *this;
	}
	Optional & operator= (Optional && other) noexcept {
		if (has_value () && other)
			replace_value_with (std::move (*other));
		else if (!has_value () && other)
			create (std::move (*other));
		else if (has_value () && !other)
			destroy ();
		return *this;
	}
	Optional & operator= (const T & t) {
		if (has_value ())
			replace_value_with (t);
		else
			create (t);
		return *this;
	}
	Optional & operator= (T && t) noexcept {
		if (has_value ())
			replace_value_with (std::move (t));
		else
			create (std::move (t));
		return *this;
	}

	// Status
	constexpr bool has_value () const noexcept { return has_value_; }
	constexpr explicit operator bool () const noexcept { return has_value (); }

	// Unchecked access
	T & value () const & noexcept {
		assert (has_value ());
		return *value_ptr ();
	}
	T && value () && noexcept {
		assert (has_value ());
		return std::move (*value_ptr ());
	}
	T * operator-> () const noexcept { return value_ptr (); }
	T & operator* () const & noexcept { return value (); }
	T && operator* () && noexcept { return std::move (*this).value (); }

	// Modifiers
	void reset () noexcept {
		if (has_value ())
			destroy ();
	}
	template <typename... Args> T & emplace (Args &&... args) {
		reset ();
		create (std::forward<Args> (args)...);
		return value ();
	}
	template <typename U, typename... Args>
	T & emplace (std::initializer_list<U> ilist, Args &&... args) {
		reset ();
		create (std::move (ilist), std::forward<Args> (args)...);
		return value ();
	}
	void swap (Optional & other) noexcept {
		using std::swap;
		if (has_value () && other) {
			swap (value (), other.value ());
		} else if (!has_value () && other) {
			create (std::move (*other));
			other.destroy ();
		} else if (has_value () && !other) {
			other.create (std::move (value ()));
			destroy ();
		}
	}

	// Access with default
	template <typename U> constexpr T value_or (U && default_value) const & {
		return has_value () ? value () : static_cast<T> (std::forward<U> (default_value));
	}
	template <typename U> T value_or (U && default_value) && {
		return has_value () ? std::move (*this).value ()
		                    : static_cast<T> (std::forward<U> (default_value));
	}

	// Access with generated default
	template <typename Callable> T value_or_generate (Callable && callable) const & {
		return has_value () ? value () : std::forward<Callable> (callable) ();
	}
	template <typename Callable> T value_or_generate (Callable && callable) && {
		return has_value () ? std::move (*this).value () : std::forward<Callable> (callable) ();
	}

	// Map : Optional<T> -> Optional<U> with f : T -> U
	template <typename Callable,
	          typename ReturnType = typename std::result_of<Callable (const T &)>::type>
	Optional<ReturnType> map (Callable && callable) const & {
		if (has_value ())
			return std::forward<Callable> (callable) (value ());
		else
			return {};
	}
	template <typename Callable, typename ReturnType = typename std::result_of<Callable (T &&)>::type>
	Optional<ReturnType> map (Callable && callable) && {
		if (has_value ())
			return std::forward<Callable> (callable) (std::move (*this).value ());
		else
			return {};
	}

	// Filter : Optional<T> -> Optional<T>, propagate value if p(value), or return NullOpt
	template <typename Predicate> Optional filter (Predicate && predicate) const & {
		if (has_value () && std::forward<Predicate> (predicate) (value ()))
			return *this;
		else
			return {};
	}
	template <typename Predicate> Optional filter (Predicate && predicate) && {
		if (has_value () && std::forward<Predicate> (predicate) (value ()))
			return std::move (*this);
		else
			return {};
	}

	// Cast : Optional<T> -> Optional<U>.
	template <typename U> Optional<U> cast () const & {
		if (has_value ())
			return static_cast<U> (value ());
		else
			return {};
	}
	template <typename U> Optional<U> cast () && {
		if (has_value ())
			return static_cast<U> (std::move (*this).value ());
		else
			return {};
	}

private:
	template <typename... Args> void create (Args &&... args) {
		assert (!has_value_);
		new (&storage_) T (std::forward<Args> (args)...);
		has_value_ = true;
	}
	void destroy () noexcept {
		assert (has_value_);
		value_ptr ()->~T ();
		has_value_ = false;
	}

	// Implement replace using assignment operators if possible, or destroy+constructor
	template <typename U> void replace_value_with (U && u) {
		replace_value_with_helper (std::forward<U> (u), std::is_assignable<T, U>{});
	}
	template <typename U> void replace_value_with_helper (U && u, std::true_type) {
		value () = std::forward<U> (u);
	}
	template <typename U> void replace_value_with_helper (U && u, std::false_type) {
		destroy ();
		create (std::forward<U> (u));
	}

	// Storage is "mutable" to support muting the object if the optional is const
	T * value_ptr () const noexcept { return reinterpret_cast<T *> (&storage_); }
	mutable typename std::aligned_storage<sizeof (T), alignof (T)>::type storage_;
	bool has_value_;
};

template <typename T> class Optional<T &> {
	/* Specialisation for references.
	 * This specialisation has the same semantics as a pointer.
	 *
	 * Copy and moves just copy the pointer/reference.
	 * Access returns the reference itself, so operator-> will modify the referenced object.
	 * The reference can be modified by assignment and the same modifiers as the main class.
	 * All modifiers bind to both T & and T * for convenience.
	 */
public:
	using value_type = T &;

	// Constructors
	constexpr Optional () = default;
	constexpr Optional (NullOpt) noexcept : Optional () {}
	constexpr Optional (const Optional &) = default;
	constexpr Optional (Optional &&) = default;
	constexpr Optional (std::nullptr_t) noexcept : Optional () {}
	constexpr Optional (T * t) noexcept : pointer_ (t) {}
	constexpr Optional (T & t) noexcept : Optional (&t) {}
	~Optional () = default;

	// Assignment
	Optional & operator= (NullOpt) noexcept {
		reset ();
		return *this;
	}
	Optional & operator= (const Optional &) = default;
	Optional & operator= (Optional &&) = default;
	Optional & operator= (std::nullptr_t) noexcept {
		pointer_ = nullptr;
		return *this;
	}
	Optional & operator= (T * t) {
		pointer_ = t;
		return *this;
	}
	Optional & operator= (T & t) { return *this = &t; }

	// Status
	constexpr bool has_value () const noexcept { return pointer_ != nullptr; }
	constexpr explicit operator bool () const noexcept { return has_value (); }

	// Unchecked access
	T & value () const noexcept {
		assert (has_value ());
		return *pointer_;
	}
	T * operator-> () const noexcept { return pointer_; }
	T & operator* () const noexcept { return value (); }

	// Modifiers
	void reset () noexcept { pointer_ = nullptr; }
	void emplace (T * t) noexcept { pointer_ = t; }
	void emplace (T & t) noexcept { emplace (&t); }
	void swap (Optional & other) noexcept {
		using std::swap;
		swap (pointer_, other.pointer_);
	}

	// Access with default
	constexpr T & value_or (T & default_ref) const { return has_value () ? value () : default_ref; }

	// Access with generated default
	template <typename Callable> constexpr T & value_or_generate (Callable && callable) const {
		return has_value () ? value () : std::forward<Callable> (callable) ();
	}

	// Map : Optional<T &> -> Optional<U> with f : T & -> U
	template <typename Callable, typename ReturnType = typename std::result_of<Callable (T &)>::type>
	Optional<ReturnType> map (Callable && callable) const {
		if (has_value ())
			return std::forward<Callable> (callable) (value ());
		else
			return {};
	}

	// Filter : Optional<T &> -> Optional<T &>, propagate ref if p(ref), or return NullOpt
	template <typename Predicate> Optional filter (Predicate && predicate) const {
		if (has_value () && std::forward<Predicate> (predicate) (value ()))
			return *this;
		else
			return {};
	}

	// Cast : Optional<T &> -> Optional<U>
	template <typename U> Optional<U> cast () const {
		if (has_value ())
			return static_cast<U> (value ());
		else
			return {};
	}

private:
	T * pointer_{nullptr};
};

/* a | b -> returns a if a, or b.
 * Both a and b are evaluated (no operator|| lazy semantics).
 * If both a and b are the same kind of reference (lvalue/rvalue), just return the reference.
 * If reference kind differs, an intermediate Optional<T> is built from the right refs.
 */
template <typename T>
const Optional<T> & operator| (const Optional<T> & a, const Optional<T> & b) noexcept {
	return a ? a : b;
}
template <typename T> Optional<T> && operator| (Optional<T> && a, Optional<T> && b) noexcept {
	return a ? std::move (a) : std::move (b);
}
template <typename T> Optional<T> operator| (const Optional<T> & a, Optional<T> && b) {
	return a ? a : std::move (b);
}
template <typename T> Optional<T> operator| (Optional<T> && a, const Optional<T> & b) {
	return a ? std::move (a) : b;
}
// Overloads ending with plain object ("true" optional)
template <typename T> T operator| (const Optional<T> & a, const T & b) {
	return a.value_or (b);
}
template <typename T> T operator| (Optional<T> && a, const T & b) {
	return std::move (a).value_or (b);
}
template <typename T> T operator| (const Optional<T> & a, T && b) {
	return a.value_or (std::move (b));
}
template <typename T> T operator| (Optional<T> && a, T && b) {
	return std::move (a).value_or (std::move (b));
}

/* Call find(key) on an associative container, return the result as an optional.
 * Intended to avoid comparing the iterator with end().
 */
template <typename Container, typename Key>
Optional<typename Container::mapped_type &> optional_find (Container & container, const Key & key) {
	auto it = container.find (key);
	if (it != container.end ())
		return it->second;
	else
		return {};
}
template <typename Container, typename Key>
Optional<const typename Container::mapped_type &> optional_find (const Container & container,
                                                                 const Key & key) {
	auto it = container.find (key);
	if (it != container.end ())
		return it->second;
	else
		return {};
}
}

#endif // BPP_NEWPHYL_OPTIONAL_H
