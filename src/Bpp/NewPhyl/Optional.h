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

// TODO add doc

struct InPlace {
	constexpr InPlace () = default;
};
constexpr InPlace in_place{};

template <typename T> class Optional {
public:
	using value_type = T;

	constexpr Optional () = default;
	Optional (const Optional & other) {
		if (other)
			create (*other);
	}
	Optional (Optional && other) noexcept {
		if (other)
			create (*std::move (other));
	}
	template <typename... Args> Optional (InPlace, Args &&... args) {
		create (std::forward<Args> (args)...);
	}

	~Optional () { reset (); }

	Optional & operator= (const Optional & other) {
		reset ();
		if (other)
			create (*other);
		return *this;
	}
	Optional & operator= (Optional && other) {
		reset ();
		if (other)
			create (std::move (other));
		return *this;
	}
	template <typename U,
	          typename = typename std::enable_if<
	              !std::is_same<Optional, typename std::decay<U>::type>::value>::type>
	Optional & operator= (U && u) {
		reset ();
		create (std::forward<U> (u));
		return *this;
	}

	T & value () & noexcept {
		assert (has_value_);
		return *value_ptr ();
	}
	const T & value () const & noexcept {
		assert (has_value_);
		return *value_ptr ();
	}
	T && value () && noexcept {
		assert (has_value_);
		return std::move (*value_ptr ());
	}
	const T && value () const && noexcept {
		assert (has_value_);
		return std::move (*value_ptr ());
	}

	T * operator-> () noexcept { return value_ptr (); }
	constexpr const T * operator-> () const noexcept { return value_ptr (); }
	T & operator* () & noexcept { return *value_ptr (); }
	constexpr const T & operator* () const & noexcept { return *value_ptr (); }
	T && operator* () && noexcept { return std::move (*value_ptr ()); }
	constexpr const T && operator* () const && noexcept { return std::move (*value_ptr ()); }

	constexpr bool has_value () const noexcept { return has_value_; }
	constexpr operator bool () const noexcept { return has_value_; }

	template <typename U> constexpr T value_or (U && default_value) const & {
		return has_value_ ? *value_ptr () : static_cast<T> (std::forward<U> (default_value));
	}
	template <typename U> T value_or (U && default_value) && {
		return has_value_ ? std::move (*value_ptr ())
		                  : static_cast<T> (std::forward<U> (default_value));
	}

	template <typename Callable> T value_or_generate (Callable && callable) const & {
		return has_value_ ? *value_ptr () : std::forward<Callable> (callable) ();
	}
	template <typename Callable> T value_or_generate (Callable && callable) && {
		return has_value_ ? std::move (*value_ptr ()) : std::forward<Callable> (callable) ();
	}

	void reset () noexcept {
		if (has_value_)
			destroy ();
	}
	template <typename... Args> T & emplace (Args &&... args) {
		reset ();
		create (std::forward<Args> (args)...);
		return *value_ptr ();
	}
	void swap (Optional & other) noexcept {
		using std::swap;
		if (has_value () && other.has_value ()) {
			swap (value (), other.value ());
		} else if (!has_value () && other.has_value ()) {
			create (std::move (*other));
			other.destroy ();
		} else if (has_value () && !other.has_value ()) {
			other.create (std::move (this->value ()));
			destroy ();
		}
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

	T * value_ptr () noexcept { return reinterpret_cast<T *> (&storage_); }
	constexpr const T * value_ptr () const noexcept {
		return reinterpret_cast<const T *> (&storage_);
	}

	typename std::aligned_storage<sizeof (T), alignof (T)>::type storage_;
	bool has_value_{false};
};
}

#endif // BPP_NEWPHYL_OPTIONAL_H
