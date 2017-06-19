//
// File: FrozenPtr.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-06-08
// Last modified: 2017-06-08
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

#ifndef BPP_NEWPHYL_FROZENPTR_H
#define BPP_NEWPHYL_FROZENPTR_H

#include <memory>
#include <type_traits>

namespace bpp {

// Forward declaration
template <typename T> class FreezableUniquePtr;
template <typename T> class FrozenSharedPtr;

template <typename T> class FreezableUniquePtr {
	/* A unique pointer to a mutable ressource.
	 *
	 * Internally it is a shared_ptr:
	 * - allow shared_from_this classes to not be UB
	 * - avoid a separate control block allocation when converting to FrozenSharedPtr.
	 * Thus copy is manually disabled.
	 */
public:
	FreezableUniquePtr () = default;
	FreezableUniquePtr (const FreezableUniquePtr &) = delete;
	FreezableUniquePtr & operator= (const FreezableUniquePtr &) = delete;
	FreezableUniquePtr (FreezableUniquePtr &&) = default;
	FreezableUniquePtr & operator= (FreezableUniquePtr &&) = default;
	~FreezableUniquePtr () = default;

	// Move construct from unique_ptr ; prefer the static make function.
	FreezableUniquePtr (std::unique_ptr<T> && ptr) noexcept : ptr_ (std::move (ptr)) {}

	// Builds "in place", less overhead than unique_ptr constructor.
	template <typename... Args> static FreezableUniquePtr make (Args &&... args) {
		return {std::make_shared<T> (std::forward<Args> (args)...)};
	}

	// Access
	constexpr explicit operator bool () const noexcept { return bool(ptr_); }
	T * get () const noexcept { return ptr_.get (); }
	T & operator* () const noexcept { return *ptr_; }
	T * operator-> () const noexcept { return get (); }

	/* Conversion:
	 * - always a move to ensure it is unique.
	 * - can upcast if convertible.
	 * - can extract the shared_ptr from here (useful for conversion)
	 */
	template <typename U,
	          typename = typename std::enable_if<std::is_convertible<U *, T *>::value>::type>
	explicit FreezableUniquePtr (FreezableUniquePtr<U> && ptr) noexcept
	    : ptr_ (static_cast<std::shared_ptr<U>> (std::move (ptr))) {}
	template <typename U,
	          typename = typename std::enable_if<std::is_convertible<U *, T *>::value>::type>
	FreezableUniquePtr & operator= (FreezableUniquePtr<U> && ptr) noexcept {
		ptr_ = static_cast<std::shared_ptr<U>> (std::move (ptr));
		return *this;
	}
	explicit operator std::shared_ptr<T> () && noexcept { return std::move (ptr_); }

	// Transform to a FrozenSharedPtr
	FrozenSharedPtr<T> freeze () &&;

private:
	// Private to avoid violating the unique contract by passing a shared shared_ptr.
	FreezableUniquePtr (std::shared_ptr<T> && ptr) noexcept : ptr_ (std::move (ptr)) {}

	std::shared_ptr<T> ptr_;
};

// Equivalent to std::make_unique
template <typename T, typename... Args> FreezableUniquePtr<T> make_freezable (Args &&... args) {
	return FreezableUniquePtr<T>::make (std::forward<Args> (args)...);
}
template <typename T> FreezableUniquePtr<T> make_freezable (T && t) {
	return FreezableUniquePtr<T>::make (std::forward<T> (t));
}

template <typename T> class FrozenSharedPtr {
	/* A shared pointer to an immutable ressource.
	 */
public:
	using ConstT = typename std::add_const<T>::type;

	// All move/copy constructor/assignemnt default

	// Move construct from FreezableUniquePtr
	FrozenSharedPtr (FreezableUniquePtr<T> && ptr) noexcept
	    : ptr_ (static_cast<std::shared_ptr<T>> (std::move (ptr))) {}

	// Build in place (skip Freezable pointer)
	template <typename... Args> static FrozenSharedPtr make (Args &&... args) {
		return std::make_shared<ConstT> (std::forward<Args> (args)...);
	}

	// Access
	constexpr explicit operator bool () const noexcept { return bool(ptr_); }
	ConstT * get () const noexcept { return ptr_.get (); }
	ConstT & operator* () const noexcept { return *ptr_; }
	ConstT * operator-> () const noexcept { return get (); }

	// Conversion: can upcast
	template <typename U,
	          typename = typename std::enable_if<std::is_convertible<U *, T *>::value>::type>
	explicit FrozenSharedPtr (const FrozenSharedPtr<U> & ptr) noexcept : ptr_ (ptr.get_shared ()) {}
	template <typename U,
	          typename = typename std::enable_if<std::is_convertible<U *, T *>::value>::type>
	explicit FrozenSharedPtr (FrozenSharedPtr<U> && ptr) noexcept
	    : ptr_ (std::move (ptr).get_shared ()) {}
	template <typename U,
	          typename = typename std::enable_if<std::is_convertible<U *, T *>::value>::type>
	FrozenSharedPtr & operator= (const FrozenSharedPtr<U> & ptr) noexcept {
		ptr_ = ptr.get_shared ();
		return *this;
	}
	template <typename U,
	          typename = typename std::enable_if<std::is_convertible<U *, T *>::value>::type>
	FrozenSharedPtr & operator= (FrozenSharedPtr<U> && ptr) noexcept {
		ptr_ = std::move (ptr).get_shared ();
		return *this;
	}

	// Allow shared_from_this functionnality, may violate the properties !
	// FIXME use a custom base enable_shared_from_this type, and std::is_base_of check
	template <typename SharedFromThisType>
	static FrozenSharedPtr shared_from_this (const SharedFromThisType & t) {
		return FrozenSharedPtr{t.shared_from_this ()};
	}

	// Impl access (considered internal)
	const std::shared_ptr<ConstT> & get_shared () const & noexcept { return ptr_; }
	std::shared_ptr<ConstT> && get_shared () && noexcept { return std::move (ptr_); }

private:
	// Only used by shared_from_this
	FrozenSharedPtr (std::shared_ptr<ConstT> && ptr) noexcept : ptr_ (std::move (ptr)) {}

	std::shared_ptr<ConstT> ptr_;
};

// Similar to std::make_shared
template <typename T, typename... Args> FrozenSharedPtr<T> make_frozen (Args &&... args) {
	return FrozenSharedPtr<T>::make (std::forward<Args> (args)...);
}
template <typename T> FrozenSharedPtr<T> make_frozen (T && t) {
	return FrozenSharedPtr<T>::make (std::forward<T> (t));
}

// Deferred freeze() impl
template <typename T> FrozenSharedPtr<T> FreezableUniquePtr<T>::freeze () && {
	return FrozenSharedPtr<T>{std::move (*this)};
}
}

#endif // BPP_NEWPHYL_FROZENPTR_H
