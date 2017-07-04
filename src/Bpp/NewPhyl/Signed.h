//
// File: Signed.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-07-04 00:00:00
// Last modified: 2017-07-04
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

#ifndef BPP_NEWPHYL_SIGNED_H
#define BPP_NEWPHYL_SIGNED_H

#include <cassert>
#include <cstddef>
#include <vector>

namespace bpp {

// Use signed values for indexes and sizes instead of unsigned (CppCoreGuidelines)
using SizeType = std::ptrdiff_t;
using IndexType = std::ptrdiff_t;

template <typename T> class Vector {
	// Vector with signed indexing
private:
	using Container = std::vector<T>;

public:
	using size_type = SizeType;
	using reference = typename Container::reference;
	using const_reference = typename Container::const_reference;
	using iterator = typename Container::iterator;
	using const_iterator = typename Container::const_iterator;

	Vector () : vec_ () {}
	explicit Vector (size_type size) : vec_ (static_cast<typename Container::size_type> (size)) {}

	reference at (size_type i) {
		assert (0 <= i);
		return vec_.at (static_cast<typename Container::size_type> (i));
	}
	const_reference at (size_type i) const {
		assert (0 <= i);
		return vec_.at (static_cast<typename Container::size_type> (i));
	}

	reference operator[] (size_type i) {
		assert (0 <= i);
		assert (i < size ());
		return vec_[static_cast<typename Container::size_type> (i)];
	}
	const_reference operator[] (size_type i) const {
		assert (0 <= i);
		assert (i < size ());
		return vec_[static_cast<typename Container::size_type> (i)];
	}

	iterator begin () noexcept { return vec_.begin (); }
	const_iterator begin () const noexcept { return vec_.begin (); }
	iterator end () noexcept { return vec_.end (); }
	const_iterator end () const noexcept { return vec_.end (); }

	size_type size () const noexcept { return static_cast<size_type> (vec_.size ()); }

private:
	Container vec_;
};
}
#endif // BPP_NEWPHYL_SIGNED_H
