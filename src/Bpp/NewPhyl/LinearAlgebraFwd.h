//
// File: LinearAlgebraFwd.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-10-24
// Last modified: 2017-10-24
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

#ifndef BPP_NEWPHYL_LINEARALGEBRAFWD_H
#define BPP_NEWPHYL_LINEARALGEBRAFWD_H

#include <Bpp/NewPhyl/Signed.h>
#include <string>

namespace bpp {
// Forward declarations of eigen types
class VectorDouble;
class MatrixDouble;

// Dimensions types
struct NoDimension {};
struct VectorDimension {
	SizeType size;
	constexpr VectorDimension (SizeType size_) noexcept : size (size_) {}
};
struct MatrixDimension {
	SizeType rows;
	SizeType cols;
	constexpr MatrixDimension (SizeType rows_, SizeType cols_) noexcept
	    : rows (rows_), cols (cols_) {}
};

// Comparisons
constexpr bool operator== (const NoDimension &, const NoDimension &) noexcept {
	return true;
}
constexpr bool operator!= (const NoDimension &, const NoDimension &) noexcept {
	return false;
}
constexpr bool operator== (const VectorDimension & lhs, const VectorDimension & rhs) noexcept {
	return lhs.size == rhs.size;
}
constexpr bool operator!= (const VectorDimension & lhs, const VectorDimension & rhs) noexcept {
	return !(lhs == rhs);
}
constexpr bool operator== (const MatrixDimension & lhs, const MatrixDimension & rhs) noexcept {
	return lhs.rows == rhs.rows && lhs.cols == rhs.cols;
}
constexpr bool operator!= (const MatrixDimension & lhs, const MatrixDimension & rhs) noexcept {
	return !(lhs == rhs);
}

// Pretty printing (overloading std::to_string out of namespace std::)
std::string to_string (const NoDimension &);
std::string to_string (const VectorDimension & dim);
std::string to_string (const MatrixDimension & dim);

// Get dimensions (overloaded)
NoDimension dimensions (const double &) noexcept;
VectorDimension dimensions (const VectorDouble & v) noexcept;
MatrixDimension dimensions (const MatrixDouble & m) noexcept;

// Get dimensions for DF nodes
namespace DF {
	template <typename T> class Value;
}
NoDimension dimensions (const DF::Value<double> &) noexcept;
VectorDimension dimensions (const DF::Value<VectorDouble> & node) noexcept;
MatrixDimension dimensions (const DF::Value<MatrixDouble> & node) noexcept;

// Exact predicates (not a fuzzy comparison)
bool isExactZero (const double & d);
bool isExactZero (const VectorDouble &);
bool isExactZero (const MatrixDouble &);
bool isExactOne (const double & d);
bool isExactOne (const VectorDouble &);
bool isExactOne (const MatrixDouble &);
bool isExactIdentity (const MatrixDouble &);

} // namespace bpp

#endif // BPP_NEWPHYL_LINEARALGEBRAFWD_H
