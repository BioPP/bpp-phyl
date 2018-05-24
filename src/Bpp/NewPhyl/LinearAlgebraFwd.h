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

/** @file Forward declaration of Eigen Vector/Matrix types.
 */

namespace bpp {
/** @name Wrappers to Eigen vector/matrix types.
 *
 * TODO doc foward decl scheme here
 */
class VectorDouble;
class MatrixDouble;

/** Store a dimension for type T.
 * Declared but undefined by default.
 *
 * Specialisations should be defined in the same header declaring the T type.
 */
template <typename T> class Dimension;

/// @name Specialisation of Dimension<T> for double, VectorDouble, MatrixDouble.
///@{
template <> class Dimension<double> {
	// Empty
};
template <> class Dimension<VectorDouble> {
public:
	SizeType size{};
	constexpr Dimension () = default;
	constexpr Dimension (SizeType size_) noexcept : size (size_) {}
};
template <> class Dimension<MatrixDouble> {
public:
	SizeType rows{};
	SizeType cols{};
	constexpr Dimension () = default;
	constexpr Dimension (SizeType rows_, SizeType cols_) noexcept : rows (rows_), cols (cols_) {}
};
///@}

/// @name Enable comparisons.
///@{
constexpr bool operator== (const Dimension<VectorDouble> & lhs,
                           const Dimension<VectorDouble> & rhs) noexcept {
	return lhs.size == rhs.size;
}
constexpr bool operator!= (const Dimension<VectorDouble> & lhs,
                           const Dimension<VectorDouble> & rhs) noexcept {
	return !(lhs == rhs);
}
constexpr bool operator== (const Dimension<MatrixDouble> & lhs,
                           const Dimension<MatrixDouble> & rhs) noexcept {
	return lhs.rows == rhs.rows && lhs.cols == rhs.cols;
}
constexpr bool operator!= (const Dimension<MatrixDouble> & lhs,
                           const Dimension<MatrixDouble> & rhs) noexcept {
	return !(lhs == rhs);
}
///@}

/// @name Specialisations of to_string (pretty printing).
///@{
std::string to_string (const Dimension<VectorDouble> & dim);
std::string to_string (const Dimension<MatrixDouble> & dim);
///@}

/// @name Get dimensions from objects.
///@{
Dimension<VectorDouble> dimension (const VectorDouble & v) noexcept;
Dimension<MatrixDouble> dimension (const MatrixDouble & m) noexcept;
///@}

/// @name Get target dimensions from objects.
///@{
inline Dimension<double> targetDimension (const double &) noexcept {
	return {};
}
Dimension<VectorDouble> targetDimension (const VectorDouble & v) noexcept;
Dimension<MatrixDouble> targetDimension (const MatrixDouble & m) noexcept;
///@}

/// @name Set target dimensions of objects.
///@{
inline void setTargetDimension (double &, const Dimension<double> &) noexcept {}
void setTargetDimension (VectorDouble & v, const Dimension<VectorDouble> & dim) noexcept;
void setTargetDimension (MatrixDouble & m, const Dimension<MatrixDouble> & dim) noexcept;
///@}

// FIXME deprecate ?
// Exact predicates (not a fuzzy comparison)
bool isExactZero (const double & d);
bool isExactZero (const VectorDouble &);
bool isExactZero (const MatrixDouble &);
bool isExactOne (const double & d);
bool isExactOne (const VectorDouble &);
bool isExactOne (const MatrixDouble &);
bool isExactIdentity (const MatrixDouble &);

/* TODO accessor.
 * These accessors cannot be inlined, and are thus of low performance.
 * If high performance access is needed, include LinearAlgebra.h (and Eigen) for inlinable access.
 *
 * If values come out of a Data Flow node, you can implement your computation as a data flow node.
 * You will get automatic caching of intermediates, derivation if supported by your node, ...
 */

} // namespace bpp

#endif // BPP_NEWPHYL_LINEARALGEBRAFWD_H
