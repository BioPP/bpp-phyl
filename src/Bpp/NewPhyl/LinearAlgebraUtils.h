//
// File: LinearAlgebraUtils.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2018-01-24
// Last modified: 2018-01-24
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

#ifndef BPP_NEWPHYL_LINEARALGEBRAUTILS_H
#define BPP_NEWPHYL_LINEARALGEBRAUTILS_H

#include <Bpp/NewPhyl/LinearAlgebra.h>

namespace bpp {
///@{
/** Create a value filled with the given constant d.
 * Versions for double, VectorDouble, MatrixDouble.
 * The right overload is selected by the Dimension<T> argument.
 * This argument also provide a size for Vector and Matrix variants.
 * This overloaded functions is used in template code (common interface for 3 types).
 */
inline double linearAlgebraValueFilledWith (const Dimension<double> &, double d) {
	return d;
}
inline auto linearAlgebraValueFilledWith (const Dimension<VectorDouble> & dim, double d)
    -> decltype (VectorDouble::Constant (dim.size, d)) {
	return VectorDouble::Constant (dim.size, d);
}
inline auto linearAlgebraValueFilledWith (const Dimension<MatrixDouble> & dim, double d)
    -> decltype (MatrixDouble::Constant (dim.rows, dim.cols, d)) {
	return MatrixDouble::Constant (dim.rows, dim.cols, d);
}
///@}

///@{
/** Create a value filled with zeroes.
 * Versions for double, VectorDouble, MatrixDouble.
 * The right overload is selected by the Dimension<T> argument.
 * This argument also provide a size for Vector and Matrix variants.
 * This overloaded functions is used in template code (common interface for 3 types).
 */
inline double linearAlgebraZeroValue (const Dimension<double> &) {
	return 0.;
}
inline auto linearAlgebraZeroValue (const Dimension<VectorDouble> & dim)
    -> decltype (VectorDouble::Zero (dim.size)) {
	return VectorDouble::Zero (dim.size);
}
inline auto linearAlgebraZeroValue (const Dimension<MatrixDouble> & dim)
    -> decltype (MatrixDouble::Zero (dim.rows, dim.cols)) {
	return MatrixDouble::Zero (dim.rows, dim.cols);
}
///@}
} // namespace bpp

#endif // BPP_NEWPHYL_LINEARALGEBRAUTILS_H
