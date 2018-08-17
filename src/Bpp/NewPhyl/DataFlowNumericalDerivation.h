//
// File: DataFlowNumericalDerivation.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-12-19 00:00:00
// Last modified: 2018-01-31
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

#ifndef BPP_NEWPHYL_DATAFLOWNUMERICALDERIVATION_H
#define BPP_NEWPHYL_DATAFLOWNUMERICALDERIVATION_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/LinearAlgebraFwd.h>
#include <vector>

namespace bpp {
namespace DF {
	///@{
	/** @brief Make a linear combination of values (constant factors), multiplied by a factor.
	 *
	 * Defined for double, VectorDouble, MatrixDouble.
	 * For composite types, delta will be added to all elements.
	 * It performs: <coeffs>(lambda, deps) -> lambda * sum_i (deps[i] * coeffs[i]).
	 *
	 * Should only be used to compute numerical derivatives of functions.
	 * This is used to combine the values from the function computed in different points.
	 * Lambda is intended to be 1/delta^n, where delta is the shift value.
	 *
	 * Lambda is a special argument (not derivable, etc).
	 * Construction will merge layers of consecutive NumericalDerivationCombineShifted.
	 */
	template <typename T> class NumericalDerivationCombineShifted;
	///@}

	enum class NumericalDerivativeType { Disabled, ThreePoints, FivePoints };

	struct NumericalDerivativeConfiguration {
		NumericalDerivativeType type{NumericalDerivativeType::Disabled};
		ValueRef<double> delta{};
	};

	template <typename T>
	ValueRef<T> generateDerivative (const NumericalDerivativeConfiguration & config,
	                                const ValueRef<T> & derivedNode, const Node & variable) {
		switch (config.type) {
		case NumericalDerivativeType::ThreePoints: {
			// Shift {-1, +1}, coeffs {-0.5, +0.5}
		} break;
		case NumericalDerivativeType::FivePoints: {
			// SHift {-2, -1, +1, +2}, coeffs {1/12, -2/3, 2/3, -1/12}
		} break;
		default:
			throw "error"; // FIXME
		}
	}

} // namespace DF
} // namespace bpp
#endif // BPP_NEWPHYL_DATAFLOWNUMERICALDERIVATION_H
