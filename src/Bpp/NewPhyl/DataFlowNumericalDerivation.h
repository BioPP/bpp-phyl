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


namespace bpp {
namespace DF {
  ///@{
	/** @brief Shift a value by a delta: <n>(delta, x) -> n * delta + x.
   *
	 * Defined for double, VectorDouble, MatrixDouble.
	 * For composite types, delta will be added to all elements.
	 *
	 * Should only be used to compute numerical derivatives of functions.
   * This is used to create the values of the different points of computation.
   * These values are then fed to duplicates of the function.
   *
	 * Delta is a special argument (not derivable, etc).
	 * Construction will merge chains of NumericalDerivationShiftDelta.
	 */
	template <typename T> class NumericalDerivationShiftDelta;

	// Factory functions
	template <> struct Builder<NumericalDerivationShiftDelta<double>> {
		static ValueRef<double> make (NodeRefVec && deps, int n,
		                              const Dimension<double> & targetDim = {});
	};
	template <> struct Builder<NumericalDerivationShiftDelta<VectorDouble>> {
		static ValueRef<VectorDouble> make (NodeRefVec && deps, int n,
		                                    const Dimension<VectorDouble> & targetDim);
	};
	template <> struct Builder<NumericalDerivationShiftDelta<MatrixDouble>> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, int n,
		                                    const Dimension<MatrixDouble> & targetDim);
	};
  ///@}

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

	template <> struct Builder<NumericalDerivationCombineShifted<double>> {
		static ValueRef<double> make (NodeRefVec && deps, const Vector<double> & coeffs,
		                              const Dimension<double> & targetDim = {});
	};
	template <> struct Builder<NumericalDerivationCombineShifted<VectorDouble>> {
		static ValueRef<VectorDouble> make (NodeRefVec && deps, const Vector<double> & coeffs,
		                                    const Dimension<VectorDouble> & targetDim);
	};
	template <> struct Builder<NumericalDerivationCombineShifted<MatrixDouble>> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, const Vector<double> & coeffs,
		                                    const Dimension<MatrixDouble> & targetDim);
	};
  ///@}

} // namespace DF
} // namespace bpp
#endif // BPP_NEWPHYL_DATAFLOWNUMERICALDERIVATION_H
