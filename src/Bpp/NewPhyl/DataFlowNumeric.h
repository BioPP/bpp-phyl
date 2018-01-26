//
// File: DataFlowNumeric.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-09-15 00:00:00
// Last modified: 2017-10-10
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

#ifndef BPP_NEWPHYL_DATAFLOWNUMERIC_H
#define BPP_NEWPHYL_DATAFLOWNUMERIC_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/LinearAlgebraFwd.h>

namespace bpp {
namespace DF {
	// Utils
	bool derivableIfAllDepsAre (const Node & toDerive, const Node & node);

	///@{
	/** @brief Constant representing a zero for type T.
	 *
	 * For vector or matrix, this means filled with zeroes.
	 * This node should be prefered to a Constant<T> initialized with zeroes when possible.
	 *
	 * The value is built lazily: starts invalid, will be created at first compute / access.
	 * This avoids memory allocations of matrices for temporary nodes.
	 * Most zero nodes will be discarded by optimizations in numeric nodes.
	 */
	template <typename T> class ConstantZero;

	// Builders
	template <> struct Builder<ConstantZero<double>> {
		static ValueRef<double> make (const Dimension<double> &);
		static ValueRef<double> make () { return make ({}); }
	};
	template <> struct Builder<ConstantZero<VectorDouble>> {
		static ValueRef<VectorDouble> make (const Dimension<VectorDouble> & dim);
	};
	template <> struct Builder<ConstantZero<MatrixDouble>> {
		static ValueRef<MatrixDouble> make (const Dimension<MatrixDouble> & dim);
	};
	///@}

	///@{
	/** @brief Constant representing a one for type T.
	 *
	 * For vector or matrix, this means filled with ones.
	 * Similar to ConstantZero.
	 */
	template <typename T> class ConstantOne;

	// Builders
	template <> struct Builder<ConstantOne<double>> {
		static ValueRef<double> make (const Dimension<double> &);
		static ValueRef<double> make () { return make ({}); }
	};
	template <> struct Builder<ConstantOne<VectorDouble>> {
		static ValueRef<VectorDouble> make (const Dimension<VectorDouble> & dim);
	};
	template <> struct Builder<ConstantOne<MatrixDouble>> {
		static ValueRef<MatrixDouble> make (const Dimension<MatrixDouble> & dim);
	};
	///@}

	// TODO
	template <typename Result, typename Dependencies> class CWiseAdd;

	template <> struct Builder<CWiseAdd<double, ReductionOfValue<double>>> {
		static ValueRef<double> make (NodeRefVec && deps, const Dimension<double> &);
	};
	template <> struct Builder<CWiseAdd<double, TupleOfValues<double, double>>> {
		static ValueRef<double> make (NodeRefVec && deps, const Dimension<double> &);
	};

	/* Double nodes.
	 */
	class AddDouble;
	template <> struct Builder<AddDouble> { static ValueRef<double> make (NodeRefVec && deps); };

	class MulDouble;
	template <> struct Builder<MulDouble> { static ValueRef<double> make (NodeRefVec && deps); };

	class NegDouble;
	template <> struct Builder<NegDouble> { static ValueRef<double> make (NodeRefVec && deps); };

	class ScalarProdDouble;
	template <> struct Builder<ScalarProdDouble> {
		static ValueRef<double> make (NodeRefVec && deps);
	};

	/* Vector nodes.
	 */
	class AddVectorDouble;
	template <> struct Builder<AddVectorDouble> {
		static ValueRef<VectorDouble> make (NodeRefVec && deps, const Dimension<VectorDouble> & dim);
	};

	class CWiseMulVectorDouble;
	template <> struct Builder<CWiseMulVectorDouble> {
		static ValueRef<VectorDouble> make (NodeRefVec && deps, const Dimension<VectorDouble> & dim);
	};

	class CWiseNegVectorDouble;
	template <> struct Builder<CWiseNegVectorDouble> {
		static ValueRef<VectorDouble> make (NodeRefVec && deps, const Dimension<VectorDouble> & dim);
	};

	class CWiseInverseVectorDouble;
	template <> struct Builder<CWiseInverseVectorDouble> {
		static ValueRef<VectorDouble> make (NodeRefVec && deps, const Dimension<VectorDouble> & dim);
	};

	class CWiseConstantPowVectorDouble;
	template <> struct Builder<CWiseConstantPowVectorDouble> {
		static ValueRef<VectorDouble> make (NodeRefVec && deps, const Dimension<VectorDouble> & dim,
		                                    double exp);
	};

	/* Matrix nodes.
	 */
	class AddMatrixDouble;
	template <> struct Builder<AddMatrixDouble> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, const Dimension<MatrixDouble> & dim);
	};

	class MulMatrixDouble;
	template <> struct Builder<MulMatrixDouble> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, const Dimension<MatrixDouble> & dim);
	};

	class CWiseMulMatrixDouble;
	template <> struct Builder<CWiseMulMatrixDouble> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, const Dimension<MatrixDouble> & dim);
	};

	/* Combinations
	 */
	class MulTransposedMatrixVectorDouble;
	template <> struct Builder<MulTransposedMatrixVectorDouble> {
		static ValueRef<VectorDouble> make (NodeRefVec && deps, const Dimension<VectorDouble> & dim);
	};

	class CWiseMulScalarVectorDouble;
	template <> struct Builder<CWiseMulScalarVectorDouble> {
		static ValueRef<VectorDouble> make (NodeRefVec && deps, const Dimension<VectorDouble> & dim);
	};

	class CWiseMulScalarMatrixDouble;
	template <> struct Builder<CWiseMulScalarMatrixDouble> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, const Dimension<MatrixDouble> & dim);
	};

} // namespace DF
} // namespace bpp
#endif // BPP_NEWPHYL_DATAFLOWNUMERIC_H
