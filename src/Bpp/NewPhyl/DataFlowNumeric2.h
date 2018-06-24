//
// File: DataFlowNumeric2.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2018-06-07
// Last modified: 2018-06-07
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

#ifndef BPP_NEWPHYL_DATAFLOWNUMERIC2_H
#define BPP_NEWPHYL_DATAFLOWNUMERIC2_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Eigen/Core>
#include <cassert>
#include <string>

namespace bpp {
/******************************************************************************
 * Dimension.
 */

/// Basic matrix dimension type
struct MatrixDimension {
	Eigen::Index rows{};
	Eigen::Index cols{};

	MatrixDimension () = default;
	MatrixDimension (Eigen::Index rows, Eigen::Index cols) : rows (rows), cols (cols) {}
};

/// Eigen vector are matrices with 1 column.
inline MatrixDimension vectorDimension (Eigen::Index size) {
	return {size, 1};
}

/** Store a dimension for type T.
 * Declared but undefined by default.
 * Specialisations should be defined in the same header declaring the T type.
 */
template <typename T> struct Dimension;

/** Specialisation of Dimension<T> for double.
 * This is a dummy empty type, required by generic code below.
 */
template <> struct Dimension<double> {};

/** Specialisation of Dimension<T> for eigen matrix types.
 * Note that in Eigen, a vector is a matrix with one column.
 * Redirect to MatrixDimension for all eigen matrix variants.
 */
template <typename T, int Rows, int Cols>
struct Dimension<Eigen::Matrix<T, Rows, Cols>> : MatrixDimension {
	using MatrixDimension::MatrixDimension; // Have the same constructors as MatrixDimension
	Dimension (const MatrixDimension & dim) : MatrixDimension (dim) {} // From MatrixDimension
};

/******************************************************************************
 * Numerical functions.
 */

// Create a zero value of the given dimension
inline double numericalZeroValue (const Dimension<double> &) {
	return 0.;
}
template <typename T, int Rows, int Cols>
auto numericalZeroValue (const Dimension<Eigen::Matrix<T, Rows, Cols>> & dim)
    -> decltype (Eigen::Matrix<T, Rows, Cols>::Zero (dim.rows, dim.cols)) {
	return Eigen::Matrix<T, Rows, Cols>::Zero (dim.rows, dim.cols);
}

// Create a one value of the given dimension
inline double numericalOneValue (const Dimension<double> &) {
	return 1.;
}
template <typename T, int Rows, int Cols>
auto numericalOneValue (const Dimension<Eigen::Matrix<T, Rows, Cols>> & dim)
    -> decltype (Eigen::Matrix<T, Rows, Cols>::Ones (dim.rows, dim.cols)) {
	return Eigen::Matrix<T, Rows, Cols>::Ones (dim.rows, dim.cols);
}

/******************************************************************************
 * Data flow nodes for those numerical functions.
 * TODO debug !
 * TODO what of rebuild ?
 */
namespace dataflow {
	/** Lazily created numeric constant equal to zero for the type.
	 */
	template <typename T> class ConstantZero : public Value<T> {
	public:
		explicit ConstantZero (const Dimension<T> & dim = {})
		    : Value<T> (noDependency), targetDimension (dim) {}

		bool hasNumericalProperty (NumericalProperty prop) const final {
			switch (prop) {
			case NumericalProperty::Constant:
				return true;
			case NumericalProperty::Zero:
				return true;
			default:
				return false;
			}
		}

		NodeRef derive (const Node &) final {
			return this->shared_from_this (); // Return handle to self, as derive (0) = 0
		}
		bool isDerivable (const Node &) const final { return true; }

	private:
		void compute () final { this->accessValueMutable () = numericalZeroValue (targetDimension); }

		Dimension<T> targetDimension;
	};

	// TODO merging builder::make ?

	/** Lazily created numeric constant equal to one for the type.
	 */
	template <typename T> class ConstantOne : public Value<T> {
	public:
		explicit ConstantOne (const Dimension<T> & dim = {})
		    : Value<T> (noDependency), targetDimension (dim) {}

		bool hasNumericalProperty (NumericalProperty prop) const final {
			switch (prop) {
			case NumericalProperty::Constant:
				return true;
			case NumericalProperty::One:
				return true;
			default:
				return false;
			}
		}

		NodeRef derive (const Node &) final { return makeNode<ConstantZero<T>> (targetDimension); }
		bool isDerivable (const Node &) const final { return true; }

	private:
		void compute () final { this->accessValueMutable () = numericalOneValue (targetDimension); }

		Dimension<T> targetDimension;
	};

	// TODO Add

	template <typename T> struct ReductionOf; // Type tag

	template <typename Result, typename From> class CWiseAdd;

	/** Addition of any number of T into R.
	 */
	template <typename R, typename T> class CWiseAdd<R, ReductionOf<T>> : public Value<R> {
	public:
		CWiseAdd (NodeRefVec && deps, const Dimension<R> & dim)
		    : Value<R> (std::move (deps)), targetDimension (dim) {}

	private:
		void compute () final {
			auto & result = this->accessValueMutable ();
			result = numericalZeroValue (targetDimension);
			for (const auto & depNodeRef : this->dependencies ()) {
				result += accessValueConstCast<T> (*depNodeRef);
			}
		}

		Dimension<R> targetDimension;
	};

	template <typename R, typename T> struct Builder<CWiseAdd<R, ReductionOf<T>>> {
		static ValueRef<R> make (NodeRefVec && deps, const Dimension<R> & dim = {}) {
			using NodeType = CWiseAdd<R, ReductionOf<T>>;
			// Check dependencies
			checkDependenciesNotNull (typeid (NodeType), deps);
			checkDependencyRangeIsValue<T> (typeid (NodeType), deps, 0, deps.size ());
			return std::make_shared<CWiseAdd<R, ReductionOf<T>>> (std::move (deps), dim);
		}
	};
} // namespace dataflow
} // namespace bpp

#endif // BPP_NEWPHYL_DATAFLOWNUMERIC2_H
