//
// File: DataFlowNumeric.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2018-06-07
// Last modified: 2018-07-11
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

#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <string>
#include <tuple>
#include <type_traits>

#include "DataFlow.h"

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

	// Get dimensions of any matrix-like eigen object.
	template <typename Derived>
	MatrixDimension (const Eigen::MatrixBase<Derived> & m) : MatrixDimension (m.rows (), m.cols ()) {}

	bool operator== (const MatrixDimension & o) const { return rows == o.rows && cols == o.cols; }
	bool operator!= (const MatrixDimension & o) const { return !(*this == o); }
};

std::string to_string (const MatrixDimension & dim);

/// Eigen vector are matrices with 1 column.
inline MatrixDimension vectorDimension (Eigen::Index size) {
	return {size, 1};
}

/** Store a dimension for type T.
 * Declared but undefined by default.
 * Specialisations should be defined in the same header declaring the T type.
 * Specialisations should define a constructor from const T & : get the dimension of a T object.
 */
template <typename T> struct Dimension;

/** Specialisation of Dimension<T> for floating point types.
 * This is a dummy empty type, required by generic code below.
 */
template <> struct Dimension<double> {
	Dimension () = default;
	Dimension (const double &) {}
};
template <> struct Dimension<float> {
	Dimension () = default;
	Dimension (const float &) {}
};

template <typename T, typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
std::string to_string (const Dimension<T> &) {
	return "()";
}

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
 * Collection of overloaded numerical functions.
 * Not documented with doxygen, as this is only intended for use in dataflow nodes.
 */
namespace numeric {
	// Error util
	void checkDimensionIsSquare (const MatrixDimension & dim);

	// Create a zero value of the given dimension
	template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	T zero (const Dimension<T> &) {
		return T (0);
	}
	template <typename T, int Rows, int Cols>
	auto zero (const Dimension<Eigen::Matrix<T, Rows, Cols>> & dim)
	    -> decltype (Eigen::Matrix<T, Rows, Cols>::Zero (dim.rows, dim.cols)) {
		return Eigen::Matrix<T, Rows, Cols>::Zero (dim.rows, dim.cols);
	}

	// Create a one value of the given dimension
	template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	T one (const Dimension<T> &) {
		return T (1);
	}
	template <typename T, int Rows, int Cols>
	auto one (const Dimension<Eigen::Matrix<T, Rows, Cols>> & dim)
	    -> decltype (Eigen::Matrix<T, Rows, Cols>::Ones (dim.rows, dim.cols)) {
		return Eigen::Matrix<T, Rows, Cols>::Ones (dim.rows, dim.cols);
	}

	// Create an identity value of the given dimension (fails if not a square matrix)
	template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	T identity (const Dimension<T> &) {
		return T (1); // Equivalent to matrix of size 1x1
	}
	template <typename T, int Rows, int Cols>
	auto identity (const Dimension<Eigen::Matrix<T, Rows, Cols>> & dim)
	    -> decltype (Eigen::Matrix<T, Rows, Cols>::Identity (dim.rows, dim.cols)) {
		checkDimensionIsSquare (dim);
		return Eigen::Matrix<T, Rows, Cols>::Identity (dim.rows, dim.cols);
	}

	// Check if value is identity itself
	template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	bool isIdentity (const T & t) {
		return t == T (1);
	}
	template <typename Derived> bool isIdentity (const Eigen::MatrixBase<Derived> & m) {
		auto dim = Dimension<Derived> (m.derived ());
		return dim.rows == dim.cols && m == identity (dim);
	}

	/* Convert from F to R (with specific dimension).
	 * scalar -> scalar: simple cast.
	 * scalar -> matrix: fill the matrix with scalar value.
	 * matrix -> matrix: copy values, size must match (conversion between eigen dynamic/fixed types).
	 *
	 * Two APIs:
	 * r = convert(f, dim); -> convert returns a "converted" result
	 * convert (r, f, dim); -> convert does the assignment directly
	 */
	template <typename R, typename F,
	          typename = typename std::enable_if<std::is_arithmetic<R>::value>::type>
	R convert (const F & from, const Dimension<R> &) {
		return R (from); // scalar -> scalar
	}
	template <typename T, int Rows, int Cols, typename F,
	          typename = typename std::enable_if<std::is_arithmetic<F>::value>::type>
	auto convert (const F & from, const Dimension<Eigen::Matrix<T, Rows, Cols>> & dim)
	    -> decltype (Eigen::Matrix<T, Rows, Cols>::Constant (dim.rows, dim.cols, from)) {
		return Eigen::Matrix<T, Rows, Cols>::Constant (dim.rows, dim.cols, from); // scalar -> matrix
	}
	template <typename T, int Rows, int Cols, typename DerivedF>
	const DerivedF & convert (const Eigen::MatrixBase<DerivedF> & from,
	                          const Dimension<Eigen::Matrix<T, Rows, Cols>> & dim) {
		return from.derived (); // matrix -> matrix, conversion will be done in the assignment
	}
	template <typename R, typename F> void convert (R & r, const F & from, const Dimension<R> & dim) {
		r = convert (from, dim);
		assert (Dimension<R> (r) == dim); // debug post check of size
	}

	// Return a reference to the object for component-wise operations
	template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	T & cwise (T & t) {
		return t; // Do nothing for basic types
	}
	template <typename Derived>
	auto cwise (const Eigen::MatrixBase<Derived> & m) -> decltype (m.array ()) {
		return m.array (); // Use Array API in Eigen
	}
	template <typename Derived> auto cwise (Eigen::MatrixBase<Derived> & m) -> decltype (m.array ()) {
		return m.array (); // Use Array API in Eigen
	}

	// 1/x
	template <typename T, typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
	T inverse (T t) {
		return T (1) / t;
	}
	using Eigen::inverse;

	// x^y
	using Eigen::pow;
	using std::pow;

	// Numerical information as text
	template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	std::string debug (const T & t) {
		// For basic arithmetic scalar types, just print the value itself
		using std::to_string;
		return "value=" + to_string (t);
	}
	template <typename Derived> std::string debug (const Eigen::MatrixBase<Derived> & m) {
		// With matrices, check some numeric properties and encode results as text
		const auto dim = Dimension<Derived> (m.derived ());
		std::string props = "dim=" + to_string (dim) + " props={";
		const auto zero_value = zero (dim);
		// Properties on all elements
		if (m == zero_value)
			props += "[0]";
		if (m == one (dim))
			props += "[1]";
		if (isIdentity (m))
			props += "[I]";
		// Properties on any element
		if (m.array ().isNaN ().any ())
			props += "N";
		if (m.array ().isInf ().any ())
			props += "i";
		if ((m.array () == zero_value.array ()).any ())
			props += "0";
		if ((m.array () > zero_value.array ()).any ())
			props += "+";
		if ((m.array () < zero_value.array ()).any ())
			props += "-";
		props += "}";
		return props;
	}
} // namespace numeric

/******************************************************************************
 * Data flow nodes for those numerical functions.
 *
 * TODO what of rebuild ?
 * TODO numerical simplification: all deps constant => return constant ?
 * TODO add nodes from Numerical derivation
 */
namespace dataflow {
	// Error utils
	[[noreturn]] void failureDeltaNotDerivable (const std::type_info & contextNodeType);

	// Type tag to indicate a reduction operation (for +,*,...).
	template <typename T> struct ReductionOf;

	// Declaration of all defined nodes, in order of implementation.
	template <typename T> class ConstantZero;
	template <typename T> class ConstantOne;
	template <typename T> class NumericConstant;
	template <typename T> class NumericMutable;
	template <typename Result, typename From> class Convert;
	template <typename Result, typename From> class CWiseAdd;
	template <typename Result, typename From> class CWiseMul;
	template <typename T> class CWiseNegate;
	template <typename T> class CWiseInverse;
	template <typename T> class CWiseConstantPow;
	template <typename T0, typename T1> class ScalarProduct;
	template <typename R, typename T0, typename T1> class MatrixProduct;
	// TODO matrix multiply with transposed variants
	template <typename T> class ShiftDelta;
	template <typename T> class CombineDeltaShifted;

	// Utilities
	template <typename Predicate> void removeDependenciesIf (NodeRefVec & deps, Predicate p) {
		auto new_end = std::remove_if (deps.begin (), deps.end (), std::move (p));
		deps.erase (new_end, deps.end ()); // Truncate vector storage
	}
	inline bool allDerivable (const NodeRefVec & deps, const Node & node) {
		return std::all_of (deps.begin (), deps.end (),
		                    [&node](const NodeRef & dep) { return dep->isDerivable (node); });
	}

	/** r = 0 for each component.
	 * r: T.
	 * Value is only created at first use (lazy).
	 */
	template <typename T> class ConstantZero : public Value<T> {
	public:
		using Self = ConstantZero;

		static std::shared_ptr<Self> create (Context &, const Dimension<T> & dim = {}) {
			return std::make_shared<Self> (dim);
		}

		explicit ConstantZero (const Dimension<T> & dim)
		    : Value<T> (NodeRefVec{}), targetDimension (dim) {}

		std::string debugInfo () const override { return "targetDim=" + to_string (targetDimension); }

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

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<T>::create (c, targetDimension);
			}
			return this->shared_from_this (); // Return handle to self, as d(0)/dx = 0
		}
		bool isDerivable (const Node &) const final { return true; }

	private:
		void compute () final {
			using namespace numeric;
			this->accessValueMutable () = zero (targetDimension);
		}

		Dimension<T> targetDimension;
	};

	/** r = 1 for each component.
	 * r: T.
	 * Value is only created at first use (lazy).
	 */
	template <typename T> class ConstantOne : public Value<T> {
	public:
		using Self = ConstantOne;

		static std::shared_ptr<Self> create (Context &, const Dimension<T> & dim = {}) {
			return std::make_shared<Self> (dim);
		}

		explicit ConstantOne (const Dimension<T> & dim)
		    : Value<T> (NodeRefVec{}), targetDimension (dim) {}

		std::string debugInfo () const override { return "targetDim=" + to_string (targetDimension); }

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

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<T>::create (c, targetDimension);
			}
			return ConstantZero<T>::create (c, targetDimension);
		}
		bool isDerivable (const Node &) const final { return true; }

	private:
		void compute () final {
			using namespace numeric;
			this->accessValueMutable () = one (targetDimension);
		}

		Dimension<T> targetDimension;
	};

	/** r = const.
	 * r: T.
	 * Value is set at construction, and cannot change.
	 * Supports derivation.
	 */
	template <typename T> class NumericConstant : public Value<T> {
	public:
		using Self = NumericConstant;

		template <typename... Args> static std::shared_ptr<Self> create (Context &, Args &&... args) {
			return std::make_shared<Self> (std::forward<Args> (args)...);
		}

		/// Sets an initial value constructed with the given arguments.
		template <typename... Args>
		explicit NumericConstant (Args &&... args)
		    : Value<T> (NodeRefVec{}, std::forward<Args> (args)...) {
			this->makeValid (); // Always valid
		}

		std::string debugInfo () const override {
			using namespace numeric;
			return debug (this->accessValueConst ());
		}

		bool hasNumericalProperty (NumericalProperty prop) const final {
			using namespace numeric;
			const auto & value = this->accessValueConst ();
			switch (prop) {
			case NumericalProperty::Constant:
				return true;
			case NumericalProperty::Zero:
				return value == zero (Dimension<T> (value));
			case NumericalProperty::One:
				return value == one (Dimension<T> (value));
			case NumericalProperty::Identity:
				return isIdentity (value);
			default:
				return false;
			}
		}

		NodeRef derive (Context & c, const Node & node) final {
			const auto dim = Dimension<T> (this->accessValueConst ());
			if (&node == this) {
				return ConstantOne<T>::create (c, dim);
			}
			return ConstantZero<T>::create (c, dim);
		}
		bool isDerivable (const Node &) const final { return true; }

	private:
		void compute () final {
			// Constant is valid from construction
			failureComputeWasCalled (typeid (*this));
		}
	};

	/** r = v.
	 * r: T.
	 * Value is set at construction, and can be changed (will invalidate all dependent values).
	 * Supports derivation.
	 */
	template <typename T> class NumericMutable : public Value<T> {
	public:
		using Self = NumericMutable;

		template <typename... Args> static std::shared_ptr<Self> create (Context &, Args &&... args) {
			return std::make_shared<Self> (std::forward<Args> (args)...);
		}

		/// Sets an initial value constructed with the given arguments.
		template <typename... Args>
		explicit NumericMutable (Args &&... args)
		    : Value<T> (NodeRefVec{}, std::forward<Args> (args)...) {
			this->makeValid (); // Initial value is valid
		}

		/** @brief General case for modification of the T object.
		 *
		 * Takes a callable object (lamda, function pointer) that performs the modification.
		 * It must take a single T& as argument, which will refer to the T object to modify.
		 * The callable is called exactly once.
		 * TODO replace with view-struct that performs invalidate on destruction ?
		 */
		template <typename Callable> void modify (Callable && modifier) {
			this->invalidateRecursively ();
			std::forward<Callable> (modifier) (this->accessValueMutable ());
			this->makeValid ();
		}

		/// Setter with invalidation.
		void setValue (const T & t) {
			modify ([&t](T & v) { v = t; });
		}
		/// Setter with invalidation (movable value version).
		void setValue (T && t) {
			modify ([&t](T & v) { v = std::move (t); });
		}

		std::string debugInfo () const override {
			using namespace numeric;
			return debug (this->accessValueConst ());
		}

		bool hasNumericalProperty (NumericalProperty prop) const final {
			using namespace numeric;
			const auto & value = this->accessValueConst ();
			switch (prop) {
			case NumericalProperty::Zero:
				return value == zero (Dimension<T> (value));
			case NumericalProperty::One:
				return value == one (Dimension<T> (value));
			case NumericalProperty::Identity:
				return isIdentity (value);
			default:
				return false;
			}
		}

		NodeRef derive (Context & c, const Node & node) final {
			const auto dim = Dimension<T> (this->accessValueConst ());
			if (&node == this) {
				return ConstantOne<T>::create (c, dim);
			}
			return ConstantZero<T>::create (c, dim);
		}
		bool isDerivable (const Node &) const final { return true; }

	private:
		void compute () final {
			// Mutable is always valid
			failureComputeWasCalled (typeid (*this));
		}
	};

	/** r = f.
	 * r: R.
	 * f: F.
	 * Convert from F to R type, semantics of numeric::convert.
	 */
	template <typename R, typename F> class Convert : public Value<R> {
	public:
		using Self = Convert;

		static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim = {}) {
			// Check dependencies
			checkDependenciesNotNull (typeid (Self), deps);
			checkDependencyVectorSize (typeid (Self), deps, 1);
			checkNthDependencyIsValue<F> (typeid (Self), deps, 0);
			// Select node
			if (std::is_same<R, F>::value) {
				return convertRef<Value<R>> (deps[0]);
			} else if (deps[0]->hasNumericalProperty (NumericalProperty::Constant) &&
			           deps[0]->hasNumericalProperty (NumericalProperty::Zero)) {
				return ConstantZero<R>::create (c, dim);
			} else if (deps[0]->hasNumericalProperty (NumericalProperty::Constant) &&
			           deps[0]->hasNumericalProperty (NumericalProperty::One)) {
				return ConstantOne<R>::create (c, dim);
			} else {
				return std::make_shared<Self> (std::move (deps), dim);
			}
		}

		Convert (NodeRefVec && deps, const Dimension<R> & dim)
		    : Value<R> (std::move (deps)), targetDimension (dim) {}

		std::string debugInfo () const override {
			using namespace numeric;
			return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension);
		}

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<R>::create (c, targetDimension);
			}
			return Self::create (c, {this->dependency (0)->derive (c, node)}, targetDimension);
		}
		bool isDerivable (const Node & node) const final {
			return this->dependency (0)->isDerivable (node);
		}

	private:
		void compute () final {
			using namespace numeric;
			auto & result = this->accessValueMutable ();
			const auto & arg = accessValueConstCast<F> (*this->dependency (0));
			convert (result, arg, targetDimension);
		}

		Dimension<R> targetDimension;
	};

	/** r = x0 + x1 for each component.
	 * r: R.
	 * x0: T0.
	 * x1: T1.
	 *
	 * Values converted to R with the semantics of numeric::convert.
	 * Only defined for N = 2 for now.
	 * The generic version is horrible in C++11 (lack of auto return).
	 * Generic simplification routine is horrible too.
	 */
	template <typename R, typename T0, typename T1>
	class CWiseAdd<R, std::tuple<T0, T1>> : public Value<R> {
	public:
		using Self = CWiseAdd;

		static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim = {}) {
			// Check dependencies
			checkDependenciesNotNull (typeid (Self), deps);
			checkDependencyVectorSize (typeid (Self), deps, 2);
			checkNthDependencyIsValue<T0> (typeid (Self), deps, 0);
			checkNthDependencyIsValue<T1> (typeid (Self), deps, 1);
			// Select node implementation
			bool dep_0_zero = deps[0]->hasNumericalProperty (NumericalProperty::Constant) &&
			                  deps[0]->hasNumericalProperty (NumericalProperty::Zero);
			bool dep_1_zero = deps[1]->hasNumericalProperty (NumericalProperty::Constant) &&
			                  deps[1]->hasNumericalProperty (NumericalProperty::Zero);
			if (dep_0_zero && dep_1_zero) {
				return ConstantZero<R>::create (c, dim);
			} else if (dep_0_zero && !dep_1_zero) {
				return Convert<R, T1>::create (c, {deps[1]}, dim);
			} else if (!dep_0_zero && dep_1_zero) {
				return Convert<R, T0>::create (c, {deps[0]}, dim);
			} else {
				return std::make_shared<Self> (std::move (deps), dim);
			}
		}

		CWiseAdd (NodeRefVec && deps, const Dimension<R> & dim)
		    : Value<R> (std::move (deps)), targetDimension (dim) {}

		std::string debugInfo () const override {
			using namespace numeric;
			return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension);
		}

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<R>::create (c, targetDimension);
			}
			constexpr std::size_t n = 2;
			NodeRefVec derivedDeps (n);
			for (std::size_t i = 0; i < n; ++i) {
				derivedDeps[i] = this->dependency (i)->derive (c, node);
			}
			return Self::create (c, std::move (derivedDeps), targetDimension);
		}
		bool isDerivable (const Node & node) const final {
			return allDerivable (this->dependencies (), node);
		}

	private:
		void compute () final {
			using namespace numeric;
			auto & result = this->accessValueMutable ();
			const auto & x0 = accessValueConstCast<T0> (*this->dependency (0));
			const auto & x1 = accessValueConstCast<T1> (*this->dependency (1));
			cwise (result) = cwise (x0) + cwise (x1);
		}

		Dimension<R> targetDimension;
	};

	/** r = sum (x_i), for each component.
	 * r: R.
	 * x_i: T.
	 *
	 * Sum of any number of T values into R.
	 * Values converted to R with the semantics of numeric::convert.
	 */
	template <typename R, typename T> class CWiseAdd<R, ReductionOf<T>> : public Value<R> {
	public:
		using Self = CWiseAdd;

		static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim = {}) {
			// Check dependencies
			checkDependenciesNotNull (typeid (Self), deps);
			checkDependencyRangeIsValue<T> (typeid (Self), deps, 0, deps.size ());
			// Remove 0s from deps
			removeDependenciesIf (deps, [](const NodeRef & ref) {
				return ref->hasNumericalProperty (NumericalProperty::Constant) &&
				       ref->hasNumericalProperty (NumericalProperty::Zero);
			});
			// Select node implementation
			if (deps.size () == 0) {
				return ConstantZero<R>::create (c, dim);
			} else if (deps.size () == 1) {
				return Convert<R, T>::create (c, std::move (deps), dim);
			} else if (deps.size () == 2) {
				return CWiseAdd<R, std::tuple<T, T>>::create (c, std::move (deps), dim);
			} else {
				return std::make_shared<Self> (std::move (deps), dim);
			}
		}

		CWiseAdd (NodeRefVec && deps, const Dimension<R> & dim)
		    : Value<R> (std::move (deps)), targetDimension (dim) {}

		std::string debugInfo () const override {
			using namespace numeric;
			return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension);
		}

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<R>::create (c, targetDimension);
			}
			const auto n = this->nbDependencies ();
			NodeRefVec derivedDeps (n);
			for (std::size_t i = 0; i < n; ++i) {
				derivedDeps[i] = this->dependency (i)->derive (c, node);
			}
			return Self::create (c, std::move (derivedDeps), targetDimension);
		}
		bool isDerivable (const Node & node) const final {
			return allDerivable (this->dependencies (), node);
		}

	private:
		void compute () final {
			using namespace numeric;
			auto & result = this->accessValueMutable ();
			result = zero (targetDimension);
			for (const auto & depNodeRef : this->dependencies ()) {
				cwise (result) += cwise (accessValueConstCast<T> (*depNodeRef));
			}
		}

		Dimension<R> targetDimension;
	};

	/** r = x0 * x1 for each component.
	 * r: R.
	 * x0: T0.
	 * x1: T1.
	 *
	 * Values converted to R with the semantics of numeric::convert.
	 * Only defined for N = 2 for now (same constraints as CWiseAdd for genericity).
	 */
	template <typename R, typename T0, typename T1>
	class CWiseMul<R, std::tuple<T0, T1>> : public Value<R> {
	public:
		using Self = CWiseMul;

		static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim = {}) {
			// Check dependencies
			checkDependenciesNotNull (typeid (Self), deps);
			checkDependencyVectorSize (typeid (Self), deps, 2);
			checkNthDependencyIsValue<T0> (typeid (Self), deps, 0);
			checkNthDependencyIsValue<T1> (typeid (Self), deps, 1);
			// Return 0 if any 0.
			if (std::any_of (deps.begin (), deps.end (), [](const NodeRef & ref) {
				    return ref->hasNumericalProperty (NumericalProperty::Constant) &&
				           ref->hasNumericalProperty (NumericalProperty::Zero);
			    })) {
				return ConstantZero<R>::create (c, dim);
			}
			// Select node implementation
			bool dep_0_one = deps[0]->hasNumericalProperty (NumericalProperty::Constant) &&
			                 deps[0]->hasNumericalProperty (NumericalProperty::One);
			bool dep_1_one = deps[1]->hasNumericalProperty (NumericalProperty::Constant) &&
			                 deps[1]->hasNumericalProperty (NumericalProperty::One);
			if (dep_0_one && dep_1_one) {
				return ConstantOne<R>::create (c, dim);
			} else if (dep_0_one && !dep_1_one) {
				return Convert<R, T1>::create (c, {deps[1]}, dim);
			} else if (!dep_0_one && dep_1_one) {
				return Convert<R, T0>::create (c, {deps[0]}, dim);
			} else {
				return std::make_shared<Self> (std::move (deps), dim);
			}
		}

		CWiseMul (NodeRefVec && deps, const Dimension<R> & dim)
		    : Value<R> (std::move (deps)), targetDimension (dim) {}

		std::string debugInfo () const override {
			using namespace numeric;
			return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension);
		}

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<R>::create (c, targetDimension);
			}
			constexpr std::size_t n = 2;
			NodeRefVec addDeps (n);
			for (std::size_t i = 0; i < n; ++i) {
				NodeRefVec ithMulDeps = this->dependencies ();
				ithMulDeps[i] = this->dependency (i)->derive (c, node);
				addDeps[i] = Self::create (c, std::move (ithMulDeps), targetDimension);
			}
			return CWiseAdd<R, std::tuple<R, R>>::create (c, std::move (addDeps), targetDimension);
		}
		bool isDerivable (const Node & node) const final {
			return allDerivable (this->dependencies (), node);
		}

	private:
		void compute () final {
			using namespace numeric;
			auto & result = this->accessValueMutable ();
			const auto & x0 = accessValueConstCast<T0> (*this->dependency (0));
			const auto & x1 = accessValueConstCast<T1> (*this->dependency (1));
			cwise (result) = cwise (x0) * cwise (x1);
		}

		Dimension<R> targetDimension;
	};

	/** r = prod (x_i), for each component.
	 * r: R.
	 * x_i: T.
	 * Product of any number of T values into R.
	 * Values converted to R with the semantics of numeric::convert.
	 */
	template <typename R, typename T> class CWiseMul<R, ReductionOf<T>> : public Value<R> {
	public:
		using Self = CWiseMul;

		static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim = {}) {
			// Check dependencies
			checkDependenciesNotNull (typeid (Self), deps);
			checkDependencyRangeIsValue<T> (typeid (Self), deps, 0, deps.size ());
			// If there is a 0 return 0.
			if (std::any_of (deps.begin (), deps.end (), [](const NodeRef & ref) {
				    return ref->hasNumericalProperty (NumericalProperty::Constant) &&
				           ref->hasNumericalProperty (NumericalProperty::Zero);
			    })) {
				return ConstantZero<R>::create (c, dim);
			}
			// Remove 1s from deps
			removeDependenciesIf (deps, [](const NodeRef & ref) {
				return ref->hasNumericalProperty (NumericalProperty::Constant) &&
				       ref->hasNumericalProperty (NumericalProperty::One);
			});
			// Select node implementation
			if (deps.size () == 0) {
				return ConstantOne<R>::create (c, dim);
			} else if (deps.size () == 1) {
				return Convert<R, T>::create (c, std::move (deps), dim);
			} else if (deps.size () == 2) {
				return CWiseMul<R, std::tuple<T, T>>::create (c, std::move (deps), dim);
			} else {
				return std::make_shared<Self> (std::move (deps), dim);
			}
		}

		CWiseMul (NodeRefVec && deps, const Dimension<R> & dim)
		    : Value<R> (std::move (deps)), targetDimension (dim) {}

		std::string debugInfo () const override {
			using namespace numeric;
			return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension);
		}

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<R>::create (c, targetDimension);
			}
			const auto n = this->nbDependencies ();
			NodeRefVec addDeps (n);
			for (std::size_t i = 0; i < n; ++i) {
				NodeRefVec ithMulDeps = this->dependencies ();
				ithMulDeps[i] = this->dependency (i)->derive (c, node);
				addDeps[i] = Self::create (c, std::move (ithMulDeps), targetDimension);
			}
			return CWiseAdd<R, ReductionOf<R>>::create (c, std::move (addDeps), targetDimension);
		}
		bool isDerivable (const Node & node) const final {
			return allDerivable (this->dependencies (), node);
		}

	private:
		void compute () final {
			using namespace numeric;
			auto & result = this->accessValueMutable ();
			result = one (targetDimension);
			for (const auto & depNodeRef : this->dependencies ()) {
				cwise (result) *= cwise (accessValueConstCast<T> (*depNodeRef));
			}
		}

		Dimension<R> targetDimension;
	};

	/** r = -x, for each component.
	 * r, x: T.
	 */
	template <typename T> class CWiseNegate : public Value<T> {
	public:
		using Self = CWiseNegate;

		static ValueRef<T> create (Context & c, NodeRefVec && deps, const Dimension<T> & dim = {}) {
			// Check dependencies
			checkDependenciesNotNull (typeid (Self), deps);
			checkDependencyVectorSize (typeid (Self), deps, 1);
			checkNthDependencyIsValue<T> (typeid (Self), deps, 0);
			// Select node
			if (deps[0]->hasNumericalProperty (NumericalProperty::Constant) &&
			    deps[0]->hasNumericalProperty (NumericalProperty::Zero)) {
				return ConstantZero<T>::create (c, dim);
			} else {
				return std::make_shared<Self> (std::move (deps), dim);
			}
		}

		CWiseNegate (NodeRefVec && deps, const Dimension<T> & dim)
		    : Value<T> (std::move (deps)), targetDimension (dim) {}

		std::string debugInfo () const override {
			using namespace numeric;
			return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension);
		}

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<T>::create (c, targetDimension);
			}
			return Self::create (c, {this->dependency (0)->derive (c, node)}, targetDimension);
		}
		bool isDerivable (const Node & node) const final {
			return this->dependency (0)->isDerivable (node);
		}

	private:
		void compute () final {
			using namespace numeric;
			auto & result = this->accessValueMutable ();
			const auto & x = accessValueConstCast<T> (*this->dependency (0));
			cwise (result) = -cwise (x);
		}

		Dimension<T> targetDimension;
	};

	/** r = 1/x for each component.
	 * r, x: T.
	 */
	template <typename T> class CWiseInverse : public Value<T> {
	public:
		using Self = CWiseInverse;

		static ValueRef<T> create (Context & c, NodeRefVec && deps, const Dimension<T> & dim = {}) {
			// Check dependencies
			checkDependenciesNotNull (typeid (Self), deps);
			checkDependencyVectorSize (typeid (Self), deps, 1);
			checkNthDependencyIsValue<T> (typeid (Self), deps, 0);
			// Select node
			if (deps[0]->hasNumericalProperty (NumericalProperty::Constant) &&
			    deps[0]->hasNumericalProperty (NumericalProperty::One)) {
				return ConstantOne<T>::create (c, dim);
			} else {
				return std::make_shared<Self> (std::move (deps), dim);
			}
		}

		CWiseInverse (NodeRefVec && deps, const Dimension<T> & dim)
		    : Value<T> (std::move (deps)), targetDimension (dim) {}

		std::string debugInfo () const override {
			using namespace numeric;
			return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension);
		}

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<T>::create (c, targetDimension);
			}
			// -1/x^2 * x'
			const auto & dep = this->dependency (0);
			return CWiseMul<T, std::tuple<T, T>>::create (
			    c,
			    {CWiseConstantPow<T>::create (c, {dep}, -2., -1., targetDimension),
			     dep->derive (c, node)},
			    targetDimension);
		}
		bool isDerivable (const Node & node) const final {
			return this->dependency (0)->isDerivable (node);
		}

	private:
		void compute () final {
			using namespace numeric;
			auto & result = this->accessValueMutable ();
			const auto & x = accessValueConstCast<T> (*this->dependency (0));
			cwise (result) = inverse (cwise (x));
		}

		Dimension<T> targetDimension;
	};

	/** r = factor * pow (x, exponent) for each component.
	 * r, x: T.
	 * exponent, factor: double (constant parameter of the node).
	 */
	template <typename T> class CWiseConstantPow : public Value<T> {
	public:
		using Self = CWiseConstantPow;

		static ValueRef<T> create (Context & c, NodeRefVec && deps, double exponent, double factor,
		                           const Dimension<T> & dim = {}) {
			// Check dependencies
			checkDependenciesNotNull (typeid (Self), deps);
			checkDependencyVectorSize (typeid (Self), deps, 1);
			checkNthDependencyIsValue<T> (typeid (Self), deps, 0);
			// Select node implementation
			if (exponent == 0. || (deps[0]->hasNumericalProperty (NumericalProperty::Constant) &&
			                       deps[0]->hasNumericalProperty (NumericalProperty::One))) {
				// pow (x, exponent) == 1
				using namespace numeric;
				return NumericConstant<T>::create (c, factor * one (dim));
			} else if (exponent == 1.) {
				// pow (x, exponent) == x
				return CWiseMul<T, std::tuple<double, T>>::create (
				    c, {NumericConstant<double>::create (c, factor), deps[0]}, dim);
			} else if (exponent == -1.) {
				// pow (x, exponent) = 1/x
				return CWiseMul<T, std::tuple<double, T>>::create (
				    c,
				    {NumericConstant<double>::create (c, factor),
				     CWiseInverse<T>::create (c, std::move (deps), dim)},
				    dim);
			} else {
				return std::make_shared<Self> (std::move (deps), exponent, factor, dim);
			}
		}

		CWiseConstantPow (NodeRefVec && deps, double exponent, double factor, const Dimension<T> & dim)
		    : Value<T> (std::move (deps)),
		      targetDimension (dim),
		      exponent (exponent),
		      factor (factor) {}

		std::string debugInfo () const override {
			using namespace numeric;
			return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension) +
			       " exponent=" + std::to_string (exponent) + " factor=" + std::to_string (factor);
		}

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<T>::create (c, targetDimension);
			}
			// factor * (exponent * x^(exponent - 1)) * x'
			const auto & dep = this->dependency (0);
			auto dpow = Self::create (c, {dep}, exponent - 1., factor * exponent, targetDimension);
			return CWiseMul<T, std::tuple<T, T>>::create (c, {dpow, dep->derive (c, node)},
			                                              targetDimension);
		}
		bool isDerivable (const Node & node) const final {
			return this->dependency (0)->isDerivable (node);
		}

	private:
		void compute () final {
			using namespace numeric;
			auto & result = this->accessValueMutable ();
			const auto & x = accessValueConstCast<T> (*this->dependency (0));
			cwise (result) = factor * pow (cwise (x), exponent);
		}

		Dimension<T> targetDimension;
		double exponent;
		double factor;
	};

	/** r = x0 * x1 (dot product)
	 * r: double.
	 * x0: T0 (vector-like).
	 * x1: T1 (vector-like).
	 */
	template <typename T0, typename T1> class ScalarProduct : public Value<double> {
	public:
		using Self = ScalarProduct;

		static ValueRef<double> create (Context & c, NodeRefVec && deps) {
			// Check dependencies
			checkDependenciesNotNull (typeid (Self), deps);
			checkDependencyVectorSize (typeid (Self), deps, 2);
			checkNthDependencyIsValue<T0> (typeid (Self), deps, 0);
			checkNthDependencyIsValue<T1> (typeid (Self), deps, 1);
			// Select node
			bool dep_0_zero = deps[0]->hasNumericalProperty (NumericalProperty::Constant) &&
			                  deps[0]->hasNumericalProperty (NumericalProperty::Zero);
			bool dep_1_zero = deps[1]->hasNumericalProperty (NumericalProperty::Constant) &&
			                  deps[1]->hasNumericalProperty (NumericalProperty::Zero);
			if (dep_0_zero || dep_1_zero) {
				return ConstantZero<double>::create (c);
			} else {
				return std::make_shared<Self> (std::move (deps));
			}
		}

		ScalarProduct (NodeRefVec && deps) : Value<double> (std::move (deps)) {}

		std::string debugInfo () const override {
			using namespace numeric;
			return debug (this->accessValueConst ());
		}

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<double>::create (c);
			}
			const auto & x0 = this->dependency (0);
			const auto & x1 = this->dependency (1);
			auto dx0_prod = Self::create (c, {x0->derive (c, node), x1});
			auto dx1_prod = Self::create (c, {x0, x1->derive (c, node)});
			return CWiseAdd<double, std::tuple<double, double>>::create (c, {dx0_prod, dx1_prod});
		}
		bool isDerivable (const Node & node) const final {
			return allDerivable (this->dependencies (), node);
		}

	private:
		void compute () final {
			auto & result = this->accessValueMutable ();
			const auto & x0 = accessValueConstCast<T0> (*this->dependency (0));
			const auto & x1 = accessValueConstCast<T1> (*this->dependency (1));
			result = x0.dot (x1); // Using lhs.dot(rhs) method from Eigen only
		}
	};

	/** r = x0 * x1 (matrix product).
	 * r: R (matrix).
	 * x0: T0 (matrix).
	 * x1: T1 (matrix).
	 */
	template <typename R, typename T0, typename T1> class MatrixProduct : public Value<R> {
	public:
		using Self = MatrixProduct;

		static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim) {
			// Check dependencies
			checkDependenciesNotNull (typeid (Self), deps);
			checkDependencyVectorSize (typeid (Self), deps, 2);
			checkNthDependencyIsValue<T0> (typeid (Self), deps, 0);
			checkNthDependencyIsValue<T1> (typeid (Self), deps, 1);
			// Return 0 if any 0.
			if (std::any_of (deps.begin (), deps.end (), [](const NodeRef & ref) {
				    return ref->hasNumericalProperty (NumericalProperty::Constant) &&
				           ref->hasNumericalProperty (NumericalProperty::Zero);
			    })) {
				return ConstantZero<R>::create (c, dim);
			}
			// Select node implementation
			bool dep_0_identity = deps[0]->hasNumericalProperty (NumericalProperty::Constant) &&
			                      deps[0]->hasNumericalProperty (NumericalProperty::Identity);
			bool dep_1_identity = deps[1]->hasNumericalProperty (NumericalProperty::Constant) &&
			                      deps[1]->hasNumericalProperty (NumericalProperty::Identity);
			if (dep_0_identity && dep_1_identity) {
				// No specific class for Identity
				using namespace numeric;
				return NumericConstant<R>::create (c, identity (dim));
			} else if (dep_0_identity && !dep_1_identity) {
				return Convert<R, T1>::create (c, {deps[1]}, dim);
			} else if (!dep_0_identity && dep_1_identity) {
				return Convert<R, T0>::create (c, {deps[0]}, dim);
			} else {
				return std::make_shared<Self> (std::move (deps), dim);
			}
		}

		MatrixProduct (NodeRefVec && deps, const Dimension<R> & dim)
		    : Value<R> (std::move (deps)), targetDimension (dim) {}

		std::string debugInfo () const override {
			using namespace numeric;
			return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension);
		}

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<R>::create (c, targetDimension);
			}
			const auto & x0 = this->dependency (0);
			const auto & x1 = this->dependency (1);
			auto dx0_prod = Self::create (c, {x0->derive (c, node), x1}, targetDimension);
			auto dx1_prod = Self::create (c, {x0, x1->derive (c, node)}, targetDimension);
			return CWiseAdd<R, std::tuple<R, R>>::create (c, {dx0_prod, dx1_prod}, targetDimension);
		}
		bool isDerivable (const Node & node) const final {
			return allDerivable (this->dependencies (), node);
		}

	private:
		void compute () final {
			auto & result = this->accessValueMutable ();
			const auto & x0 = accessValueConstCast<T0> (*this->dependency (0));
			const auto & x1 = accessValueConstCast<T1> (*this->dependency (1));
			result.noalias () = x0 * x1;
		}

		Dimension<R> targetDimension;
	};

	// TODO transposed variants ?  r.noalias () = lhs.transpose () * rhs;

	/** r = n * delta + x.
	 * r: T.
	 * delta: double.
	 * x: T.
	 * n: constant int.
	 * Order of dependencies: (delta, x).
	 *
	 * Adds n * delta to all values (component wise) of x.
	 * Used to generate x +/- delta values for numerical derivation.
	 */
	template <typename T> class ShiftDelta : public Value<T> {
	public:
		using Self = ShiftDelta;

		static ValueRef<T> create (Context & c, NodeRefVec && deps, int n,
		                           const Dimension<T> & dim = {}) {
			// Check dependencies
			checkDependenciesNotNull (typeid (Self), deps);
			checkDependencyVectorSize (typeid (Self), deps, 2);
			checkNthDependencyIsValue<double> (typeid (Self), deps, 0);
			checkNthDependencyIsValue<T> (typeid (Self), deps, 1);
			// Detect if we have a chain of ShiftDelta with the same delta.
			auto & delta = deps[0];
			auto & x = deps[1];
			auto * xAsShiftDelta = dynamic_cast<const ShiftDelta<T> *> (x.get ());
			if (xAsShiftDelta != nullptr && xAsShiftDelta->dependency (0) == delta) {
				// Merge with ShiftDelta dependency by summing the n.
				return Self::create (c, NodeRefVec{x->dependencies ()}, n + xAsShiftDelta->getN (), dim);
			}
			// Not a merge, select node implementation.
			if (n == 0 || (delta->hasNumericalProperty (NumericalProperty::Constant) &&
			               delta->hasNumericalProperty (NumericalProperty::Zero))) {
				return convertRef<Value<T>> (x);
			} else {
				return std::make_shared<Self> (std::move (deps), n, dim);
			}
		}

		ShiftDelta (NodeRefVec && deps, int n, const Dimension<T> & dim)
		    : Value<T> (std::move (deps)), targetDimension (dim), n (n) {}

		std::string debugInfo () const override {
			using namespace numeric;
			return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension) +
			       " n=" + std::to_string (n);
		}

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<T>::create (c, targetDimension);
			}
			auto & delta = this->dependency (0);
			auto & x = this->dependency (1);
			return Self::create (c, {delta->derive (c, node), x->derive (c, node)}, n, targetDimension);
		}
		bool isDerivable (const Node & node) const final {
			return allDerivable (this->dependencies (), node);
		}

		int getN () const { return n; }

	private:
		void compute () final {
			using namespace numeric;
			auto & result = this->accessValueMutable ();
			const auto & delta = accessValueConstCast<double> (*this->dependency (0));
			const auto & x = accessValueConstCast<T> (*this->dependency (1));
			cwise (result) = n * delta + cwise (x);
		}

		Dimension<T> targetDimension;
		int n;
	};

	/** r = (1/delta)^n * sum_i coeffs_i * x_i.
	 * r: T.
	 * delta: double.
	 * x_i: T.
	 * n: constant int.
	 * coeffs_i: constant double.
	 * Order of dependencies: (delta, x_i).
	 *
	 * Weighted sum of dependencies, multiplied by a double.
	 * Used to combine f(x+n*delta) values in numerical derivation.
	 * Lambda represents the 1/delta^n for a nth-order numerical derivation.
	 */
	template <typename T> class CombineDeltaShifted : public Value<T> {
	public:
		using Self = CombineDeltaShifted;

		static ValueRef<T> create (Context & c, NodeRefVec && deps, int n,
		                           std::vector<double> && coeffs, const Dimension<T> & dim = {}) {
			// Check dependencies
			checkDependenciesNotNull (typeid (Self), deps);
			checkDependencyVectorSize (typeid (Self), deps, 1 + coeffs.size ());
			checkNthDependencyIsValue<double> (typeid (Self), deps, 0);
			checkDependencyRangeIsValue<T> (typeid (Self), deps, 1, deps.size ());
			// TODO merge with dependencies if CombineDeltaShifted with same delta.
			// TODO remove constant 0, or stuff with 0 factor.
			return std::make_shared<Self> (std::move (deps), n, std::move (coeffs), dim);
		}

		CombineDeltaShifted (NodeRefVec && deps, int n, std::vector<double> && coeffs,
		                     const Dimension<T> & dim)
		    : Value<T> (std::move (deps)), targetDimension (dim), coeffs (std::move (coeffs)), n (n) {
			assert (this->coeffs.size () + 1 == this->nbDependencies ());
		}

		std::string debugInfo () const override {
			using namespace numeric;
			std::string s = debug (this->accessValueConst ()) +
			                " targetDim=" + to_string (targetDimension) + " n=" + std::to_string (n) +
			                " coeffs={";
			if (!coeffs.empty ()) {
				s += std::to_string (coeffs[0]);
				for (std::size_t i = 1; i < coeffs.size (); ++i) {
					s += ';' + std::to_string (coeffs[i]);
				}
			}
			s += '}';
			return s;
		}

		NodeRef derive (Context & c, const Node & node) final {
			if (&node == this) {
				return ConstantOne<T>::create (c, targetDimension);
			}
			// For simplicity, we assume delta is a constant with respect to derivation node.
			auto & delta = this->dependency (0);
			if (isTransitivelyDependentOn (node, *delta)) {
				// Fail if delta is not constant for node.
				failureDeltaNotDerivable (typeid (Self));
			}
			// Derivation is a simple weighted sums of derivatives.
			const auto nbDeps = this->nbDependencies ();
			NodeRefVec derivedDeps (nbDeps);
			derivedDeps[0] = delta;
			for (std::size_t i = 1; i < nbDeps; ++i) {
				derivedDeps[i] = this->dependency (i)->derive (c, node);
			}
			return Self::create (c, std::move (derivedDeps), n, std::vector<double>{coeffs},
			                     targetDimension);
		}
		bool isDerivable (const Node & node) const final {
			return allDerivable (this->dependencies (), node); // FIXME delta arg ?
		}

		const std::vector<double> & getCoeffs () const { return coeffs; }
		int getN () const { return n; }

	private:
		void compute () final {
			using namespace numeric;
			auto & result = this->accessValueMutable ();
			const auto & delta = accessValueConstCast<double> (*this->dependency (0));
			const double lambda = pow (delta, -n);
			result = zero (targetDimension);
			for (std::size_t i = 0; i < coeffs.size (); ++i) {
				const auto & x = accessValueConstCast<T> (*this->dependency (1 + i));
				cwise (result) += (lambda * coeffs[i]) * cwise (x);
			}
		}

		Dimension<T> targetDimension;
		std::vector<double> coeffs;
		int n;
	};

	// Precompiled instantiations
	extern template class ConstantZero<double>;
	extern template class ConstantZero<Eigen::VectorXd>;
	extern template class ConstantZero<Eigen::MatrixXd>;

	extern template class ConstantOne<double>;
	extern template class ConstantOne<Eigen::VectorXd>;
	extern template class ConstantOne<Eigen::MatrixXd>;

	extern template class NumericConstant<double>;
	extern template class NumericConstant<Eigen::VectorXd>;
	extern template class NumericConstant<Eigen::MatrixXd>;

	extern template class NumericMutable<double>;
	extern template class NumericMutable<Eigen::VectorXd>;
	extern template class NumericMutable<Eigen::MatrixXd>;

	extern template class Convert<double, double>;
	extern template class Convert<Eigen::VectorXd, Eigen::VectorXd>;
	extern template class Convert<Eigen::MatrixXd, Eigen::MatrixXd>;
	extern template class Convert<Eigen::VectorXd, double>;
	extern template class Convert<Eigen::MatrixXd, double>;

	extern template class CWiseAdd<double, std::tuple<double, double>>;
	extern template class CWiseAdd<Eigen::VectorXd, std::tuple<Eigen::VectorXd, Eigen::VectorXd>>;
	extern template class CWiseAdd<Eigen::MatrixXd, std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>>;

	extern template class CWiseAdd<double, ReductionOf<double>>;
	extern template class CWiseAdd<Eigen::VectorXd, ReductionOf<Eigen::VectorXd>>;
	extern template class CWiseAdd<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

	extern template class CWiseMul<double, std::tuple<double, double>>;
	extern template class CWiseMul<Eigen::VectorXd, std::tuple<Eigen::VectorXd, Eigen::VectorXd>>;
	extern template class CWiseMul<Eigen::MatrixXd, std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>>;
	extern template class CWiseMul<Eigen::VectorXd, std::tuple<double, Eigen::VectorXd>>;
	extern template class CWiseMul<Eigen::MatrixXd, std::tuple<double, Eigen::MatrixXd>>;

	extern template class CWiseMul<double, ReductionOf<double>>;
	extern template class CWiseMul<Eigen::VectorXd, ReductionOf<Eigen::VectorXd>>;
	extern template class CWiseMul<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

	extern template class CWiseNegate<double>;
	extern template class CWiseNegate<Eigen::VectorXd>;
	extern template class CWiseNegate<Eigen::MatrixXd>;

	extern template class CWiseInverse<double>;
	extern template class CWiseInverse<Eigen::VectorXd>;
	extern template class CWiseInverse<Eigen::MatrixXd>;

	extern template class CWiseConstantPow<double>;
	extern template class CWiseConstantPow<Eigen::VectorXd>;
	extern template class CWiseConstantPow<Eigen::MatrixXd>;

	extern template class ScalarProduct<Eigen::VectorXd, Eigen::VectorXd>;

	extern template class MatrixProduct<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>;

	extern template class ShiftDelta<double>;
	extern template class ShiftDelta<Eigen::VectorXd>;
	extern template class ShiftDelta<Eigen::MatrixXd>;

	extern template class CombineDeltaShifted<double>;
	extern template class CombineDeltaShifted<Eigen::VectorXd>;
	extern template class CombineDeltaShifted<Eigen::MatrixXd>;
} // namespace dataflow
} // namespace bpp

#endif // BPP_NEWPHYL_DATAFLOWNUMERIC_H
