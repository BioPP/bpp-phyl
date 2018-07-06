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
#include <algorithm>
#include <cassert>
#include <string>
#include <tuple>
#include <type_traits>

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
};

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
 */
namespace numeric {
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

	// Convert from F to R (with specific dimension)
	template <typename R, typename F,
	          typename = typename std::enable_if<std::is_convertible<F, R>::value>::type>
	R convert (const F & from, const Dimension<R> &) {
		return R (from); // For simple arithmetic types
	}
	template <typename T, int Rows, int Cols, typename F>
	auto convert (const F & from, const Dimension<Eigen::Matrix<T, Rows, Cols>> & dim)
	    -> decltype (Eigen::Matrix<T, Rows, Cols>::Constant (T (from), dim.rows, dim.cols)) {
		// Eigen case, build a matrix/vector filled with constant value.
		return Eigen::Matrix<T, Rows, Cols>::Constant (dim.rows, dim.cols, T (from));
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
} // namespace numeric

/******************************************************************************
 * Data flow nodes for those numerical functions.
 * TODO debug !
 * TODO what of rebuild ?
 * TODO add other nodes from DFN1
 * TODO add nodes from Numerical derivation
 * TODO rename to DFN1
 * TODO rm CPP14, DFN1, integer range uses
 */
namespace dataflow {
	template <typename T> struct ReductionOf; // Type tag

	// Declaration of all defined nodes, in order of implementation.
	template <typename T> class ConstantZero;
	template <typename T> class ConstantOne;
	template <typename T> class NumericConstant;
	template <typename T> class NumericMutable;
	template <typename Result, typename From> class Convert;
	template <typename Result, typename From> class CWiseAdd;
	template <typename Result, typename From> class CWiseMul;

	// Utilities
	template <typename Predicate> void removeDependenciesIf (NodeRefVec & deps, Predicate p) {
		auto new_end = std::remove_if (deps.begin (), deps.end (), std::move (p));
		deps.erase (new_end, deps.end ()); // Truncate vector storage
	}
	inline bool allDerivable (const NodeRefVec & deps, const Node & node) {
		return std::all_of (deps.begin (), deps.end (),
		                    [&node](const NodeRef & dep) { return dep->isDerivable (node); });
	}

	/** Lazily created numeric constant equal to zero for the type.
	 */
	template <typename T> class ConstantZero : public Value<T> {
	public:
		using Self = ConstantZero;

		static std::shared_ptr<Self> create (Context &, const Dimension<T> & dim = {}) {
			return std::make_shared<Self> (dim);
		}

		explicit ConstantZero (const Dimension<T> & dim)
		    : Value<T> (NodeRefVec{}), targetDimension (dim) {}

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

		NodeRef derive (Context &, const Node &) final {
			return this->shared_from_this (); // Return handle to self, as derive (0) = 0
		}
		bool isDerivable (const Node &) const final { return true; }

	private:
		void compute () final {
			using namespace numeric;
			this->accessValueMutable () = zero (targetDimension);
		}

		Dimension<T> targetDimension;
	};

	/** Lazily created numeric constant equal to one for the type.
	 */
	template <typename T> class ConstantOne : public Value<T> {
	public:
		using Self = ConstantOne;

		static std::shared_ptr<Self> create (Context &, const Dimension<T> & dim = {}) {
			return std::make_shared<Self> (dim);
		}

		explicit ConstantOne (const Dimension<T> & dim)
		    : Value<T> (NodeRefVec{}), targetDimension (dim) {}

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

		NodeRef derive (Context & c, const Node &) final {
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

	/** Numeric constant of type T.
	 * Supports derivation.
	 * Eagerly constructed (value is created at construction).
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
			default:
				return false;
			}
		}

		NodeRef derive (Context & c, const Node &) final {
			return ConstantZero<T>::create (c, Dimension<T> (this->accessValueConst ()));
		}
		bool isDerivable (const Node &) const final { return true; }

	private:
		void compute () final {
			// Constant is valid from construction
			failureComputeWasCalled (typeid (*this));
		}
	};

	/** Numeric mutable value of type T.
	 * Supports derivation.
	 * Eagerly constructed (value is created at construction).
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

		bool hasNumericalProperty (NumericalProperty prop) const final {
			using namespace numeric;
			const auto & value = this->accessValueConst ();
			switch (prop) {
			case NumericalProperty::Zero:
				return value == zero (Dimension<T> (value));
			case NumericalProperty::One:
				return value == one (Dimension<T> (value));
			default:
				return false;
			}
		}

		NodeRef derive (Context & c, const Node & node) final {
			const auto dim = Dimension<T> (this->accessValueConst ());
			if (&node == static_cast<const Node *> (this)) {
				return ConstantOne<T>::create (c, dim);
			} else {
				return ConstantZero<T>::create (c, dim);
			}
		}
		bool isDerivable (const Node &) const final { return true; }

	private:
		void compute () final {
			// Mutable is always valid
			failureComputeWasCalled (typeid (*this));
		}
	};

	/** Convert a F value into a R value, with specified dimension.
	 */
	template <typename R, typename F> class Convert : public Value<R> {
	public:
		using Self = Convert;

		static ValueRef<R> create (Context &, NodeRefVec && deps, const Dimension<R> & dim = {}) {
			// Check dependencies
			checkDependenciesNotNull (typeid (Self), deps);
			checkDependencyVectorSize (typeid (Self), deps, 1);
			checkNthDependencyIsValue<F> (typeid (Self), deps, 0);
			// Select node
			if (std::is_same<R, F>::value) {
				return convertRef<Value<R>> (deps[0]);
			} else {
				return std::make_shared<Self> (std::move (deps), dim);
			}
		}

		Convert (NodeRefVec && deps, const Dimension<R> & dim)
		    : Value<R> (std::move (deps)), targetDimension (dim) {}

		NodeRef derive (Context & c, const Node & node) final {
			return Self::create (c, NodeRefVec{this->dependency (0)->derive (c, node)}, targetDimension);
		}
		bool isDerivable (const Node & node) const final {
			return this->dependency (0)->isDerivable (node);
		}

	private:
		void compute () final {
			using namespace numeric;
			const auto & arg = accessValueConstCast<F> (*this->dependency (0));
			this->accessValueMutable () = convert (arg, targetDimension);
		}

		Dimension<R> targetDimension;
	};

	/** Addition of a fixed set of values into R.
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
				return Convert<R, T1>::create (c, NodeRefVec{deps[1]}, dim);
			} else if (!dep_0_zero && dep_0_zero) {
				return Convert<R, T0>::create (c, NodeRefVec{deps[0]}, dim);
			} else {
				return std::make_shared<Self> (std::move (deps), dim);
			}
		}

		CWiseAdd (NodeRefVec && deps, const Dimension<R> & dim)
		    : Value<R> (std::move (deps)), targetDimension (dim) {}

		NodeRef derive (Context & c, const Node & node) {
			constexpr std::size_t n = 2;
			NodeRefVec derivedDeps (n);
			for (std::size_t i = 0; i < n; ++i) {
				derivedDeps[i] = this->dependency (i)->derive (c, node);
			}
			return Self::create (c, std::move (derivedDeps), targetDimension);
		}
		bool isDerivable (const Node & node) const {
			return allDerivable (this->dependencies (), node);
		}

	private:
		void compute () final {
			using namespace numeric;
			auto & result = this->accessValueMutable ();
			cwise (result) = cwise (accessValueConstCast<T0> (*this->dependency (0))) +
			                 cwise (accessValueConstCast<T1> (*this->dependency (1)));
		}

		Dimension<R> targetDimension;
	}; // namespace dataflow

	// Pre compiled instantiations
	extern template class CWiseAdd<double, std::tuple<double, double>>;
	extern template class CWiseAdd<Eigen::VectorXd, std::tuple<Eigen::VectorXd, Eigen::VectorXd>>;
	extern template class CWiseAdd<Eigen::MatrixXd, std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>>;

	/** Addition of any number of T into R.
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

		NodeRef derive (Context & c, const Node & node) {
			const auto n = this->nbDependencies ();
			NodeRefVec derivedDeps (n);
			for (std::size_t i = 0; i < n; ++i) {
				derivedDeps[i] = this->dependency (i)->derive (c, node);
			}
			return Self::create (c, std::move (derivedDeps), targetDimension);
		}
		bool isDerivable (const Node & node) const {
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

	// Pre compiled instantiations
	extern template class CWiseAdd<double, ReductionOf<double>>;
	extern template class CWiseAdd<Eigen::VectorXd, ReductionOf<Eigen::VectorXd>>;
	extern template class CWiseAdd<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

	/** Multiplication of a fixed set of values into R.
	 * Only defined for N = 2 for now. Same problems as CWiseAdd.
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
				return Convert<R, T1>::create (c, NodeRefVec{deps[1]}, dim);
			} else if (!dep_0_one && dep_0_one) {
				return Convert<R, T0>::create (c, NodeRefVec{deps[0]}, dim);
			} else {
				return std::make_shared<Self> (std::move (deps), dim);
			}
		}

		CWiseMul (NodeRefVec && deps, const Dimension<R> & dim)
		    : Value<R> (std::move (deps)), targetDimension (dim) {}

		NodeRef derive (Context & c, const Node & node) {
			constexpr std::size_t n = 2;
			NodeRefVec addDeps (n);
			for (std::size_t i = 0; i < n; ++i) {
				NodeRefVec ithMulDeps = this->dependencies ();
				ithMulDeps[i] = this->dependency (i)->derive (c, node);
				addDeps[i] = Self::create (c, std::move (ithMulDeps), targetDimension);
			}
			return CWiseAdd<R, std::tuple<R, R>>::create (c, std::move (addDeps), targetDimension);
		}
		bool isDerivable (const Node & node) const {
			return allDerivable (this->dependencies (), node);
		}

	private:
		void compute () final {
			using namespace numeric;
			auto & result = this->accessValueMutable ();
			cwise (result) = cwise (accessValueConstCast<T0> (*this->dependency (0))) *
			                 cwise (accessValueConstCast<T1> (*this->dependency (1)));
		}

		Dimension<R> targetDimension;
	}; // namespace dataflow

	// Pre compiled instantiations
	extern template class CWiseMul<double, std::tuple<double, double>>;
	extern template class CWiseMul<Eigen::VectorXd, std::tuple<Eigen::VectorXd, Eigen::VectorXd>>;
	extern template class CWiseMul<Eigen::MatrixXd, std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>>;
	extern template class CWiseMul<Eigen::VectorXd, std::tuple<double, Eigen::VectorXd>>;
	extern template class CWiseMul<Eigen::MatrixXd, std::tuple<double, Eigen::MatrixXd>>;

	/** Multiplication of any number of T into R.
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

		NodeRef derive (Context & c, const Node & node) {
			const auto n = this->nbDependencies ();
			NodeRefVec addDeps (n);
			for (std::size_t i = 0; i < n; ++i) {
				NodeRefVec ithMulDeps = this->dependencies ();
				ithMulDeps[i] = this->dependency (i)->derive (c, node);
				addDeps[i] = Self::create (c, std::move (ithMulDeps), targetDimension);
			}
			return CWiseAdd<R, ReductionOf<R>>::create (c, std::move (addDeps), targetDimension);
		}
		bool isDerivable (const Node & node) const {
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

	// Pre compiled instantiations
	extern template class CWiseMul<double, ReductionOf<double>>;
	extern template class CWiseMul<Eigen::VectorXd, ReductionOf<Eigen::VectorXd>>;
	extern template class CWiseMul<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

	/** Negation of T.
	 */
	template <typename T> class Negate : public Value<T> {
	public:
		using Self = Negate;

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

		Negate (NodeRefVec && deps, const Dimension<T> & dim)
		    : Value<T> (std::move (deps)), targetDimension (dim) {}

		NodeRef derive (Context & c, const Node & node) final {
			return Self::create (c, NodeRefVec{this->dependency (0)->derive (c, node)}, targetDimension);
		}
		bool isDerivable (const Node & node) const final {
			return this->dependency (0)->isDerivable (node);
		}

	private:
		void compute () final {
			this->accessValueMutable () = -accessValueConstCast<T> (*this->dependency (0));
		}

		Dimension<T> targetDimension;
	};

	// Pre compiled instantiations
	extern template class Negate<double>;
	extern template class Negate<Eigen::VectorXd>;
	extern template class Negate<Eigen::MatrixXd>;
} // namespace dataflow
} // namespace bpp

#endif // BPP_NEWPHYL_DATAFLOWNUMERIC2_H
