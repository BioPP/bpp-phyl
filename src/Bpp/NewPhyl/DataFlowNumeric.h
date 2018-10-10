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

#include <Bpp/NewPhyl/ExtendedFloat.h>
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

  /// Empty type representing no dimensions.
  struct NoDimension {};

  std::string to_string (const NoDimension &);
  inline std::size_t hash (const NoDimension &) { return 0; }
  inline bool operator== (const NoDimension &, const NoDimension &) { return true; }
  inline bool operator!= (const NoDimension &, const NoDimension &) { return false; }

  /// Basic matrix dimension type
  struct MatrixDimension {
    Eigen::Index rows{};
    Eigen::Index cols{};

    MatrixDimension () = default;
    MatrixDimension (Eigen::Index rows_, Eigen::Index cols_) : rows (rows_), cols (cols_) {}

    // Get dimensions of any matrix-like eigen object.
    template <typename Derived>
    MatrixDimension (const Eigen::MatrixBase<Derived> & m) : MatrixDimension (m.rows (), m.cols ()) {}
  };

  std::string to_string (const MatrixDimension & dim);
  std::size_t hash (const MatrixDimension & dim);
  inline bool operator== (const MatrixDimension & lhs, const MatrixDimension & rhs) {
    return lhs.rows == rhs.rows && lhs.cols == rhs.cols;
  }
  inline bool operator!= (const MatrixDimension & lhs, const MatrixDimension & rhs) { return !(lhs == rhs); }

  /// Eigen vector are matrices with 1 column.
  inline MatrixDimension vectorDimension (Eigen::Index size) { return {size, 1}; }
  /// Eigen RowVector are matrices with 1 row.
  inline MatrixDimension rowVectorDimension (Eigen::Index size) { return {1, size}; }

  /** Store a dimension for type T.
   * Declared but undefined by default.
   * Specialisations should be defined in the same header declaring the T type.
   * Specialisations should define a constructor from const T & : get the dimension of a T object.
   */
  template <typename T> struct Dimension;

  /** Specialisation of Dimension<T> for floating point types.
   * This is a dummy empty type, required by generic code below.
   */
  template <> struct Dimension<double> : NoDimension {
    Dimension () = default;
    Dimension (const double &) {}
  };
  template <> struct Dimension<float> : NoDimension {
    Dimension () = default;
    Dimension (const float &) {}
  };

  /** Specialisation of Dimension<T> for eigen matrix types.
   * Note that in Eigen, a vector is a matrix with one column.
   * Redirect to MatrixDimension for all eigen matrix variants.
   */
  template <typename T, int Rows, int Cols> struct Dimension<Eigen::Matrix<T, Rows, Cols>> : MatrixDimension {
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
     * matrix -> matrix: copy values, size must match (conversion between eigen dynamic/fixed)
     *
     * Two APIs:
     * r = convert(f, dim); -> convert returns a "converted" result
     * convert (r, f, dim); -> convert does the assignment directly
     */
    template <typename R, typename F, typename = typename std::enable_if<std::is_arithmetic<R>::value>::type>
    R convert (const F & from, const Dimension<R> &) {
      return R (from); // scalar -> scalar
    }
    template <typename T, int Rows, int Cols, typename F,
              typename = typename std::enable_if<std::is_arithmetic<F>::value>::type>
    auto convert (const F & from, const Dimension<Eigen::Matrix<T, Rows, Cols>> & dim)
      -> decltype (Eigen::Matrix<T, Rows, Cols>::Constant (dim.rows, dim.cols, from)) {
      // scalar -> matrix
      return Eigen::Matrix<T, Rows, Cols>::Constant (dim.rows, dim.cols, from);
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
    template <typename Derived> auto cwise (const Eigen::MatrixBase<Derived> & m) -> decltype (m.array ()) {
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

    // Hash of numerical value
    template <typename T, typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
    std::size_t hash (T t) {
      return std::hash<T>{}(t);
    }
    template <typename Derived> std::size_t hash (const Eigen::MatrixBase<Derived> & m) {
      std::size_t seed = 0;
      for (Eigen::Index j = 0; j < m.cols (); ++j) {
        for (Eigen::Index i = 0; i < m.rows (); ++i) {
          combineHash (seed, m (i, j));
        }
      }
      return seed;
    }
  } // namespace numeric

  /******************************************************************************
   * Data flow nodes for those numerical functions.
   *
   * TODO numerical simplification:
   * add(x,x) -> 2*x ? (and similar for mul, ...)
   * all deps constant => return constant ?
   */
  namespace dataflow {
    // Error utils
    [[noreturn]] void failureDeltaNotDerivable (const std::type_info & contextNodeType);
    [[noreturn]] void failureNumericalDerivationNotConfigured ();
    void checkRecreateWithoutDependencies (const std::type_info & contextNodeType, const NodeRefVec & deps);

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
    template <typename F> class SumOfLogarithms;
    template <typename R, typename T0, typename T1> class MatrixProduct;
    template <typename T> class ShiftDelta;
    template <typename T> class CombineDeltaShifted;

    // Utilities
    template <typename Predicate> void removeDependenciesIf (NodeRefVec & deps, Predicate p) {
      auto new_end = std::remove_if (deps.begin (), deps.end (), std::move (p));
      deps.erase (new_end, deps.end ()); // Truncate vector storage
    }

    /** Template struct used to describe a dependency transformation before compute().
     *
     * Transforms allow to generate variants of computation nodes with a readable syntax:
     * MatrixProduct<R, T0, T1> will perform R = T0 * T1.
     * MatrixProduct<R, Transposed<T0>, T1> will perform R = transpose(T0) * T1.
     *
     * This struct is used to implement the transformation.
     * A default case is provided that performs no transformation.
     * Adding a new transformation type consists of declaring a type tag (like Transposed<T>),
     * and specialising this struct for the type tag to implement the transformation.
     * Specialisations should be recursive to allow nesting of transformations (see Transposed<T>).
     *
     * The DepType field gives the real type of the dependency, removing any type tag.
     * For Transposed<T0>, DepType = T0.
     * DepType is used to check the dependency types (Value<DepType>), and casting during computation.
     *
     * transform(const DepType & d) -> ? "performs" the transformation.
     * By default, this just forwards the d reference, doing nothing.
     * For Transposed<T>, this returns an Eigen expression template declaring a transposition.
     *
     * A transformation must not change the result of hasNumericalProperty(), because it is not applied to it,
     * and would make the simplifications wrong.
     * This makes it unlikely that something other than transposition can be implemented with this model.
     *
     * For now, this is only used to implement transposed variants of MatrixProduct, and Convert.
     * This can be extended to other nodes if useful.
     */
    template <typename T> struct NumericalDependencyTransform {
      /// Real type of the dependency: T in the default case.
      using DepType = T;
      /// Transform the DepType value: do nothing in the default case.
      static const DepType & transform (const DepType & d) { return d; }
    };

    /// The T dependency should be transposed before computation.
    template <typename T> struct Transposed;
    /// Implementation for a dependency transposition.
    template <typename T> struct NumericalDependencyTransform<Transposed<T>> {
      // Get DepType by recursion
      using DepType = typename NumericalDependencyTransform<T>::DepType;

      // Perform inner transformation for T, then transpose. Only for Eigen types.
      static auto transform (const DepType & d)
        -> decltype (NumericalDependencyTransform<T>::transform (d).transpose ()) {
        return NumericalDependencyTransform<T>::transform (d).transpose ();
      }
    };

    /** r = 0 for each component.
     * r: T.
     * Value is only created at first use (lazy).
     */
    template <typename T> class ConstantZero : public Value<T> {
    public:
      using Self = ConstantZero;

      static std::shared_ptr<Self> create (Context & c, const Dimension<T> & dim) {
        return cachedAs<Self> (c, std::make_shared<Self> (dim));
      }

      explicit ConstantZero (const Dimension<T> & dim) : Value<T> (NodeRefVec{}), targetDimension_ (dim) {}

      std::string debugInfo () const override { return "targetDim=" + to_string (targetDimension_); }

      bool hasNumericalProperty (NumericalProperty prop) const final {
        switch (prop) {
        case NumericalProperty::Constant:
          return true;
        case NumericalProperty::ConstantZero:
          return true;
        default:
          return false;
        }
      }

      // ConstantZero<T> additional arguments = (targetDimension_).
      bool compareAdditionalArguments (const Node & other) const final {
        const auto * derived = dynamic_cast<const Self *> (&other);
        return derived != nullptr && targetDimension_ == derived->targetDimension_;
      }
      std::size_t hashAdditionalArguments () const final { return hash (targetDimension_); }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<T>::create (c, targetDimension_);
        }
        return this->shared_from_this (); // Return handle to self, as d(0)/dx = 0
      }

      NodeRef recreate (Context &, NodeRefVec && deps) final {
        checkRecreateWithoutDependencies (typeid (Self), deps);
        return this->shared_from_this ();
      }

    private:
      void compute () final {
        using namespace numeric;
        this->accessValueMutable () = zero (targetDimension_);
      }

      Dimension<T> targetDimension_;
    };

    /** r = 1 for each component.
     * r: T.
     * Value is only created at first use (lazy).
     */
    template <typename T> class ConstantOne : public Value<T> {
    public:
      using Self = ConstantOne;

      static std::shared_ptr<Self> create (Context & c, const Dimension<T> & dim) {
        return cachedAs<Self> (c, std::make_shared<Self> (dim));
      }

      explicit ConstantOne (const Dimension<T> & dim) : Value<T> (NodeRefVec{}), targetDimension_ (dim) {}

      std::string debugInfo () const override { return "targetDim=" + to_string (targetDimension_); }

      bool hasNumericalProperty (NumericalProperty prop) const final {
        switch (prop) {
        case NumericalProperty::Constant:
          return true;
        case NumericalProperty::ConstantOne:
          return true;
        default:
          return false;
        }
      }

      // ConstantOne<T> additional arguments = (targetDimension_).
      bool compareAdditionalArguments (const Node & other) const final {
        const auto * derived = dynamic_cast<const Self *> (&other);
        return derived != nullptr && targetDimension_ == derived->targetDimension_;
      }
      std::size_t hashAdditionalArguments () const final { return hash (targetDimension_); }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<T>::create (c, targetDimension_);
        }
        return ConstantZero<T>::create (c, targetDimension_);
      }

      NodeRef recreate (Context &, NodeRefVec && deps) final {
        checkRecreateWithoutDependencies (typeid (Self), deps);
        return this->shared_from_this ();
      }

    private:
      void compute () final {
        using namespace numeric;
        this->accessValueMutable () = one (targetDimension_);
      }

      Dimension<T> targetDimension_;
    };

    /** r = const.
     * r: T.
     * Value is set at construction, and cannot change.
     * Supports derivation.
     */
    template <typename T> class NumericConstant : public Value<T> {
    public:
      using Self = NumericConstant;

      template <typename... Args> static std::shared_ptr<Self> create (Context & c, Args &&... args) {
        return cachedAs<Self> (c, std::make_shared<Self> (std::forward<Args> (args)...));
      }

      /// Sets an initial value constructed with the given arguments.
      template <typename... Args>
      explicit NumericConstant (Args &&... args) : Value<T> (NodeRefVec{}, std::forward<Args> (args)...) {
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
        case NumericalProperty::ConstantZero:
          return value == zero (Dimension<T> (value));
        case NumericalProperty::ConstantOne:
          return value == one (Dimension<T> (value));
        case NumericalProperty::ConstantIdentity:
          return isIdentity (value);
        default:
          return false;
        }
      }

      // NumericConstant<T> additional arguments = (value).
      bool compareAdditionalArguments (const Node & other) const final {
        const auto * derived = dynamic_cast<const Self *> (&other);
        return derived != nullptr && this->accessValueConst () == derived->accessValueConst ();
      }
      std::size_t hashAdditionalArguments () const final {
        using namespace numeric;
        return hash (this->accessValueConst ());
      }

      NodeRef derive (Context & c, const Node & node) final {
        const auto dim = Dimension<T> (this->accessValueConst ());
        if (&node == this) {
          return ConstantOne<T>::create (c, dim);
        }
        return ConstantZero<T>::create (c, dim);
      }

      NodeRef recreate (Context &, NodeRefVec && deps) final {
        checkRecreateWithoutDependencies (typeid (Self), deps);
        return this->shared_from_this ();
      }

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
     * This node has no Context merging support: mutable nodes are always different.
     */
    template <typename T> class NumericMutable : public Value<T> {
    public:
      using Self = NumericMutable;

      template <typename... Args> static std::shared_ptr<Self> create (Context &, Args &&... args) {
        return std::make_shared<Self> (std::forward<Args> (args)...);
      }

      /// Sets an initial value constructed with the given arguments.
      template <typename... Args>
      explicit NumericMutable (Args &&... args) : Value<T> (NodeRefVec{}, std::forward<Args> (args)...) {
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

      NodeRef derive (Context & c, const Node & node) final {
        const auto dim = Dimension<T> (this->accessValueConst ());
        if (&node == this) {
          return ConstantOne<T>::create (c, dim);
        }
        return ConstantZero<T>::create (c, dim);
      }

      NodeRef recreate (Context &, NodeRefVec && deps) final {
        checkRecreateWithoutDependencies (typeid (Self), deps);
        return this->shared_from_this ();
      }

    private:
      void compute () final {
        // Mutable is always valid
        failureComputeWasCalled (typeid (*this));
      }
    };

    /** r = f.
     * r: R.
     * f: F, allows NumericalDependencyTransform.
     * Convert from F to R type, semantics of numeric::convert.
     */
    template <typename R, typename F> class Convert : public Value<R> {
    public:
      using Self = Convert;
      using DepF = typename NumericalDependencyTransform<F>::DepType;

      static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim) {
        // Check dependencies
        checkDependenciesNotNull (typeid (Self), deps);
        checkDependencyVectorSize (typeid (Self), deps, 1);
        checkNthDependencyIsValue<DepF> (typeid (Self), deps, 0);
        // Select node
        if (std::is_same<R, F>::value) {
          return convertRef<Value<R>> (deps[0]);
        } else if (deps[0]->hasNumericalProperty (NumericalProperty::ConstantZero)) {
          return ConstantZero<R>::create (c, dim);
        } else if (deps[0]->hasNumericalProperty (NumericalProperty::ConstantOne)) {
          return ConstantOne<R>::create (c, dim);
        } else {
          return cachedAs<Value<R>> (c, std::make_shared<Self> (std::move (deps), dim));
        }
      }

      Convert (NodeRefVec && deps, const Dimension<R> & dim)
        : Value<R> (std::move (deps)), targetDimension_ (dim) {}

      std::string debugInfo () const override {
        using namespace numeric;
        return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
      }

      // Convert<T> additional arguments = ().
      bool compareAdditionalArguments (const Node & other) const final {
        return dynamic_cast<const Self *> (&other) != nullptr;
      }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<R>::create (c, targetDimension_);
        }
        return Self::create (c, {this->dependency (0)->derive (c, node)}, targetDimension_);
      }

      NodeRef recreate (Context & c, NodeRefVec && deps) final {
        return Self::create (c, std::move (deps), targetDimension_);
      }

    private:
      void compute () final {
        using namespace numeric;
        auto & result = this->accessValueMutable ();
        const auto & arg = accessValueConstCast<DepF> (*this->dependency (0));
        convert (result, NumericalDependencyTransform<F>::transform (arg), targetDimension_);
      }

      Dimension<R> targetDimension_;
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
    template <typename R, typename T0, typename T1> class CWiseAdd<R, std::tuple<T0, T1>> : public Value<R> {
    public:
      using Self = CWiseAdd;

      static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim) {
        // Check dependencies
        checkDependenciesNotNull (typeid (Self), deps);
        checkDependencyVectorSize (typeid (Self), deps, 2);
        checkNthDependencyIsValue<T0> (typeid (Self), deps, 0);
        checkNthDependencyIsValue<T1> (typeid (Self), deps, 1);
        // Select node implementation
        bool zeroDep0 = deps[0]->hasNumericalProperty (NumericalProperty::ConstantZero);
        bool zeroDep1 = deps[1]->hasNumericalProperty (NumericalProperty::ConstantZero);
        if (zeroDep0 && zeroDep1) {
          return ConstantZero<R>::create (c, dim);
        } else if (zeroDep0 && !zeroDep1) {
          return Convert<R, T1>::create (c, {deps[1]}, dim);
        } else if (!zeroDep0 && zeroDep1) {
          return Convert<R, T0>::create (c, {deps[0]}, dim);
        } else {
          return cachedAs<Value<R>> (c, std::make_shared<Self> (std::move (deps), dim));
        }
      }

      CWiseAdd (NodeRefVec && deps, const Dimension<R> & dim)
        : Value<R> (std::move (deps)), targetDimension_ (dim) {}

      std::string debugInfo () const override {
        using namespace numeric;
        return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
      }

      // Convert<T> additional arguments = ().
      bool compareAdditionalArguments (const Node & other) const final {
        return dynamic_cast<const Self *> (&other) != nullptr;
      }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<R>::create (c, targetDimension_);
        }
        constexpr std::size_t n = 2;
        NodeRefVec derivedDeps (n);
        for (std::size_t i = 0; i < n; ++i) {
          derivedDeps[i] = this->dependency (i)->derive (c, node);
        }
        return Self::create (c, std::move (derivedDeps), targetDimension_);
      }

      NodeRef recreate (Context & c, NodeRefVec && deps) final {
        return Self::create (c, std::move (deps), targetDimension_);
      }

    private:
      void compute () final {
        using namespace numeric;
        auto & result = this->accessValueMutable ();
        const auto & x0 = accessValueConstCast<T0> (*this->dependency (0));
        const auto & x1 = accessValueConstCast<T1> (*this->dependency (1));
        cwise (result) = cwise (x0) + cwise (x1);
      }

      Dimension<R> targetDimension_;
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

      static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim) {
        // Check dependencies
        checkDependenciesNotNull (typeid (Self), deps);
        checkDependencyRangeIsValue<T> (typeid (Self), deps, 0, deps.size ());
        // Remove 0s from deps
        removeDependenciesIf (deps, [](const NodeRef & dep) {
          return dep->hasNumericalProperty (NumericalProperty::ConstantZero);
        });
        // Select node implementation
        if (deps.size () == 0) {
          return ConstantZero<R>::create (c, dim);
        } else if (deps.size () == 1) {
          return Convert<R, T>::create (c, std::move (deps), dim);
        } else if (deps.size () == 2) {
          return CWiseAdd<R, std::tuple<T, T>>::create (c, std::move (deps), dim);
        } else {
          return cachedAs<Value<R>> (c, std::make_shared<Self> (std::move (deps), dim));
        }
      }

      CWiseAdd (NodeRefVec && deps, const Dimension<R> & dim)
        : Value<R> (std::move (deps)), targetDimension_ (dim) {}

      std::string debugInfo () const override {
        using namespace numeric;
        return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
      }

      // CWiseAdd additional arguments = ().
      bool compareAdditionalArguments (const Node & other) const final {
        return dynamic_cast<const Self *> (&other) != nullptr;
      }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<R>::create (c, targetDimension_);
        }
        const auto n = this->nbDependencies ();
        NodeRefVec derivedDeps (n);
        for (std::size_t i = 0; i < n; ++i) {
          derivedDeps[i] = this->dependency (i)->derive (c, node);
        }
        return Self::create (c, std::move (derivedDeps), targetDimension_);
      }

      NodeRef recreate (Context & c, NodeRefVec && deps) final {
        return Self::create (c, std::move (deps), targetDimension_);
      }

    private:
      void compute () final {
        using namespace numeric;
        auto & result = this->accessValueMutable ();
        result = zero (targetDimension_);
        for (const auto & depNodeRef : this->dependencies ()) {
          cwise (result) += cwise (accessValueConstCast<T> (*depNodeRef));
        }
      }

      Dimension<R> targetDimension_;
    };

    /** r = x0 * x1 for each component.
     * r: R.
     * x0: T0.
     * x1: T1.
     *
     * Values converted to R with the semantics of numeric::convert.
     * Only defined for N = 2 for now (same constraints as CWiseAdd for genericity).
     */
    template <typename R, typename T0, typename T1> class CWiseMul<R, std::tuple<T0, T1>> : public Value<R> {
    public:
      using Self = CWiseMul;

      static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim) {
        // Check dependencies
        checkDependenciesNotNull (typeid (Self), deps);
        checkDependencyVectorSize (typeid (Self), deps, 2);
        checkNthDependencyIsValue<T0> (typeid (Self), deps, 0);
        checkNthDependencyIsValue<T1> (typeid (Self), deps, 1);
        // Return 0 if any 0.
        if (std::any_of (deps.begin (), deps.end (), [](const NodeRef & dep) {
              return dep->hasNumericalProperty (NumericalProperty::ConstantZero);
            })) {
          return ConstantZero<R>::create (c, dim);
        }
        // Select node implementation
        bool oneDep0 = deps[0]->hasNumericalProperty (NumericalProperty::ConstantOne);
        bool oneDep1 = deps[1]->hasNumericalProperty (NumericalProperty::ConstantOne);
        if (oneDep0 && oneDep1) {
          return ConstantOne<R>::create (c, dim);
        } else if (oneDep0 && !oneDep1) {
          return Convert<R, T1>::create (c, {deps[1]}, dim);
        } else if (!oneDep0 && oneDep1) {
          return Convert<R, T0>::create (c, {deps[0]}, dim);
        } else {
          return cachedAs<Value<R>> (c, std::make_shared<Self> (std::move (deps), dim));
        }
      }

      CWiseMul (NodeRefVec && deps, const Dimension<R> & dim)
        : Value<R> (std::move (deps)), targetDimension_ (dim) {}

      std::string debugInfo () const override {
        using namespace numeric;
        return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
      }

      // CWiseAdd additional arguments = ().
      bool compareAdditionalArguments (const Node & other) const final {
        return dynamic_cast<const Self *> (&other) != nullptr;
      }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<R>::create (c, targetDimension_);
        }
        constexpr std::size_t n = 2;
        NodeRefVec addDeps (n);
        for (std::size_t i = 0; i < n; ++i) {
          NodeRefVec ithMulDeps = this->dependencies ();
          ithMulDeps[i] = this->dependency (i)->derive (c, node);
          addDeps[i] = Self::create (c, std::move (ithMulDeps), targetDimension_);
        }
        return CWiseAdd<R, std::tuple<R, R>>::create (c, std::move (addDeps), targetDimension_);
      }

      NodeRef recreate (Context & c, NodeRefVec && deps) final {
        return Self::create (c, std::move (deps), targetDimension_);
      }

    private:
      void compute () final {
        using namespace numeric;
        auto & result = this->accessValueMutable ();
        const auto & x0 = accessValueConstCast<T0> (*this->dependency (0));
        const auto & x1 = accessValueConstCast<T1> (*this->dependency (1));
        cwise (result) = cwise (x0) * cwise (x1);
      }

      Dimension<R> targetDimension_;
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

      static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim) {
        // Check dependencies
        checkDependenciesNotNull (typeid (Self), deps);
        checkDependencyRangeIsValue<T> (typeid (Self), deps, 0, deps.size ());
        // If there is a 0 return 0.
        if (std::any_of (deps.begin (), deps.end (), [](const NodeRef & dep) {
              return dep->hasNumericalProperty (NumericalProperty::ConstantZero);
            })) {
          return ConstantZero<R>::create (c, dim);
        }
        // Remove 1s from deps
        removeDependenciesIf (deps, [](const NodeRef & dep) {
          return dep->hasNumericalProperty (NumericalProperty::ConstantOne);
        });
        // Select node implementation
        if (deps.size () == 0) {
          return ConstantOne<R>::create (c, dim);
        } else if (deps.size () == 1) {
          return Convert<R, T>::create (c, std::move (deps), dim);
        } else if (deps.size () == 2) {
          return CWiseMul<R, std::tuple<T, T>>::create (c, std::move (deps), dim);
        } else {
          return cachedAs<Value<R>> (c, std::make_shared<Self> (std::move (deps), dim));
        }
      }

      CWiseMul (NodeRefVec && deps, const Dimension<R> & dim)
        : Value<R> (std::move (deps)), targetDimension_ (dim) {}

      std::string debugInfo () const override {
        using namespace numeric;
        return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
      }

      // CWiseAdd additional arguments = ().
      bool compareAdditionalArguments (const Node & other) const final {
        return dynamic_cast<const Self *> (&other) != nullptr;
      }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<R>::create (c, targetDimension_);
        }
        const auto n = this->nbDependencies ();
        NodeRefVec addDeps (n);
        for (std::size_t i = 0; i < n; ++i) {
          NodeRefVec ithMulDeps = this->dependencies ();
          ithMulDeps[i] = this->dependency (i)->derive (c, node);
          addDeps[i] = Self::create (c, std::move (ithMulDeps), targetDimension_);
        }
        return CWiseAdd<R, ReductionOf<R>>::create (c, std::move (addDeps), targetDimension_);
      }

      NodeRef recreate (Context & c, NodeRefVec && deps) final {
        return Self::create (c, std::move (deps), targetDimension_);
      }

    private:
      void compute () final {
        using namespace numeric;
        auto & result = this->accessValueMutable ();
        result = one (targetDimension_);
        for (const auto & depNodeRef : this->dependencies ()) {
          cwise (result) *= cwise (accessValueConstCast<T> (*depNodeRef));
        }
      }

      Dimension<R> targetDimension_;
    };

    /** r = -x, for each component.
     * r, x: T.
     */
    template <typename T> class CWiseNegate : public Value<T> {
    public:
      using Self = CWiseNegate;

      static ValueRef<T> create (Context & c, NodeRefVec && deps, const Dimension<T> & dim) {
        // Check dependencies
        checkDependenciesNotNull (typeid (Self), deps);
        checkDependencyVectorSize (typeid (Self), deps, 1);
        checkNthDependencyIsValue<T> (typeid (Self), deps, 0);
        // Select node
        if (deps[0]->hasNumericalProperty (NumericalProperty::ConstantZero)) {
          return ConstantZero<T>::create (c, dim);
        } else {
          return cachedAs<Value<T>> (c, std::make_shared<Self> (std::move (deps), dim));
        }
      }

      CWiseNegate (NodeRefVec && deps, const Dimension<T> & dim)
        : Value<T> (std::move (deps)), targetDimension_ (dim) {}

      std::string debugInfo () const override {
        using namespace numeric;
        return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
      }

      // CWiseNegate additional arguments = ().
      bool compareAdditionalArguments (const Node & other) const final {
        return dynamic_cast<const Self *> (&other) != nullptr;
      }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<T>::create (c, targetDimension_);
        }
        return Self::create (c, {this->dependency (0)->derive (c, node)}, targetDimension_);
      }

      NodeRef recreate (Context & c, NodeRefVec && deps) final {
        return Self::create (c, std::move (deps), targetDimension_);
      }

    private:
      void compute () final {
        using namespace numeric;
        auto & result = this->accessValueMutable ();
        const auto & x = accessValueConstCast<T> (*this->dependency (0));
        cwise (result) = -cwise (x);
      }

      Dimension<T> targetDimension_;
    };

    /** r = 1/x for each component.
     * r, x: T.
     */
    template <typename T> class CWiseInverse : public Value<T> {
    public:
      using Self = CWiseInverse;

      static ValueRef<T> create (Context & c, NodeRefVec && deps, const Dimension<T> & dim) {
        // Check dependencies
        checkDependenciesNotNull (typeid (Self), deps);
        checkDependencyVectorSize (typeid (Self), deps, 1);
        checkNthDependencyIsValue<T> (typeid (Self), deps, 0);
        // Select node
        if (deps[0]->hasNumericalProperty (NumericalProperty::ConstantOne)) {
          return ConstantOne<T>::create (c, dim);
        } else {
          return cachedAs<Value<T>> (c, std::make_shared<Self> (std::move (deps), dim));
        }
      }

      CWiseInverse (NodeRefVec && deps, const Dimension<T> & dim)
        : Value<T> (std::move (deps)), targetDimension_ (dim) {}

      std::string debugInfo () const override {
        using namespace numeric;
        return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
      }

      // CWiseInverse additional arguments = ().
      bool compareAdditionalArguments (const Node & other) const final {
        return dynamic_cast<const Self *> (&other) != nullptr;
      }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<T>::create (c, targetDimension_);
        }
        // -1/x^2 * x'
        const auto & dep = this->dependency (0);
        return CWiseMul<T, std::tuple<T, T>>::create (
          c, {CWiseConstantPow<T>::create (c, {dep}, -2., -1., targetDimension_), dep->derive (c, node)},
          targetDimension_);
      }

      NodeRef recreate (Context & c, NodeRefVec && deps) final {
        return Self::create (c, std::move (deps), targetDimension_);
      }

    private:
      void compute () final {
        using namespace numeric;
        auto & result = this->accessValueMutable ();
        const auto & x = accessValueConstCast<T> (*this->dependency (0));
        cwise (result) = inverse (cwise (x));
      }

      Dimension<T> targetDimension_;
    };

    /** r = factor * pow (x, exponent) for each component.
     * r, x: T.
     * exponent, factor: double (constant parameter of the node).
     */
    template <typename T> class CWiseConstantPow : public Value<T> {
    public:
      using Self = CWiseConstantPow;

      static ValueRef<T> create (Context & c, NodeRefVec && deps, double exponent, double factor,
                                 const Dimension<T> & dim) {
        // Check dependencies
        checkDependenciesNotNull (typeid (Self), deps);
        checkDependencyVectorSize (typeid (Self), deps, 1);
        checkNthDependencyIsValue<T> (typeid (Self), deps, 0);
        // Select node implementation
        if (exponent == 0. || deps[0]->hasNumericalProperty (NumericalProperty::ConstantOne)) {
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
            {NumericConstant<double>::create (c, factor), CWiseInverse<T>::create (c, std::move (deps), dim)},
            dim);
        } else {
          return cachedAs<Value<T>> (c, std::make_shared<Self> (std::move (deps), exponent, factor, dim));
        }
      }

      CWiseConstantPow (NodeRefVec && deps, double exponent, double factor, const Dimension<T> & dim)
        : Value<T> (std::move (deps)), targetDimension_ (dim), exponent_ (exponent), factor_ (factor) {}

      std::string debugInfo () const override {
        using namespace numeric;
        return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_) +
               " exponent=" + std::to_string (exponent_) + " factor=" + std::to_string (factor_);
      }

      // CWiseConstantPow additional arguments = (exponent_, factor_).
      bool compareAdditionalArguments (const Node & other) const final {
        const auto * derived = dynamic_cast<const Self *> (&other);
        return derived != nullptr && exponent_ == derived->exponent_ && factor_ == derived->factor_;
      }
      std::size_t hashAdditionalArguments () const final {
        std::size_t seed = 0;
        combineHash (seed, exponent_);
        combineHash (seed, factor_);
        return seed;
      }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<T>::create (c, targetDimension_);
        }
        // factor * (exponent * x^(exponent - 1)) * x'
        const auto & dep = this->dependency (0);
        auto dpow = Self::create (c, {dep}, exponent_ - 1., factor_ * exponent_, targetDimension_);
        return CWiseMul<T, std::tuple<T, T>>::create (c, {dpow, dep->derive (c, node)}, targetDimension_);
      }

      NodeRef recreate (Context & c, NodeRefVec && deps) final {
        return Self::create (c, std::move (deps), exponent_, factor_, targetDimension_);
      }

    private:
      void compute () final {
        using namespace numeric;
        auto & result = this->accessValueMutable ();
        const auto & x = accessValueConstCast<T> (*this->dependency (0));
        cwise (result) = factor_ * pow (cwise (x), exponent_);
      }

      Dimension<T> targetDimension_;
      double exponent_;
      double factor_;
    };

    /** r = x0 * x1 (dot product).
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
        if (deps[0]->hasNumericalProperty (NumericalProperty::ConstantZero) &&
            deps[1]->hasNumericalProperty (NumericalProperty::ConstantZero)) {
          return ConstantZero<double>::create (c, Dimension<double> ());
        } else {
          return cachedAs<Value<double>> (c, std::make_shared<Self> (std::move (deps)));
        }
      }

      ScalarProduct (NodeRefVec && deps) : Value<double> (std::move (deps)) {}

      std::string debugInfo () const override {
        using namespace numeric;
        return debug (this->accessValueConst ());
      }

      // ScalarProduct additional arguments = ().
      bool compareAdditionalArguments (const Node & other) const final {
        return dynamic_cast<const Self *> (&other) != nullptr;
      }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<double>::create (c, Dimension<double> ());
        }
        const auto & x0 = this->dependency (0);
        const auto & x1 = this->dependency (1);
        auto dx0_prod = Self::create (c, {x0->derive (c, node), x1});
        auto dx1_prod = Self::create (c, {x0, x1->derive (c, node)});
        return CWiseAdd<double, std::tuple<double, double>>::create (c, {dx0_prod, dx1_prod},
                                                                     Dimension<double> ());
      }

      NodeRef recreate (Context & c, NodeRefVec && deps) final { return Self::create (c, std::move (deps)); }

    private:
      void compute () final {
        auto & result = this->accessValueMutable ();
        const auto & x0 = accessValueConstCast<T0> (*this->dependency (0));
        const auto & x1 = accessValueConstCast<T1> (*this->dependency (1));
        result = x0.dot (x1); // Using lhs.dot(rhs) method from Eigen only
      }
    };

    /** r = sum_{v in m} log (v).
     * r: double.
     * m: F (matrix-like type).
     *
     * The node has no dimension (double).
     * The dimension of m should be provided for derivation.
     */
    template <typename F> class SumOfLogarithms : public Value<double> {
    public:
      using Self = SumOfLogarithms;

      static ValueRef<double> create (Context & c, NodeRefVec && deps, const Dimension<F> & mDim) {
        checkDependenciesNotNull (typeid (Self), deps);
        checkDependencyVectorSize (typeid (Self), deps, 1);
        checkNthDependencyIsValue<F> (typeid (Self), deps, 0);
        return cachedAs<Value<double>> (c, std::make_shared<Self> (std::move (deps), mDim));
      }

      SumOfLogarithms (NodeRefVec && deps, const Dimension<F> & mDim)
        : Value<double> (std::move (deps)), mTargetDimension_ (mDim) {}

      std::string debugInfo () const override {
        using namespace numeric;
        return debug (this->accessValueConst ());
      }

      // SumOfLogarithms additional arguments = ().
      bool compareAdditionalArguments (const Node & other) const final {
        return dynamic_cast<const Self *> (&other) != nullptr;
      }

      NodeRef derive (Context & c, const Node & node) final {
        const auto & m = this->dependency (0);
        auto dm_dn = m->derive (c, node);
        auto m_inverse = CWiseInverse<F>::create (c, {m}, mTargetDimension_);
        return ScalarProduct<F, F>::create (c, {std::move (dm_dn), std::move (m_inverse)});
      }

      NodeRef recreate (Context & c, NodeRefVec && deps) final {
        return Self::create (c, std::move (deps), mTargetDimension_);
      }

    private:
      void compute () final {
        auto & result = this->accessValueMutable ();
        const auto & m = accessValueConstCast<F> (*this->dependency (0));
        const ExtendedFloat product = m.unaryExpr ([](double d) {
                                         ExtendedFloat ef{d};
                                         ef.normalize_small ();
                                         return ef;
                                       })
                                        .redux ([](const ExtendedFloat & lhs, const ExtendedFloat & rhs) {
                                          auto r = denorm_mul (lhs, rhs);
                                          r.normalize_small ();
                                          return r;
                                        });
        result = log (product);
      }

      Dimension<F> mTargetDimension_;
    };

    /** r = x0 * x1 (matrix product).
     * r: R (matrix).
     * x0: T0 (matrix), allows NumericalDependencyTransform.
     * x1: T1 (matrix), allows NumericalDependencyTransform.
     */
    template <typename R, typename T0, typename T1> class MatrixProduct : public Value<R> {
    public:
      using Self = MatrixProduct;
      using DepT0 = typename NumericalDependencyTransform<T0>::DepType;
      using DepT1 = typename NumericalDependencyTransform<T1>::DepType;

      static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim) {
        // Check dependencies
        checkDependenciesNotNull (typeid (Self), deps);
        checkDependencyVectorSize (typeid (Self), deps, 2);
        checkNthDependencyIsValue<DepT0> (typeid (Self), deps, 0);
        checkNthDependencyIsValue<DepT1> (typeid (Self), deps, 1);
        // Return 0 if any 0.
        if (std::any_of (deps.begin (), deps.end (), [](const NodeRef & dep) {
              return dep->hasNumericalProperty (NumericalProperty::ConstantZero);
            })) {
          return ConstantZero<R>::create (c, dim);
        }
        // Select node implementation
        bool identityDep0 = deps[0]->hasNumericalProperty (NumericalProperty::ConstantIdentity);
        bool identityDep1 = deps[1]->hasNumericalProperty (NumericalProperty::ConstantIdentity);
        if (identityDep0 && identityDep1) {
          // No specific class for Identity
          using namespace numeric;
          return NumericConstant<R>::create (c, identity (dim));
        } else if (identityDep0 && !identityDep1) {
          return Convert<R, T1>::create (c, {deps[1]}, dim);
        } else if (!identityDep0 && identityDep1) {
          return Convert<R, T0>::create (c, {deps[0]}, dim);
        } else {
          return cachedAs<Value<R>> (c, std::make_shared<Self> (std::move (deps), dim));
        }
      }

      MatrixProduct (NodeRefVec && deps, const Dimension<R> & dim)
        : Value<R> (std::move (deps)), targetDimension_ (dim) {}

      std::string debugInfo () const override {
        using namespace numeric;
        return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
      }

      // MatrixProduct additional arguments = ().
      bool compareAdditionalArguments (const Node & other) const final {
        return dynamic_cast<const Self *> (&other) != nullptr;
      }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<R>::create (c, targetDimension_);
        }
        const auto & x0 = this->dependency (0);
        const auto & x1 = this->dependency (1);
        auto dx0_prod = Self::create (c, {x0->derive (c, node), x1}, targetDimension_);
        auto dx1_prod = Self::create (c, {x0, x1->derive (c, node)}, targetDimension_);
        return CWiseAdd<R, std::tuple<R, R>>::create (c, {dx0_prod, dx1_prod}, targetDimension_);
      }

      NodeRef recreate (Context & c, NodeRefVec && deps) final {
        return Self::create (c, std::move (deps), targetDimension_);
      }

    private:
      void compute () final {
        auto & result = this->accessValueMutable ();
        const auto & x0 = accessValueConstCast<DepT0> (*this->dependency (0));
        const auto & x1 = accessValueConstCast<DepT1> (*this->dependency (1));
        result.noalias () =
          NumericalDependencyTransform<T0>::transform (x0) * NumericalDependencyTransform<T1>::transform (x1);
      }

      Dimension<R> targetDimension_;
    };

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

      static ValueRef<T> create (Context & c, NodeRefVec && deps, int n, const Dimension<T> & dim) {
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
        if (n == 0 || delta->hasNumericalProperty (NumericalProperty::ConstantZero)) {
          return convertRef<Value<T>> (x);
        } else {
          return cachedAs<Value<T>> (c, std::make_shared<Self> (std::move (deps), n, dim));
        }
      }

      ShiftDelta (NodeRefVec && deps, int n, const Dimension<T> & dim)
        : Value<T> (std::move (deps)), targetDimension_ (dim), n_ (n) {}

      std::string debugInfo () const override {
        using namespace numeric;
        return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_) +
               " n=" + std::to_string (n_);
      }

      // ShiftDelta additional arguments = (n_).
      bool compareAdditionalArguments (const Node & other) const final {
        const auto * derived = dynamic_cast<const Self *> (&other);
        return derived != nullptr && n_ == derived->n_;
      }
      std::size_t hashAdditionalArguments () const final {
        std::size_t seed = 0;
        combineHash (seed, n_);
        return seed;
      }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<T>::create (c, targetDimension_);
        }
        auto & delta = this->dependency (0);
        auto & x = this->dependency (1);
        return Self::create (c, {delta->derive (c, node), x->derive (c, node)}, n_, targetDimension_);
      }

      NodeRef recreate (Context & c, NodeRefVec && deps) final {
        return Self::create (c, std::move (deps), n_, targetDimension_);
      }

      int getN () const { return n_; }

    private:
      void compute () final {
        using namespace numeric;
        auto & result = this->accessValueMutable ();
        const auto & delta = accessValueConstCast<double> (*this->dependency (0));
        const auto & x = accessValueConstCast<T> (*this->dependency (1));
        cwise (result) = n_ * delta + cwise (x);
      }

      Dimension<T> targetDimension_;
      int n_;
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
     * Note that in the whole class, coeffs[i] is the coefficient of deps[i + 1] !
     */
    template <typename T> class CombineDeltaShifted : public Value<T> {
    public:
      using Self = CombineDeltaShifted;

      static ValueRef<T> create (Context & c, NodeRefVec && deps, int n, std::vector<double> && coeffs,
                                 const Dimension<T> & dim) {
        // Check dependencies
        checkDependenciesNotNull (typeid (Self), deps);
        checkDependencyVectorSize (typeid (Self), deps, 1 + coeffs.size ());
        checkNthDependencyIsValue<double> (typeid (Self), deps, 0);
        checkDependencyRangeIsValue<T> (typeid (Self), deps, 1, deps.size ());
        //
        auto cleanAndCreateNode = [&c, &dim](NodeRefVec && deps2, int n2,
                                             std::vector<double> && coeffs2) -> ValueRef<T> {
          // Clean deps : remove constant 0, of deps with a 0 factor.
          for (std::size_t i = 0; i < coeffs2.size ();) {
            if (coeffs2[i] == 0. || deps2[i + 1]->hasNumericalProperty (NumericalProperty::ConstantZero)) {
              coeffs2.erase (coeffs2.begin () + std::ptrdiff_t (i));
              deps2.erase (deps2.begin () + std::ptrdiff_t (i + 1));
            } else {
              ++i;
            }
          }
          // Final node selection
          if (coeffs2.empty ()) {
            return ConstantZero<T>::create (c, dim);
          } else {
            return cachedAs<Value<T>> (
              c, std::make_shared<Self> (std::move (deps2), n2, std::move (coeffs2), dim));
          }
        };
        // Detect if we can merge this node with its dependencies
        const auto & delta = deps[0];
        auto isSelfWithSameDelta = [&delta](const NodeRef & dep) {
          return dynamic_cast<const Self *> (dep.get ()) != nullptr && dep->dependency (0) == delta;
        };
        if (!coeffs.empty () && std::all_of (deps.begin () + 1, deps.end (), isSelfWithSameDelta)) {
          const auto depN = static_cast<const Self &> (*deps[1]).getN ();
          auto useSameNasDep1 = [depN](const NodeRef & dep) {
            return static_cast<const Self &> (*dep).getN () == depN;
          };
          if (std::all_of (deps.begin () + 2, deps.end (), useSameNasDep1)) {
            /* Merge with dependencies because they use the same delta and a common N.
             *
             * V(this) = 1/delta^n * sum_i V(deps[i]) * coeffs[i].
             * V(deps[i]) = 1/delta^depN * sum_j V(depDeps[i][j]) * depCoeffs[i][j].
             * V(this) = 1/delta^(n + depN) * sum_ij V(depDeps[i][j]) * coeffs[i] * depCoeffs[i][j].
             * And simplify the sum by summing coeffs for each unique depDeps[i][j].
             */
            NodeRefVec mergedDeps{delta};
            std::vector<double> mergedCoeffs;
            for (std::size_t i = 0; i < coeffs.size (); ++i) {
              const auto & dep = static_cast<const Self &> (*deps[i + 1]);
              const auto & depCoeffs = dep.getCoeffs ();
              for (std::size_t j = 0; j < depCoeffs.size (); ++j) {
                const auto & subDep = dep.dependency (j + 1);
                auto it = std::find (mergedDeps.begin () + 1, mergedDeps.end (), subDep);
                if (it != mergedDeps.end ()) {
                  // Found
                  const auto subDepIndexInMerged =
                    static_cast<std::size_t> (std::distance (mergedDeps.begin () + 1, it));
                  mergedCoeffs[subDepIndexInMerged] += coeffs[i] * depCoeffs[j];
                } else {
                  // Not found
                  mergedDeps.emplace_back (subDep);
                  mergedCoeffs.emplace_back (coeffs[i] * depCoeffs[j]);
                }
              }
            }
            return cleanAndCreateNode (std::move (mergedDeps), n + depN, std::move (mergedCoeffs));
          }
        }
        // If not able to merge
        return cleanAndCreateNode (std::move (deps), n, std::move (coeffs));
      }

      CombineDeltaShifted (NodeRefVec && deps, int n, std::vector<double> && coeffs, const Dimension<T> & dim)
        : Value<T> (std::move (deps)), targetDimension_ (dim), coeffs_ (std::move (coeffs)), n_ (n) {
        assert (this->coeffs_.size () + 1 == this->nbDependencies ());
      }

      std::string debugInfo () const override {
        using namespace numeric;
        std::string s = debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_) +
                        " n=" + std::to_string (n_) + " coeffs={";
        if (!coeffs_.empty ()) {
          s += std::to_string (coeffs_[0]);
          for (std::size_t i = 1; i < coeffs_.size (); ++i) {
            s += ';' + std::to_string (coeffs_[i]);
          }
        }
        s += '}';
        return s;
      }

      // CombineDeltaShifted additional arguments = (n_, coeffs_).
      bool compareAdditionalArguments (const Node & other) const final {
        const auto * derived = dynamic_cast<const Self *> (&other);
        return derived != nullptr && n_ == derived->n_ && coeffs_ == derived->coeffs_;
      }
      std::size_t hashAdditionalArguments () const final {
        std::size_t seed = 0;
        combineHash (seed, n_);
        for (const auto d : coeffs_) {
          combineHash (seed, d);
        }
        return seed;
      }

      NodeRef derive (Context & c, const Node & node) final {
        if (&node == this) {
          return ConstantOne<T>::create (c, targetDimension_);
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
        return Self::create (c, std::move (derivedDeps), n_, std::vector<double>{coeffs_}, targetDimension_);
      }

      NodeRef recreate (Context & c, NodeRefVec && deps) final {
        return Self::create (c, std::move (deps), n_, std::vector<double>{coeffs_}, targetDimension_);
      }

      const std::vector<double> & getCoeffs () const { return coeffs_; }
      int getN () const { return n_; }

    private:
      void compute () final {
        using namespace numeric;
        auto & result = this->accessValueMutable ();
        const auto & delta = accessValueConstCast<double> (*this->dependency (0));
        const double lambda = pow (delta, -n_);
        result = zero (targetDimension_);
        for (std::size_t i = 0; i < coeffs_.size (); ++i) {
          const auto & x = accessValueConstCast<T> (*this->dependency (1 + i));
          cwise (result) += (lambda * coeffs_[i]) * cwise (x);
        }
      }

      Dimension<T> targetDimension_;
      std::vector<double> coeffs_;
      int n_;
    };

    // Precompiled instantiations
    extern template class ConstantZero<double>;
    extern template class ConstantZero<Eigen::VectorXd>;
    extern template class ConstantZero<Eigen::RowVectorXd>;
    extern template class ConstantZero<Eigen::MatrixXd>;

    extern template class ConstantOne<double>;
    extern template class ConstantOne<Eigen::VectorXd>;
    extern template class ConstantOne<Eigen::RowVectorXd>;
    extern template class ConstantOne<Eigen::MatrixXd>;

    extern template class NumericConstant<double>;
    extern template class NumericConstant<Eigen::VectorXd>;
    extern template class NumericConstant<Eigen::RowVectorXd>;
    extern template class NumericConstant<Eigen::MatrixXd>;

    extern template class NumericMutable<double>;
    extern template class NumericMutable<Eigen::VectorXd>;
    extern template class NumericMutable<Eigen::RowVectorXd>;
    extern template class NumericMutable<Eigen::MatrixXd>;

    extern template class Convert<double, double>;
    extern template class Convert<Eigen::VectorXd, Eigen::VectorXd>;
    extern template class Convert<Eigen::RowVectorXd, Eigen::RowVectorXd>;
    extern template class Convert<Eigen::MatrixXd, Eigen::MatrixXd>;
    extern template class Convert<Eigen::VectorXd, double>;
    extern template class Convert<Eigen::RowVectorXd, double>;
    extern template class Convert<Eigen::MatrixXd, double>;
    extern template class Convert<Eigen::MatrixXd, Transposed<Eigen::MatrixXd>>;
    extern template class Convert<Eigen::RowVectorXd, Transposed<Eigen::VectorXd>>;
    extern template class Convert<Eigen::VectorXd, Transposed<Eigen::RowVectorXd>>;

    extern template class CWiseAdd<double, std::tuple<double, double>>;
    extern template class CWiseAdd<Eigen::VectorXd, std::tuple<Eigen::VectorXd, Eigen::VectorXd>>;
    extern template class CWiseAdd<Eigen::RowVectorXd, std::tuple<Eigen::RowVectorXd, Eigen::RowVectorXd>>;
    extern template class CWiseAdd<Eigen::MatrixXd, std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>>;

    extern template class CWiseAdd<double, ReductionOf<double>>;
    extern template class CWiseAdd<Eigen::VectorXd, ReductionOf<Eigen::VectorXd>>;
    extern template class CWiseAdd<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>>;
    extern template class CWiseAdd<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

    extern template class CWiseMul<double, std::tuple<double, double>>;
    extern template class CWiseMul<Eigen::VectorXd, std::tuple<Eigen::VectorXd, Eigen::VectorXd>>;
    extern template class CWiseMul<Eigen::RowVectorXd, std::tuple<Eigen::RowVectorXd, Eigen::RowVectorXd>>;
    extern template class CWiseMul<Eigen::MatrixXd, std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>>;
    extern template class CWiseMul<Eigen::VectorXd, std::tuple<double, Eigen::VectorXd>>;
    extern template class CWiseMul<Eigen::RowVectorXd, std::tuple<double, Eigen::RowVectorXd>>;
    extern template class CWiseMul<Eigen::MatrixXd, std::tuple<double, Eigen::MatrixXd>>;

    extern template class CWiseMul<double, ReductionOf<double>>;
    extern template class CWiseMul<Eigen::VectorXd, ReductionOf<Eigen::VectorXd>>;
    extern template class CWiseMul<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>>;
    extern template class CWiseMul<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

    extern template class CWiseNegate<double>;
    extern template class CWiseNegate<Eigen::VectorXd>;
    extern template class CWiseNegate<Eigen::RowVectorXd>;
    extern template class CWiseNegate<Eigen::MatrixXd>;

    extern template class CWiseInverse<double>;
    extern template class CWiseInverse<Eigen::VectorXd>;
    extern template class CWiseInverse<Eigen::RowVectorXd>;
    extern template class CWiseInverse<Eigen::MatrixXd>;

    extern template class CWiseConstantPow<double>;
    extern template class CWiseConstantPow<Eigen::VectorXd>;
    extern template class CWiseConstantPow<Eigen::RowVectorXd>;
    extern template class CWiseConstantPow<Eigen::MatrixXd>;

    extern template class ScalarProduct<Eigen::VectorXd, Eigen::VectorXd>;
    extern template class ScalarProduct<Eigen::RowVectorXd, Eigen::RowVectorXd>;

    extern template class SumOfLogarithms<Eigen::VectorXd>;
    extern template class SumOfLogarithms<Eigen::RowVectorXd>;

    extern template class MatrixProduct<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>;
    extern template class MatrixProduct<Eigen::RowVectorXd, Eigen::RowVectorXd, Eigen::MatrixXd>;
    extern template class MatrixProduct<Eigen::MatrixXd, Transposed<Eigen::MatrixXd>, Eigen::MatrixXd>;

    extern template class ShiftDelta<double>;
    extern template class ShiftDelta<Eigen::VectorXd>;
    extern template class ShiftDelta<Eigen::RowVectorXd>;
    extern template class ShiftDelta<Eigen::MatrixXd>;

    extern template class CombineDeltaShifted<double>;
    extern template class CombineDeltaShifted<Eigen::VectorXd>;
    extern template class CombineDeltaShifted<Eigen::RowVectorXd>;
    extern template class CombineDeltaShifted<Eigen::MatrixXd>;

    /*****************************************************************************
     * Numerical derivation helpers.
     */
    enum class NumericalDerivativeType { Disabled, ThreePoints, FivePoints };

    /// Configuration for a numerical derivation: what delta to use, and type of derivation.
    struct NumericalDerivativeConfiguration {
      NumericalDerivativeType type{NumericalDerivativeType::Disabled};
      ValueRef<double> delta{};
    };

    /** Helper used to generate data flow expressions computing a numerical derivative.
     *
     * For an expression e = f(x0,...,xn), which may be composed of multiple nodes.
     * dep is one of the expression dependencies: exists i, dep is node xi.
     * buildNodeWithDep(new_xi) should create the f(x0,...,new_xi,...,xn) node.
     * It must duplicate the expression by replacing dep by new_xi.
     * e must be a Value<NodeT>, and xi a Value<DepT>.
     *
     * This function generates a new node de_ddep.
     * de_ddep = g(delta,x0,...xn) = df/ddep value at x0,...,xn.
     * delta is the shift between points of computation of f (numerical derivation).
     * The pattern of computation points and delta are given by config.
     * After creation, the pattern and node for delta cannot change, but delta value can.
     */
    template <typename NodeT, typename DepT, typename B>
    ValueRef<NodeT> generateNumericalDerivative (Context & c, const NumericalDerivativeConfiguration & config,
                                                 NodeRef dep, const Dimension<DepT> & depDim,
                                                 const Dimension<NodeT> & nodeDim, B buildNodeWithDep) {
      if (config.delta == nullptr) {
        failureNumericalDerivationNotConfigured ();
      }
      switch (config.type) {
      case NumericalDerivativeType::ThreePoints: {
        // Shift {-1, +1}, coeffs {-0.5, +0.5}
        auto shift_m1 = ShiftDelta<DepT>::create (c, {config.delta, dep}, -1, depDim);
        auto shift_p1 = ShiftDelta<DepT>::create (c, {config.delta, dep}, 1, depDim);
        NodeRefVec combineDeps (3);
        combineDeps[0] = config.delta;
        combineDeps[1] = buildNodeWithDep (std::move (shift_m1));
        combineDeps[2] = buildNodeWithDep (std::move (shift_p1));
        return CombineDeltaShifted<NodeT>::create (c, std::move (combineDeps), 1, {-0.5, 0.5}, nodeDim);
      } break;
      case NumericalDerivativeType::FivePoints: {
        // Shift {-2, -1, +1, +2}, coeffs {1/12, -2/3, 2/3, -1/12}
        auto shift_m2 = ShiftDelta<DepT>::create (c, {config.delta, dep}, -2, depDim);
        auto shift_m1 = ShiftDelta<DepT>::create (c, {config.delta, dep}, -1, depDim);
        auto shift_p1 = ShiftDelta<DepT>::create (c, {config.delta, dep}, 1, depDim);
        auto shift_p2 = ShiftDelta<DepT>::create (c, {config.delta, dep}, 2, depDim);
        NodeRefVec combineDeps (5);
        combineDeps[0] = config.delta;
        combineDeps[1] = buildNodeWithDep (std::move (shift_m2));
        combineDeps[2] = buildNodeWithDep (std::move (shift_m1));
        combineDeps[3] = buildNodeWithDep (std::move (shift_p1));
        combineDeps[4] = buildNodeWithDep (std::move (shift_p2));
        return CombineDeltaShifted<NodeT>::create (c, std::move (combineDeps), 1,
                                                   {1. / 12., -2. / 3., 2. / 3., -1. / 12.}, nodeDim);
      } break;
      default:
        failureNumericalDerivationNotConfigured ();
      }
    }
  } // namespace dataflow
} // namespace bpp

#endif // BPP_NEWPHYL_DATAFLOWNUMERIC_H
