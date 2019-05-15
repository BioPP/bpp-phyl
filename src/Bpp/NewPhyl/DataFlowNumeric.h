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
#include <Bpp/Numeric/Parameter.h>
  
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

  /** @brief Store a dimension for type T.
   *
   * Declared but undefined by default.
   * Specialisations should be defined in the same header declaring the T type.
   * Specialisations should define a constructor from const T & : get the dimension of a T object.
   * If used in dataflow numeric nodes, it should be comparable and hashable.
   */
  template <typename T> struct Dimension;

  /// Specialisation of Dimension<T> for floating point types.
  template <> struct Dimension<double> : NoDimension {
    Dimension () = default;
    Dimension (const double &) {}
  };
  template <> struct Dimension<float> : NoDimension {
    Dimension () = default;
    Dimension (const float &) {}
  };

  template <> struct Dimension<Parameter> : NoDimension {
    Dimension () = default;
    Dimension (const Parameter &) {}
  };

  /** @brief Specialisation of Dimension<T> for eigen matrix types.
   *
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

    template <typename T = void>
    Parameter& zero (const Dimension<Parameter> &) {
      return *std::make_shared<Parameter>("Zero",0);
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
    
    template <typename T = void>
    Parameter& one (const Dimension<Parameter> &) {
      return *std::make_shared<Parameter>("One",1);
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
    
    inline Eigen::RowVectorXd convert (const Eigen::RowVectorXi& from, const Dimension<Eigen::RowVectorXd>& dim)
    {
      return from.template cast<double>();
    }

    inline Eigen::RowVectorXd convert (Eigen::RowVectorXi& from, const Dimension<Eigen::RowVectorXd>& dim)
    {
      return from.template cast<double>();
    }

    template <typename T, int Rows, int Cols, typename DerivedF>
    const DerivedF & convert (const Eigen::MatrixBase<DerivedF> & from,
                              const Dimension<Eigen::Matrix<T, Rows, Cols>> & dim,
                              typename std::enable_if<! std::is_same<DerivedF, Eigen::RowVectorXi>::value>::type* = 0) {
      return from.derived (); // matrix -> matrix, conversion will be done in the assignment
    }
    
    template <typename R, typename F> void convert (R & r, const F & from, const Dimension<R> & dim) {
      r = convert (from, dim);
      assert (Dimension<R> (r) == dim); // debug post check of size
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
    template <typename T>
    std::string debug (const T & t, typename std::enable_if<std::is_arithmetic<T>::value>::type* = 0) {
      // For basic arithmetic scalar types, just print the value itself
      using std::to_string;
      return "value=" + to_string (t);
    }
    
    template <typename T>
    std::string debug (const Parameter& t){
      // For basic arithmetic scalar types, just print the value itself
      using std::to_string;
      return "value=" + to_string(t.getValue());
    }

    template <typename Derived>
    std::string debug (const Eigen::MatrixBase<Derived> & m){//, typename std::enable_if<!std::is_same<Derived, Parameter const&>::value>::type* = 0) {
      // With matrices, check some numeric properties and encode results as text
      using std::to_string;
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

    // Utilities
    template <typename Predicate> void removeDependenciesIf (NodeRefVec & deps, Predicate p) {
      auto new_end = std::remove_if (deps.begin (), deps.end (), std::move (p));
      deps.erase (new_end, deps.end ()); // Truncate vector storage
    }

    /** @brief Template struct used to describe a dependency transformation before compute().
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

    /** @brief r = 0 for each component.
     * - r: T.
     *
     * Node construction should be done with the create static method.
     * Value is only created at first use (lazy).
     */
    template <typename T> class ConstantZero : public Value<T> {
    public:
      using Self = ConstantZero;

      /// Build a new ConstantZero node of the given dimension.
      static std::shared_ptr<Self> create (Context & c, const Dimension<T> & dim) {
        return cachedAs<Self> (c, std::make_shared<Self> (dim));
      }

      explicit ConstantZero (const Dimension<T> & dim) : Value<T> (NodeRefVec{}), targetDimension_ (dim) {
      }

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

    /** @brief r = 1 for each component.
     * - r: T.
     *
     * Node construction should be done with the create static method.
     * Value is only created at first use (lazy).
     */
    template <typename T> class ConstantOne : public Value<T> {
    public:
      using Self = ConstantOne;

      /// Build a new ConstantOne node of the given dimension.
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

    /** @brief r = constant_value.
     * - r: T.
     *
     * Node construction should be done with the create static method.
     * Value is set at construction, and cannot change.
     * Supports derivation.
     */
    template <typename T> class NumericConstant : public Value<T> {
    public:
      using Self = NumericConstant;

      /// Build a new NumericConstant node with T(args...) value.
      template <typename... Args> static std::shared_ptr<Self> create (Context & c, Args &&... args) {
        return cachedAs<Self> (c, std::make_shared<Self> (std::forward<Args> (args)...));
      }

      template <typename... Args>
      explicit NumericConstant (Args &&... args) : Value<T> (NodeRefVec{}, std::forward<Args> (args)...) {
        this->makeValid (); // Always valid
      }

      std::string debugInfo () const override {
        using namespace numeric;
        return debug (this->accessValueConst ());
      }

      std::string description() const override
      {
        using namespace numeric;
        return Node::description() + "\n"+ debug (this->accessValueConst ());
      }

      std::string color () const override {
        return "grey";
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

    /** @brief r = variable_value.
     * - r: T.
     *
     * Value is set at construction, and can be changed (will invalidate all dependent values).
     * Node construction should be done with the create static method.
     * Supports derivation.
     * This node has no Context merging support: mutable nodes are always different.
     */
    template <typename T> class NumericMutable : public Value<T> {
    public:
      using Self = NumericMutable;

      /// Build a new NumericMutable node with T(args...) value.
      template <typename... Args> static std::shared_ptr<Self> create (Context &, Args &&... args) {
        return std::make_shared<Self> (std::forward<Args> (args)...);
      }

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

      std::string description() const override
      {
        using namespace numeric;
        return Node::description() + "\n"+ debug (this->accessValueConst ());
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

    /** @brief r = convert(f).
     * - r: R.
     * - f: F, allows NumericalDependencyTransform.
     *
     * Convert from F to R type, semantics of numeric::convert.
     * Node construction should be done with the create static method.
     */
    template <typename R, typename F> class Convert : public Value<R> {
    public:
      using Self = Convert;
      using DepF = typename NumericalDependencyTransform<F>::DepType;

      /// Build a new Convert node with the given output dimensions.
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


    
  } // namespace dataflow
} // namespace bpp

#endif // BPP_NEWPHYL_DATAFLOWNUMERIC_H
