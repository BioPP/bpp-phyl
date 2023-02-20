//
// File: DataFlowCWiseComputing.h
// Authors:
//   Francois Gindraud (2017), Laurent GuÃÂ©guen (2019)
// Created: 2018-06-07 00:00:00
// Last modified: 2018-07-11 00:00:00
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_DATAFLOWCWISECOMPUTING_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_DATAFLOWCWISECOMPUTING_H

#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Phyl/Likelihood/DataFlow/ExtendedFloat.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Parameter.h>
#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <list>
#include <string>
#include <tuple>
#include <type_traits>

#include "DataFlowCWise.h"
#include "DataFlowNumeric.h"
#include "Definitions.h"

namespace bpp
{
/******************************************************************************
 * Data flow nodes for those numerical functions.
 *
 * TODO numerical simplification:
 * add(x,x) -> 2*x ? (and similar for mul, ...)
 * all deps constant => return constant ?
 */

template<typename Result, typename From> class CWiseAdd;
template<typename Result, typename From, typename Prop> class CWiseMean;
template<typename Result, typename From> class CWiseSub;
template<typename Result, typename From> class CWiseMul;
template<typename Result, typename From> class CWiseDiv;

template<typename R, typename T0, typename T1> class MatrixProduct;
// Return same type as given
template<typename T> class CWiseNegate;
template<typename T> class CWiseInverse;
template<typename T> class CWiseLog;
template<typename T> class CWiseExp;
template<typename T> class CWiseConstantPow;

// Return DataLik (R)
template<typename R, typename T0, typename T1> class ScalarProduct;
template<typename R, typename T0, typename T1> class LogSumExp;

// Return double
template<typename F> class SumOfLogarithms;

// Return same type as given

template<typename T> class ShiftDelta;
template<typename T> class CombineDeltaShifted;


/*************************************************************************
 * @brief r = f(x0) for each component or column
 * - r: R.
 * - x0: R
 * - f : function between columns or components of R (type F)
 *
 */

template<typename R, typename T, typename F> class CWiseApply : public Value<R>
{
public:
  using Self = CWiseApply;

  /// Build a new CWiseApply node.
  static ValueRef<R> create (Context& c, NodeRefVec&& deps, const Dimension<R>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 2);
    checkNthDependencyIsValue<T>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<F>(typeid (Self), deps, 1);

    return cachedAs<Value<R> >(c, std::make_shared<Self>(std::move (deps), dim));
  }

  CWiseApply (NodeRefVec&& deps, const Dimension<R>& dim)
    : Value<R>(std::move (deps)), bppLik_(size_t(dim.cols)), targetDimension_ (dim)
  {
    this->accessValueMutable().resize(dim.rows, dim.cols);
  }

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // CWiseApply additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  // df(X(theta))/dtheta|X(theta) = df/dtheta|X(theta) + df/dX.dX/dtheta|X(theta)
  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<R>::create (c, targetDimension_);
    }

    NodeRefVec dep(2);
    NodeRef df = this->dependency(1)->derive (c, node);

    if (df->hasNumericalProperty (NumericalProperty::ConstantZero))
      dep[0] = ConstantZero<R>::create (c, targetDimension_);
    else
      dep[0] = Self::create (c, {this->dependency(0), df}, targetDimension_);

    NodeRef dX = this->dependency(0)->derive (c, node);

    if (dX->hasNumericalProperty (NumericalProperty::ConstantZero))
      dep[1] = ConstantZero<R>::create (c, targetDimension_);
    else  // product df/dX.dX/dtheta|X(theta)
    {
      NodeRef dfX = this->dependency(1)->derive(c, NodeX);
      NodeRef dfXX = Self::create(c, {this->dependency(0), dfX}, targetDimension_);
      dep[1] = CWiseMul<R, std::tuple<VectorLik, R> >::create (c, {dX, dfXX}, targetDimension_);
    }

    return CWiseAdd<R, std::tuple<R, R> >::create (c, {std::move(dep)}, targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

  std::string shape() const override
  {
    return "doubleoctagon";
  }

  std::string color() const override
  {
    return "#5e5e5e";
  }

  std::string description() const override
  {
    return "Function Apply";
  }

private:
  void compute() override { compute<T>();}

  template<class U>
  typename std::enable_if<std::is_base_of<U, MatrixLik>::value && std::is_same<F, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<T>(*this->dependency (0));
    const auto& func = accessValueConstCast<F>(*this->dependency (1));

    for (auto i = 0; i < x0.cols(); i++)
      bppLik_[(size_t)i]=func(x0.col(i));

    copyBppToEigen(bppLik_, result);
    
#ifdef DEBUG
    std::cerr << "=== Function Apply === " << this << std::endl;
    std::cerr << "x0= " << x0 << std::endl;
    std::cerr << "res=" << result << std::endl;
    std::cerr << "=== end Function Apply === " << this << std::endl << std::endl;
#endif

  }

  /**
   *@brief For computation purpose
   *
   */
  
  std::vector<VectorLik> bppLik_;
  
  Dimension<R> targetDimension_;
};


/*************************************************************************
 * @brief r = x0 + x1 for each component.
 * - r: R.
 * - x0: T0.
 * - x1: T1.
 *
 * Values converted to R with the semantics of numeric::convert.
 * Node construction should be done with the create static method.
 *
 * Only defined for N = 2 for now.
 * The generic version is horrible in C++11 (lack of auto return).
 * Generic simplification routine is horrible too.
 */

template<typename R, typename T0, typename T1> class CWiseAdd<R, std::tuple<T0, T1> > : public Value<R>
{
public:
  using Self = CWiseAdd;

  /// Build a new CWiseAdd node with the given output dimensions.
  static ValueRef<R> create (Context& c, NodeRefVec&& deps, const Dimension<R>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 2);
    checkNthDependencyIsValue<T0>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<T1>(typeid (Self), deps, 1);
    // Select node implementation
    bool zeroDep0 = deps[0]->hasNumericalProperty (NumericalProperty::ConstantZero);
    bool zeroDep1 = deps[1]->hasNumericalProperty (NumericalProperty::ConstantZero);
    if (zeroDep0 && zeroDep1)
    {
      return ConstantZero<R>::create (c, dim);
    }
    else if (zeroDep0 && !zeroDep1)
    {
      return Convert<R, T1>::create (c, {deps[1]}, dim);
    }
    else if (!zeroDep0 && zeroDep1)
    {
      return Convert<R, T0>::create (c, {deps[0]}, dim);
    }
    else
    {
      return cachedAs<Value<R> >(c, std::make_shared<Self>(std::move (deps), dim));
    }
  }

  CWiseAdd (NodeRefVec&& deps, const Dimension<R>& dim)
    : Value<R>(std::move (deps)), targetDimension_ (dim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  std::string shape() const override
  {
    return "triangle";
  }

  std::string color() const override
  {
    return "#ff6060";
  }

  std::string description() const override
  {
    return "Add";
  }

  // Convert<T> additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<R>::create (c, targetDimension_);
    }
    constexpr std::size_t n = 2;
    NodeRefVec derivedDeps (n);
    for (std::size_t i = 0; i < n; ++i)
    {
      derivedDeps[i] = this->dependency (i)->derive (c, node);
    }
    return Self::create (c, std::move (derivedDeps), targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute() override { compute<R>();}

  template<class U>
  typename std::enable_if<!std::is_same<U, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<T0>(*this->dependency (0));
    const auto& x1 = accessValueConstCast<T1>(*this->dependency (1));
    result = cwise(x0) + cwise(x1);

#ifdef DEBUG
    std::cerr << "=== Add === " << this << std::endl;
    std::cerr << "x0= " << x0 << std::endl;
    std::cerr << "x1= " << x1 << std::endl;
    std::cerr << "res=" << result << std::endl;
    std::cerr << "=== end Add === " << this << std::endl << std::endl;
#endif
  }

  template<class U>
  typename std::enable_if<std::is_same<U, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<T0>(*this->dependency (0));
    const auto& x1 = accessValueConstCast<T1>(*this->dependency (1));

    result = [x0, x1](const VectorLik& x)->VectorLik {return x0(x) + x1(x);};
  }

  Dimension<R> targetDimension_;
};


/*************************************************************************
 * @brief r = x0 - x1 for each component.
 * - r: R.
 * - x0: T0.
 * - x1: T1.
 *
 * Values converted to R with the semantics of numeric::convert.
 * Node construction should be done with the create static method.
 *
 * Only defined for N = 2 for now.
 * The generic version is horrible in C++11 (lack of auto return).
 * Generic simplification routine is horrible too.
 */
template<typename R, typename T0, typename T1> class CWiseSub<R, std::tuple<T0, T1> > : public Value<R>
{
public:
  using Self = CWiseSub;

  /// Build a new CWiseSub node with the given output dimensions.
  static ValueRef<R> create (Context& c, NodeRefVec&& deps, const Dimension<R>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 2);
    checkNthDependencyIsValue<T0>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<T1>(typeid (Self), deps, 1);
    // Select node implementation
    bool zeroDep0 = deps[0]->hasNumericalProperty (NumericalProperty::ConstantZero);
    bool zeroDep1 = deps[1]->hasNumericalProperty (NumericalProperty::ConstantZero);
    if (zeroDep0 && zeroDep1)
    {
      return ConstantZero<R>::create (c, dim);
    }
    else if (zeroDep0 && !zeroDep1)
    {
      return Convert<R, T1>::create (c, {deps[1]}, dim);
    }
    else if (!zeroDep0 && zeroDep1)
    {
      return Convert<R, T0>::create (c, {deps[0]}, dim);
    }
    else
    {
      return cachedAs<Value<R> >(c, std::make_shared<Self>(std::move (deps), dim));
    }
  }

  CWiseSub (NodeRefVec&& deps, const Dimension<R>& dim)
    : Value<R>(std::move (deps)), targetDimension_ (dim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // Convert<T> additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<R>::create (c, targetDimension_);
    }
    constexpr std::size_t n = 2;
    NodeRefVec derivedDeps (n);
    for (std::size_t i = 0; i < n; ++i)
    {
      derivedDeps[i] = this->dependency (i)->derive (c, node);
    }
    return Self::create (c, std::move (derivedDeps), targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute () final
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<T0>(*this->dependency (0));
    const auto& x1 = accessValueConstCast<T1>(*this->dependency (1));
    result = cwise(x0) - cwise(x1);
#ifdef DEBUG
    std::cerr << "=== Sub === " << this << std::endl;
    std::cerr << "x0= " << x0 << std::endl;
    std::cerr << "x1= " << x1 << std::endl;
    std::cerr << "result= " << result << std::endl;
    std::cerr << "=== end Sub === " << this << std::endl << std::endl;
#endif
  }

  Dimension<R> targetDimension_;
};


/*************************************************************************
 * @brief r = sum (x_i), for each component.
 * - r: R.
 * - x_i: T.
 *
 * Sum of any number of T values into R.
 * Values converted to R with the semantics of numeric::convert.
 * Node construction should be done with the create static method.
 */
template<typename R, typename T> class CWiseAdd<R, ReductionOf<T> > : public Value<R>
{
public:
  using Self = CWiseAdd;

  /// Build a new CWiseAdd node with the given output dimensions.
  static ValueRef<R> create (Context& c, NodeRefVec&& deps, const Dimension<R>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyRangeIsValue<T>(typeid (Self), deps, 0, deps.size ());
    // Remove 0s from deps
    removeDependenciesIf (deps, [](const NodeRef& dep) -> bool {
        return dep->hasNumericalProperty (NumericalProperty::ConstantZero);
      });
    // Select node implementation
    if (deps.size () == 0)
    {
      return ConstantZero<R>::create (c, dim);
    }
    else if (deps.size () == 1)
    {
      return Convert<R, T>::create (c, std::move (deps), dim);
    }
    else if (deps.size () == 2)
    {
      return CWiseAdd<R, std::tuple<T, T> >::create (c, std::move (deps), dim);
    }
    else
    {
      return cachedAs<Value<R> >(c, std::make_shared<Self>(std::move (deps), dim));
    }
  }

  CWiseAdd (NodeRefVec&& deps, const Dimension<R>& dim)
    : Value<R>(std::move (deps)), targetDimension_ (dim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // CWiseAdd additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<R>::create (c, targetDimension_);
    }
    const auto n = this->nbDependencies ();
    NodeRefVec derivedDeps (n);
    for (std::size_t i = 0; i < n; ++i)
    {
      derivedDeps[i] = this->dependency (i)->derive (c, node);
    }
    return Self::create (c, std::move (derivedDeps), targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute() override { compute<T>();}

  template<class U>
  typename std::enable_if<!std::is_same<U, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    result = zero (targetDimension_);
    for (const auto& depNodeRef : this->dependencies ())
    {
      result += accessValueConstCast<T>(*depNodeRef);
    }
  }

  template<class U>
  typename std::enable_if<std::is_same<U, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    std::list<const T*> lT;
    for (const auto& depNodeRef : this->dependencies ())
    {
      lT.push_back(&accessValueConstCast<T>(*depNodeRef));
    }

    result = [lT, this](const VectorLik& x)->VectorLik {
               VectorLik r = zero (Dimension<VectorLik>(targetDimension_.cols, (Eigen::Index)1));

               for (const auto f:lT)
               {
                 cwise(r) += cwise((*f)(x));
               }
               return r;
             };
  }

  Dimension<R> targetDimension_;
};


/*************************************************************************
 * @brief r = sum (x_i), for each component either column-wise or row-wise, depending of entry & return types.
 *
 * - r: R.
 * - (x_i): T.
 *
 * Sum of any number of values of T into bits of R
 *.
 * Node construction should be done with the create static method.
 */

template<typename R, typename T> class CWiseAdd : public Value<R>
{
public:
  using Self = CWiseAdd;

  /// Build a new CWiseAdd node with the given output dimensions.
  static ValueRef<R> create (Context& c, NodeRefVec&& deps, const Dimension<R>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 1);

    if (deps[0]->hasNumericalProperty (NumericalProperty::ConstantZero))
      return ConstantZero<R>::create (c, dim);
    else
    {
      return cachedAs<Value<R> >(c, std::make_shared<Self>(std::move (deps), dim));
    }
  }

  CWiseAdd (NodeRefVec&& deps, const Dimension<R>& dim)
    : Value<R>(std::move (deps)), targetDimension_ (dim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // CWiseAdd additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<R>::create (c, targetDimension_);
    }
    return Self::create (c, {this->dependency(0)->derive (c, node)}, targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute() override { compute<R, T>();}

  template<class U, class V>
  typename std::enable_if<(std::is_base_of<V, MatrixLik>::value) && (std::is_same<U, RowLik>::value || std::is_same<U, Eigen::RowVectorXd>::value), void>::type
  compute ()
  {
    const auto& mat = accessValueConstCast<T>(*this->dependency(0));
    this->accessValueMutable () = mat.colwise().sum();
  }

  template<class U, class V>
  typename std::enable_if<(std::is_base_of<V, MatrixLik>::value) && (std::is_same<U, VectorLik>::value || std::is_same<U, Eigen::VectorXd>::value), void>::type
  compute ()
  {
    const auto& mat = accessValueConstCast<T>(*this->dependency(0));
    this->accessValueMutable () = mat.rowwise().sum();
  }

  template<class U, class V>
  typename std::enable_if<(std::is_base_of<V, VectorLik>::value) || (std::is_same<V, RowLik>::value || std::is_same<V, Eigen::RowVectorXd>::value), void>::type
  compute ()
  {
    const auto& vec = accessValueConstCast<T>(*this->dependency(0));
    this->accessValueMutable () = vec.sum();
  }


  Dimension<R> targetDimension_;
};


/*************************************************************************
 * @brief r = sum (p_i * x_i), for each component.
 * - r: R.
 * - x_i: T.
 * - p_i: P
 *
 * Sum of any number of T values multiplied per P values into R.
 * Values converted to R with the semantics of numeric::convert.
 * Node construction should be done with the create static method.
 */

template<typename R, typename T, typename P> class CWiseMean<R, ReductionOf<T>, ReductionOf<P> > : public Value<R>
{
public:
  using Self = CWiseMean;

  /// Build a new CWiseMean node with the given output dimensions.
  //  Last dependency is for P component
  static ValueRef<R> create (Context& c, NodeRefVec&& deps, const Dimension<R>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    if (deps.size() % 2 == 1)
      Exception ("CWiseMean needs an even number of dependencies, got " + std::to_string(deps.size()));

    size_t half = deps.size () / 2;

    checkDependencyRangeIsValue<T>(typeid (Self), deps, 0, half);
    checkDependencyRangeIsValue<P>(typeid (Self), deps, half, deps.size());

    // Remove both p_i and x_i from deps if one is null

    for (size_t i = 0; i < half; i++)
    {
      if (deps[i]->hasNumericalProperty (NumericalProperty::ConstantZero) ||
          deps[i + half]->hasNumericalProperty (NumericalProperty::ConstantZero))
      {
        deps[i] = 0;
        deps[i + half] = 0;
      }
    }

    removeDependenciesIf (deps, [](const NodeRef& dep)->bool {
        return dep == 0;
      });

    // if p_i * x_i are all Null
    if (deps.size() == 0)
      return ConstantZero<R>::create (c, dim);

    // Select node implementation
    return cachedAs<Value<R> >(c, std::make_shared<Self>(std::move (deps), dim));
  }

  CWiseMean (NodeRefVec&& deps, const Dimension<R>& dim)
    : Value<R>(std::move (deps)), targetDimension_ (dim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  std::string shape() const override
  {
    return "trapezium";
  }

  std::string color() const override
  {
    return "#ffd0d0";
  }

  std::string description() const override
  {
    return "Mean";
  }

  // CWiseMean additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<R>::create (c, targetDimension_);
    }
    const auto n = this->nbDependencies ();
    size_t half = n / 2;

    NodeRefVec derivedDeps_T (n);
    for (std::size_t i = 0; i < half; ++i)
    {
      derivedDeps_T[i] = this->dependency (i)->derive (c, node);
    }
    for (std::size_t i = half; i < n; ++i)
    {
      derivedDeps_T[i] = this->dependency (i);
    }
    NodeRef dR_dT = Self::create (c, std::move (derivedDeps_T), targetDimension_);

    NodeRefVec derivedDeps_P(n);
    for (std::size_t i = 0; i < half; ++i)
    {
      derivedDeps_P[i] = this->dependency (i);
    }
    for (std::size_t i = half; i < n; ++i)
    {
      derivedDeps_P[i] = this->dependency (i)->derive (c, node);
    }

    NodeRef dR_dP = Self::create (c, std::move (derivedDeps_P), targetDimension_);

    return CWiseAdd<R, std::tuple<R, R> >::create(c, {dR_dT, dR_dP}, targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute () final
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    result = zero (targetDimension_);
    size_t half = this->nbDependencies() / 2;
    for (size_t i = 0; i < half; i++)
    {
      cwise(result) += cwise (accessValueConstCast<P>(*this->dependency(i + half))) * cwise (accessValueConstCast<T>(*this->dependency(i)));
    }
  }

  Dimension<R> targetDimension_;
};

template<typename R, typename T, typename P> class CWiseMean<R, ReductionOf<T>, P> : public Value<R>
{
public:
  using Self = CWiseMean;

  /// Build a new CWiseMean node with the given output dimensions.
  //  Last dependency is for P component
  static ValueRef<R> create (Context& c, NodeRefVec&& deps, const Dimension<R>& dim)
  {
    // Check dependencies
    if (deps.size () <= 1)
      return ConstantZero<R>::create (c, dim);
    checkDependencyRangeIsValue<T>(typeid (Self), deps, 0, deps.size () - 1);
    checkNthDependencyIsValue<P>(typeid (Self), deps, deps.size () - 1);
    // if p_i are all Null
    if (deps[deps.size() - 1]->hasNumericalProperty (NumericalProperty::ConstantZero))
      return ConstantZero<R>::create (c, dim);
    // if x_i are all Null
    if (std::all_of (deps.begin (), deps.end () - 1, [](const NodeRef& dep)->bool {
        return dep->hasNumericalProperty (NumericalProperty::ConstantZero);
      }))
    {
      return ConstantZero<R>::create (c, dim);
    }

    // Select node implementation
    return cachedAs<Value<R> >(c, std::make_shared<Self>(std::move (deps), dim));
  }

  CWiseMean (NodeRefVec&& deps, const Dimension<R>& dim)
    : Value<R>(std::move (deps)), targetDimension_ (dim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  std::string shape() const override
  {
    return "trapezium";
  }

  std::string color() const override
  {
    return "#ffd0d0";
  }

  std::string description() const override
  {
    return "Mean";
  }

  // CWiseAdd additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<R>::create (c, targetDimension_);
    }
    const auto n = this->nbDependencies ();
    NodeRefVec derivedDeps_T (n);
    for (std::size_t i = 0; i < n - 1; ++i)
    {
      derivedDeps_T[i] = this->dependency (i)->derive (c, node);
    }
    derivedDeps_T[n - 1] = this->dependency (n - 1);
    NodeRef dR_dT = Self::create (c, std::move (derivedDeps_T), targetDimension_);

    NodeRefVec derivedDeps_P(n);
    for (std::size_t i = 0; i < n - 1; ++i)
    {
      derivedDeps_P[i] = this->dependency (i);
    }
    derivedDeps_P[n - 1] = this->dependency (n - 1)->derive (c, node);
    NodeRef dR_dP = Self::create (c, std::move (derivedDeps_P), targetDimension_);
    return CWiseAdd<R, std::tuple<R, R> >::create(c, {dR_dT, dR_dP}, targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute () final
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    auto& p = accessValueConstCast<P>(*this->dependency(this->nbDependencies() - 1));

#ifdef DEBUG
    std::cerr << "=== CWiseMean === " << this << std::endl;
    std::cerr << this->nbDependencies() << std::endl;
    std::cerr << this->dependency(this->nbDependencies() - 1) << std::endl;
    std::cerr << "p= " << p << std::endl;
    std::cerr << "=== end CWiseMean === " << this << std::endl << std::endl;
#endif

    result = cwise(p)[0] * cwise (accessValueConstCast<T>(*this->dependency(0)));
    for (Eigen::Index i = 1; i < Eigen::Index(this->nbDependencies() - 1); i++)
    {
      cwise(result) += cwise(p)[i] * cwise (accessValueConstCast<T>(*this->dependency(size_t(i))));
    }
  }

  Dimension<R> targetDimension_;
};


/*************************************************************************
 * @brief r = x0 * x1 for each component.
 * - r: R.
 * - x0: T0.
 * - x1: T1.
 *
 * Values converted to R with the semantics of numeric::convert.
 * Node construction should be done with the create static method.
 *
 * Only defined for N = 2 for now (same constraints as CWiseAdd for genericity).
 */

template<typename R, typename T0, typename T1> class CWiseMul<R, std::tuple<T0, T1> > : public Value<R>
{
public:
  using Self = CWiseMul;

  /// Build a new CWiseMul node with the given output dimensions.
  static ValueRef<R> create (Context& c, NodeRefVec&& deps, const Dimension<R>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 2);
    checkNthDependencyIsValue<T0>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<T1>(typeid (Self), deps, 1);
    // Return 0 if any 0.
    if (std::any_of (deps.begin (), deps.end (), [](const NodeRef& dep) -> bool {
        return dep->hasNumericalProperty (NumericalProperty::ConstantZero);
      }))
    {
      return ConstantZero<R>::create (c, dim);
    }
    // Select node implementation
    bool oneDep0 = deps[0]->hasNumericalProperty (NumericalProperty::ConstantOne);
    bool oneDep1 = deps[1]->hasNumericalProperty (NumericalProperty::ConstantOne);
    if (oneDep0 && oneDep1)
    {
      return ConstantOne<R>::create (c, dim);
    }
    else if (oneDep0 && !oneDep1)
    {
      return Convert<R, T1>::create (c, {deps[1]}, dim);
    }
    else if (!oneDep0 && oneDep1)
    {
      return Convert<R, T0>::create (c, {deps[0]}, dim);
    }
    else
    {
      return cachedAs<Value<R> >(c, std::make_shared<Self>(std::move (deps), dim));
    }
  }

  CWiseMul (NodeRefVec&& deps, const Dimension<R>& dim)
    : Value<R>(std::move (deps)), targetDimension_ (dim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  std::string shape() const override
  {
    return "house";
  }

  std::string color() const override
  {
    return "#ff9090";
  }

  std::string description() const override
  {
    return "Mul";
  }

  // CWiseMul additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<R>::create (c, targetDimension_);
    }
    constexpr std::size_t n = 2;
    NodeRefVec addDeps (n);
    for (std::size_t i = 0; i < n; ++i)
    {
      NodeRefVec ithMulDeps = this->dependencies ();
      ithMulDeps[i] = this->dependency (i)->derive (c, node);
      addDeps[i] = Self::create (c, std::move (ithMulDeps), targetDimension_);
    }
    return CWiseAdd<R, std::tuple<R, R> >::create (c, std::move (addDeps), targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute() override { compute<T0, T1>();}

  template<class U, class V>
  typename std::enable_if<!std::is_same<U, TransitionFunction>::value && !std::is_same<V, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<U>(*this->dependency (0));
    const auto& x1 = accessValueConstCast<V>(*this->dependency (1));
#ifdef DEBUG
    std::cerr << "=== Mul === " << this << std::endl;
    std::cerr << "x0= "     << x0 << std::endl;
    std::cerr << "x1= "     << x1 << std::endl;
#endif
    result = cwise (x0) * cwise (x1);

#ifdef DEBUG
    std::cerr << "result= " << result << std::endl;
    std::cerr << "=== end Mul === " << this << std::endl << std::endl;
#endif
  }

  template<class U, class V>
  typename std::enable_if<std::is_same<U, TransitionFunction>::value && std::is_same<V, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<U>(*this->dependency (0));
    const auto& x1 = accessValueConstCast<V>(*this->dependency (1));

    result = [x0, x1](const VectorLik& x)->VectorLik {return cwise(x0(x)) * cwise(x1(x));};

#ifdef DEBUG
    std::cerr << "=== Mul Transition Function X Transition Function === " << this << std::endl;
    std::cerr << "=== end Mul Transition Function X Transition Function  === " << this << std::endl << std::endl;
#endif
  }

  template<class U, class V>
  typename std::enable_if<std::is_same<U, TransitionFunction>::value && !std::is_same<V, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<U>(*this->dependency (0));
    const auto& x1 = accessValueConstCast<V>(*this->dependency (1));

    result = [x0, x1](const VectorLik& x)->VectorLik {return cwise(x0(x)) * cwise(x1);};

#ifdef DEBUG
    std::cerr << "=== Mul Transition Function X Normal === " << this << std::endl;
    std::cerr << "x1= "     << x1 << std::endl;
    std::cerr << "=== end Mul Transition Function  X Normal  === " << this << std::endl << std::endl;
#endif
}

  template<class U, class V>
  typename std::enable_if<!std::is_same<U, TransitionFunction>::value && std::is_same<V, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<U>(*this->dependency (0));
    const auto& x1 = accessValueConstCast<V>(*this->dependency (1));

    result = [x0, x1](const VectorLik& x)->VectorLik {return cwise(x1(x)) * cwise(x0);};

#ifdef DEBUG
    std::cerr << "=== Mul Normal X Transition Function === " << this << std::endl;
    std::cerr << "x0= "     << x0 << std::endl;
    std::cerr << "=== end Mul Normal X Transition Function  === " << this << std::endl << std::endl;
#endif
}

  Dimension<R> targetDimension_;
};

/*************************************************************************
 * @brief r = prod (x_i), for each component.
 * - r: R.
 * - x_i: T.
 *
 * Product of any number of T values into R.
 * Values converted to R with the semantics of numeric::convert.
 * Node construction should be done with the create static method.
 */
template<typename R, typename T> class CWiseMul<R, ReductionOf<T> > : public Value<R>
{
public:
  using Self = CWiseMul;

  /// Build a new CWiseMul node with the given output dimensions.
  static ValueRef<R> create (Context& c, NodeRefVec&& deps, const Dimension<R>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyRangeIsValue<T>(typeid (Self), deps, 0, deps.size ());
    // If there is a 0 return 0.
    if (std::any_of (deps.begin (), deps.end (), [](const NodeRef& dep) -> bool {
        return dep->hasNumericalProperty (NumericalProperty::ConstantZero);
      }))
    {
      return ConstantZero<R>::create (c, dim);
    }
    // Remove 1s from deps
    removeDependenciesIf (deps, [](const NodeRef& dep)->bool {
        return dep->hasNumericalProperty (NumericalProperty::ConstantOne);
      });
    // Select node implementation
    if (deps.size () == 0)
    {
      return ConstantOne<R>::create (c, dim);
    }
    else if (deps.size () == 1)
    {
      return Convert<R, T>::create (c, std::move (deps), dim);
    }
    else if (deps.size () == 2)
    {
      return CWiseMul<R, std::tuple<T, T> >::create (c, std::move (deps), dim);
    }
    else
    {
      return cachedAs<Value<R> >(c, std::make_shared<Self>(std::move (deps), dim));
    }
  }

  CWiseMul (NodeRefVec&& deps, const Dimension<R>& dim)
    : Value<R>(std::move (deps)), targetDimension_ (dim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // CWiseMul additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<R>::create (c, targetDimension_);
    }
    const auto n = this->nbDependencies ();
    NodeRefVec addDeps (n);
    for (std::size_t i = 0; i < n; ++i)
    {
      NodeRefVec ithMulDeps = this->dependencies ();
      ithMulDeps[i] = this->dependency (i)->derive (c, node);
      addDeps[i] = Self::create (c, std::move (ithMulDeps), targetDimension_);
    }
    return CWiseAdd<R, ReductionOf<R> >::create (c, std::move (addDeps), targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute () final
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    result = one (targetDimension_);
    for (const auto& depNodeRef : this->dependencies ())
    {
      cwise(result) *= cwise (accessValueConstCast<T>(*depNodeRef));
    }
  }

  Dimension<R> targetDimension_;
};



/*************************************************************************
 * @brief r = x0 / x1 for each component.
 * - r: R.
 * - x0: T0.
 * - x1: T1.
 *
 * Values converted to R with the semantics of numeric::convert.
 * Node construction should be done with the create static method.
 *
 * Only defined for N = 2 for now (same constraints as CWiseAdd for genericity).
 */
template<typename R, typename T0, typename T1> class CWiseDiv<R, std::tuple<T0, T1> > : public Value<R>
{
public:
  using Self = CWiseDiv;

  /// Build a new CWiseDiv node with the given output dimensions.
  static ValueRef<R> create (Context& c, NodeRefVec&& deps, const Dimension<R>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 2);
    checkNthDependencyIsValue<T0>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<T1>(typeid (Self), deps, 1);
    // Select node
    if (deps[1]->hasNumericalProperty (NumericalProperty::ConstantOne))
      return Convert<R, T0>::create (c, {deps[0]}, dim);
    else if (deps[0]->hasNumericalProperty (NumericalProperty::ConstantOne))
      return CWiseInverse<R>::create (c, {deps[1]}, dim);
    else
      return cachedAs<Value<R> >(c, std::make_shared<Self>(std::move (deps), dim));
  }

  CWiseDiv (NodeRefVec&& deps, const Dimension<R>& dim)
    : Value<R>(std::move (deps)), targetDimension_ (dim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  std::string shape() const override
  {
    return "house";
  }

  std::string color() const override
  {
    return "#ff9090";
  }

  std::string description() const override
  {
    return "Div";
  }

  // CWiseDiv additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<R>::create (c, targetDimension_);
    }

    auto f0 = this->dependency (0);
    auto f1 = this->dependency (1);
    auto fp0 = this->dependency (0)->derive (c, node);
    auto fp1 = this->dependency (1)->derive (c, node);

    auto upv = CWiseMul<R, std::tuple<T0, T1> >::create(c, {fp0, f1}, targetDimension_);
    auto vpu = CWiseMul<R, std::tuple<T0, T1> >::create(c, {fp1, f0}, targetDimension_);
    auto diff = CWiseSub<R, std::tuple<R, R> >::create(c, {upv, vpu}, targetDimension_);

    return CWiseMul<R, std::tuple<R, R> >::create (
      c, {CWiseConstantPow<R>::create(c, {f1}, -2., 1., targetDimension_), diff}, targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute() override { compute<T0, T1>();}

  template<class U, class V>
  typename std::enable_if<!std::is_same<U, TransitionFunction>::value && !std::is_same<V, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<U>(*this->dependency (0));
    const auto& x1 = accessValueConstCast<V>(*this->dependency (1));
    result = cwise (x0) / cwise (x1);
  }

  template<class U, class V>
  typename std::enable_if<std::is_same<U, TransitionFunction>::value && std::is_same<V, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<U>(*this->dependency (0));
    const auto& x1 = accessValueConstCast<V>(*this->dependency (1));

    result = [x0, x1](const VectorLik& x)->VectorLik {return cwise(x0(x)) / cwise(x1(x));};
  }

  template<class U, class V>
  typename std::enable_if<std::is_same<U, TransitionFunction>::value && !std::is_same<V, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<U>(*this->dependency (0));
    const auto& x1 = accessValueConstCast<V>(*this->dependency (1));

    result = [x0, x1](const VectorLik& x)->VectorLik {return cwise(x0(x)) / cwise(x1);};
  }

  template<class U, class V>
  typename std::enable_if<!std::is_same<U, TransitionFunction>::value && std::is_same<V, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<U>(*this->dependency (0));
    const auto& x1 = accessValueConstCast<V>(*this->dependency (1));

    result = [x0, x1](const VectorLik& x)->VectorLik {return cwise(x1(x)) / cwise(x0);};
  }

  Dimension<R> targetDimension_;
};


/*************************************************************************
 * @brief r = -x, for each component.
 * - r, x: T.
 *
 * Node construction should be done with the create static method.
 */

template<typename T> class CWiseNegate : public Value<T>
{
public:
  using Self = CWiseNegate;

  /// Build a new CWiseNegate node with the given output dimensions.
  static ValueRef<T> create (Context& c, NodeRefVec&& deps, const Dimension<T>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 1);
    checkNthDependencyIsValue<T>(typeid (Self), deps, 0);
    // Select node
    if (deps[0]->hasNumericalProperty (NumericalProperty::ConstantZero))
    {
      return ConstantZero<T>::create (c, dim);
    }
    else
    {
      return cachedAs<Value<T> >(c, std::make_shared<Self>(std::move (deps), dim));
    }
  }

  CWiseNegate (NodeRefVec&& deps, const Dimension<T>& dim)
    : Value<T>(std::move (deps)), targetDimension_ (dim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // CWiseNegate additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<T>::create (c, targetDimension_);
    }
    return Self::create (c, {this->dependency (0)->derive (c, node)}, targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute () final
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x = accessValueConstCast<T>(*this->dependency (0));
    result = -x;
  }

  Dimension<T> targetDimension_;
};


/*************************************************************************
 * @brief r = 1/x for each component.
 * - r, x: T.
 *
 * Node construction should be done with the create static method.
 */

template<typename T> class CWiseInverse : public Value<T>
{
public:
  using Self = CWiseInverse;

  /// Build a new CWiseInverse node with the given output dimensions.
  static ValueRef<T> create (Context& c, NodeRefVec&& deps, const Dimension<T>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 1);
    checkNthDependencyIsValue<T>(typeid (Self), deps, 0);
    // Select node
    if (deps[0]->hasNumericalProperty (NumericalProperty::ConstantOne))
    {
      return ConstantOne<T>::create (c, dim);
    }
    else
    {
      return cachedAs<Value<T> >(c, std::make_shared<Self>(std::move (deps), dim));
    }
  }

  CWiseInverse (NodeRefVec&& deps, const Dimension<T>& dim)
    : Value<T>(std::move (deps)), targetDimension_ (dim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // CWiseInverse additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<T>::create (c, targetDimension_);
    }
    // -1/x^2 * x'
    const auto& dep = this->dependency (0);
    return CWiseMul<T, std::tuple<T, T> >::create (
      c, {CWiseConstantPow<T>::create (c, {dep}, -2., -1., targetDimension_), dep->derive (c, node)},
      targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute () final
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x = accessValueConstCast<T>(*this->dependency (0));
    result = inverse (cwise (x));
  }

  Dimension<T> targetDimension_;
};


/*************************************************************************
 * @brief r = log(x) for each component.
 * - r, x: T.
 *
 * Node construction should be done with the create static method.
 */

using std::log;

template<typename T> class CWiseLog : public Value<T>
{
public:
  using Self = CWiseLog;

  /// Build a new CWiseLog node with the given output dimensions.
  static ValueRef<T> create (Context& c, NodeRefVec&& deps, const Dimension<T>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 1);
    checkNthDependencyIsValue<T>(typeid (Self), deps, 0);
    // Select node
    if (deps[0]->hasNumericalProperty (NumericalProperty::ConstantOne))
    {
      return ConstantZero<T>::create (c, dim);
    }
    else
    {
      return cachedAs<Value<T> >(c, std::make_shared<Self>(std::move (deps), dim));
    }
  }

  CWiseLog (NodeRefVec&& deps, const Dimension<T>& dim)
    : Value<T>(std::move (deps)), targetDimension_ (dim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // CWiseLog additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<T>::create (c, targetDimension_);
    }
    // x'/x
    const auto& dep = this->dependency (0);
    return CWiseMul<T, std::tuple<T, T> >::create (
      c, {CWiseInverse<T>::create (c, {dep}, targetDimension_), dep->derive (c, node)}, targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute () final
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x = accessValueConstCast<T>(*this->dependency (0));
    result = log (cwise (x));
  }

  Dimension<T> targetDimension_;
};


/*************************************************************************
 * @brief r = exp(x) for each component.
 * - r, x: T.
 *
 * Node construction should be done with the create static method.
 */

using std::exp;

template<typename T> class CWiseExp : public Value<T>
{
public:
  using Self = CWiseExp;

  /// Build a new CWiseExp node with the given output dimensions.
  static ValueRef<T> create (Context& c, NodeRefVec&& deps, const Dimension<T>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 1);
    checkNthDependencyIsValue<T>(typeid (Self), deps, 0);
    // Select node
    if (deps[0]->hasNumericalProperty (NumericalProperty::ConstantZero))
    {
      return ConstantOne<T>::create (c, dim);
    }
    else
    {
      return cachedAs<Value<T> >(c, std::make_shared<Self>(std::move (deps), dim));
    }
  }

  CWiseExp (NodeRefVec&& deps, const Dimension<T>& dim)
    : Value<T>(std::move (deps)), targetDimension_ (dim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // CWiseExp additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<T>::create (c, targetDimension_);
    }
    // x'* exp(x)
    const auto& dep = this->dependency (0);
    return CWiseMul<T, std::tuple<T, T> >::create (
      c, {this->shared_from_this(), dep->derive (c, node)}, targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute () final
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x = accessValueConstCast<T>(*this->dependency (0));
    result = exp (cwise (x));
  }

  Dimension<T> targetDimension_;
};


/*************************************************************************
 * @brief r = factor * pow (x, exponent) for each component.
 * - r, x: T.
 * - exponent, factor: double (constant parameter of the node).
 *
 * Node construction should be done with the create static method.
 */

template<typename T> class CWiseConstantPow : public Value<T>
{
public:
  using Self = CWiseConstantPow;

  /// Build a new CWiseConstantPow node with the given output dimensions and factors.
  static ValueRef<T> create (Context& c, NodeRefVec&& deps, double exponent, double factor,
                             const Dimension<T>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 1);
    checkNthDependencyIsValue<T>(typeid (Self), deps, 0);
    // Select node implementation
    if (exponent == 0. || deps[0]->hasNumericalProperty (NumericalProperty::ConstantOne))
    {
      // pow (x, exponent) == 1
      using namespace numeric;
      return NumericConstant<T>::create (c, factor * one (dim));
    }
    else if (exponent == 1.)
    {
      // pow (x, exponent) == x
      return CWiseMul<T, std::tuple<double, T> >::create (
        c, {NumericConstant<double>::create (c, factor), deps[0]}, dim);
    }
    else if (exponent == -1.)
    {
      // pow (x, exponent) = 1/x
      return CWiseMul<T, std::tuple<double, T> >::create (
        c,
        {NumericConstant<double>::create (c, factor), CWiseInverse<T>::create (c, std::move (deps), dim)},
        dim);
    }
    else
    {
      return cachedAs<Value<T> >(c, std::make_shared<Self>(std::move (deps), exponent, factor, dim));
    }
  }

  CWiseConstantPow (NodeRefVec&& deps, double exponent, double factor, const Dimension<T>& dim)
    : Value<T>(std::move (deps)), targetDimension_ (dim), exponent_ (exponent), factor_ (factor) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_) +
           " exponent=" + std::to_string (exponent_) + " factor=" + std::to_string (factor_);
  }

  // CWiseConstantPow additional arguments = (exponent_, factor_).
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    const auto* derived = dynamic_cast<const Self*>(&other);
    return derived != nullptr && exponent_ == derived->exponent_ && factor_ == derived->factor_;
  }
  std::size_t hashAdditionalArguments () const final
  {
    std::size_t seed = 0;
    combineHash (seed, exponent_);
    combineHash (seed, factor_);
    return seed;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<T>::create (c, targetDimension_);
    }
    // factor * (exponent * x^(exponent - 1)) * x'
    const auto& dep = this->dependency (0);
    auto dpow = Self::create (c, {dep}, exponent_ - 1., factor_ * exponent_, targetDimension_);
    return CWiseMul<T, std::tuple<T, T> >::create (c, {dpow, dep->derive (c, node)}, targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), exponent_, factor_, targetDimension_);
  }

private:
  void compute () final
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& x = accessValueConstCast<T>(*this->dependency (0));
    result = factor_ * pow (cwise (x), exponent_);
  }

  Dimension<T> targetDimension_;
  double exponent_;
  double factor_;
};


/*************************************************************************
 * @brief r = x0 * x1 (dot product).
 * - r: double.
 * - x0: T0 (vector-like).
 * - x1: T1 (vector-like).
 *
 * Node construction should be done with the create static method.
 */
template<typename R, typename T0, typename T1> class ScalarProduct : public Value<R>
{
public:
  using Self = ScalarProduct;

  /// Build a new ScalarProduct node.
  static ValueRef<R> create (Context& c, NodeRefVec&& deps)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 2);
    checkNthDependencyIsValue<T0>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<T1>(typeid (Self), deps, 1);
    // Select node
    if (deps[0]->hasNumericalProperty (NumericalProperty::ConstantZero) ||
        deps[1]->hasNumericalProperty (NumericalProperty::ConstantZero))
    {
      return ConstantZero<R>::create (c, Dimension<R> ());
    }
    else
    {
      return cachedAs<Value<R> >(c, std::make_shared<Self>(std::move (deps)));
    }
  }

  ScalarProduct (NodeRefVec&& deps) : Value<R>(std::move (deps)) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ());
  }

  // ScalarProduct additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<R>::create (c, Dimension<R> ());
    }
    const auto& x0 = this->dependency (0);
    const auto& x1 = this->dependency (1);
    auto dx0_prod = Self::create (c, {x0->derive (c, node), x1});
    auto dx1_prod = Self::create (c, {x0, x1->derive (c, node)});
    return CWiseAdd<R, std::tuple<R, R> >::create (c, {dx0_prod, dx1_prod},
                                                   Dimension<R> ());
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps));
  }

private:
  void compute () final
  {
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<T0>(*this->dependency (0));
    const auto& x1 = accessValueConstCast<T1>(*this->dependency (1));
    auto d = x0.dot (x1); // Using lhs.dot(rhs) method from Eigen only
    result = d;
  }
};


/*************************************************************************
 * @brief r = sum_{v in m} log (v).
 * - r: double.
 * - m: F (matrix-like type).
 *
 * or  r = sum_{i in 0:len(m)-1} log (p[i]*m[i]).
 * - r: DataLik
 * - m: F (matrix-like type).
 * - p: <Eigen::RowVectorXd>
 *
 * The node has no dimension (double).
 * The dimension of m should be provided for derivation.
 * Node construction should be done with the create static method.
 */

template<typename F> class SumOfLogarithms : public Value<DataLik>
{
public:
  using Self = SumOfLogarithms;

  /// Build a new SumOfLogarithms node with the given input matrix dimensions.
  static ValueRef<DataLik> create (Context& c, NodeRefVec&& deps, const Dimension<F>& mDim)
  {
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorMinSize (typeid (Self), deps, 1);
    checkNthDependencyIsValue<F>(typeid (Self), deps, 0);
    if (deps.size() == 2)
      checkNthDependencyIsValue<Eigen::RowVectorXi>(typeid (Self), deps, 1);
    return cachedAs<Value<DataLik> >(c, std::make_shared<Self>(std::move (deps), mDim));
  }

  SumOfLogarithms (NodeRefVec&& deps, const Dimension<F>& mDim)
    : Value<DataLik>(std::move (deps)), mTargetDimension_ (mDim) // , temp_() {
//         if (dependencies().size()==2)
//         {
//           const auto & p = accessValueConstCast<Eigen::VectorXi> (*this->dependency (1));
//          temp_.resize(p.size());
//         }
  {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ());
  }

  // SumOfLogarithms additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    const auto& m = this->dependency (0);
    auto dm_dn = m->derive (c, node);
    auto m_inverse = CWiseInverse<F>::create (c, {m}, mTargetDimension_);
    if (nbDependencies() == 2 && dependency(1))
    {
      auto m2 = CWiseMul<F, std::tuple<F, Eigen::RowVectorXi> >::create (c, {m_inverse, dependency (1)}, mTargetDimension_);
      return ScalarProduct<typename F::Scalar, F, F>::create (c, {std::move (dm_dn), std::move (m2)});
    }
    else
      return ScalarProduct<typename F::Scalar, F, F>::create (c, {std::move (dm_dn), std::move (m_inverse)});
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), mTargetDimension_);
  }

private:
  void compute() override { compute<F>();}

  template<class G = F>
  typename std::enable_if<std::is_convertible<G*, Eigen::DenseBase<G>*>::value, void>::type
  compute ()
  {
#ifdef DEBUG
    std::cerr << "=== SumOfLogarithms === " << this << std::endl;
#endif

    auto& result = this->accessValueMutable ();

    const auto& m = accessValueConstCast<F>(*this->dependency (0));

    if (nbDependencies() == 1)
    {
      const ExtendedFloat product = m.unaryExpr ([](double d) {
          ExtendedFloat ef{d};
          ef.normalize_small ();
          return ef;
        }).redux ([](const ExtendedFloat& lhs, const ExtendedFloat& rhs) {
          auto r = ExtendedFloat::denorm_mul (lhs, rhs);
          r.normalize_small ();
          return r;
        });
      result = product.log();
#ifdef DEBUG
      std::cerr << "product= " << product << std::endl;
      std::cerr << "result log= " << result << std::endl;
#endif
    }
    else
    {
      const auto& p = accessValueConstCast<Eigen::RowVectorXi>(*this->dependency (1));
      // Old version:
      // double resold  = (numeric::cwise(m).log() * numeric::cwise(p)).sum();

      temp_ = m.unaryExpr ([](double d) {
          ExtendedFloat ef{d};
          ef.normalize ();
          return ef;
        });

      for (Eigen::Index i = 0; i < Eigen::Index(p.size()); i++)
      {
        temp_[i] = temp_[i].pow(p[i]);
      }

      const ExtendedFloat product = temp_.redux ([](const ExtendedFloat& lhs, const ExtendedFloat& rhs) {
          auto r = ExtendedFloat::denorm_mul (lhs, rhs);
          r.normalize ();
          return r;
        });

      result = product.log ();
#ifdef DEBUG
      std::cerr << "PRODUCT= " << product << std::endl;
      std::cerr << "RESULT log= " << result << std::endl;
#endif
    }
#ifdef DEBUG
    std::cerr << "=== end SumOfLogarithms === " << this << std::endl;
#endif
  }

  template<class G = F>
  typename std::enable_if<std::is_convertible<G*, ExtendedFloatEigenBase<G>*>::value, void>::type
  compute ()
  {
#ifdef DEBUG
    std::cerr << "=== SumOfLogarithms === " << this << std::endl;
#endif

    auto& result = this->accessValueMutable ();

    const auto& m = accessValueConstCast<F>(*this->dependency (0));

    if (nbDependencies() == 1)
    {
      const ExtendedFloat product = m.float_part().unaryExpr ([](double d) {
          ExtendedFloat ef{d};
          ef.normalize ();
          return ef;
        }).redux ([](const ExtendedFloat& lhs, const ExtendedFloat& rhs) {
          auto r = ExtendedFloat::denorm_mul (lhs, rhs);
          r.normalize_small ();
          return r;
        });
      result = product.log() + double(m.exponent_part() * m.size()) * ExtendedFloat::ln_radix;
#ifdef DEBUG
      std::cerr << "product= " << product << "* 2^" << m.size() * m.exponent_part() << std::endl;
      std::cerr << "result log= " << result << std::endl;
#endif
    }
    else
    {
      const auto& p = accessValueConstCast<Eigen::RowVectorXi>(*this->dependency (1));
      // Old version:
      // double resold  = (numeric::cwise(m).log() * numeric::cwise(p)).sum();

      temp_ = m.float_part().unaryExpr ([](double d) {
          ExtendedFloat ef{d};
          ef.normalize ();
          return ef;
        });

      for (Eigen::Index i = 0; i < Eigen::Index(p.size()); i++)
      {
        temp_[i] = temp_[i].pow(p[i]);
      }

      const ExtendedFloat product = temp_.redux ([](const ExtendedFloat& lhs, const ExtendedFloat& rhs) {
          auto r = ExtendedFloat::denorm_mul (lhs, rhs);
          r.normalize ();
          return r;
        });

      result = product.log () + double(m.exponent_part() * p.sum()) * ExtendedFloat::ln_radix;
#ifdef DEBUG
      std::cerr << "RESULT log= " << result << std::endl;
#endif
    }
#ifdef DEBUG
    std::cerr << "=== end SumOfLogarithms === " << this << std::endl;
#endif
  }

  Dimension<F> mTargetDimension_;

  Eigen::Matrix<ExtendedFloat, 1, Eigen::Dynamic> temp_;
};


/*************************************************************************
 * @brief r = log(sum_i p_i * exp (v_i))
 * - r: R.
 * - v: T0 (vector like type).
 * - p: T1 (vector like type).
 *
 * The node has no dimension (double).
 * Node construction should be done with the create static method.
 */

template<typename R, typename T0, typename T1> class LogSumExp : public Value<R>
{
public:
  using Self = LogSumExp;

  /// Build a new LogSumExp node with the given input matrix dimensions.
  static ValueRef<R> create (Context& c, NodeRefVec&& deps, const Dimension<T0>& mDim)
  {
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 2);
    checkNthDependencyIsValue<T0>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<T1>(typeid (Self), deps, 1);
    // Select node
    return cachedAs<Value<R> >(c, std::make_shared<Self>(std::move (deps), mDim));
  }

  LogSumExp (NodeRefVec&& deps, const Dimension<T0>& mDim)
    : Value<R>(std::move (deps)), mTargetDimension_ (mDim) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ());
  }

  // LogSumExp additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<R>::create (c, Dimension<R>());
    }
    const auto& v = this->dependency (0);
    const auto& p = this->dependency (1);
    auto diffvL = CWiseSub<T0, std::tuple<R, T0> >::create(c, {this->shared_from_this(), v}, mTargetDimension_);
    auto expdiffvL = CWiseExp<T0>::create(c, {std::move(diffvL)}, mTargetDimension_);
    auto dp_prod = ScalarProduct<R, T0, T1>::create (c, {expdiffvL, p->derive (c, node)});
    auto pexpdiffvL = CWiseMul<T0, std::tuple<T0, T1> >::create(c, {expdiffvL, p}, mTargetDimension_);
    auto dv_prod = ScalarProduct<R, T0, T1>::create (c, {std::move(pexpdiffvL), v->derive (c, node)});
    return CWiseAdd<R, std::tuple<R, R> >::create (c, {std::move (dv_prod), std::move (dp_prod)}, Dimension<R>());
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), mTargetDimension_);
  }

private:
  void compute () final
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& v = accessValueConstCast<T0>(*this->dependency (0));
    const auto& p = accessValueConstCast<T1>(*this->dependency (1));

    auto M = v.maxCoeff();
    if (isinf(M))
      result = M;
    else
    {
      T0 v2;
      v2 = exp(cwise(v - T0::Constant(mTargetDimension_.rows, mTargetDimension_.cols, M)));
      auto ve = v2.dot(p);
      result = log(ve) + M;
    }
  }

  Dimension<T0> mTargetDimension_;
};


/*************************************************************************
 * @brief r = x0 * x1 (matrix product).
 * - r: R (matrix).
 * - x0: T0 (matrix), allows NumericalDependencyTransform.
 * - x1: T1 (matrix), allows NumericalDependencyTransform.
 *
 * Node construction should be done with the create static method.
 */
template<typename R, typename T0, typename T1> class MatrixProduct : public Value<R>
{
public:
  using Self = MatrixProduct;
  using DepT0 = typename NumericalDependencyTransform<T0>::DepType;
  using DepT1 = typename NumericalDependencyTransform<T1>::DepType;

  /// Build a new MatrixProduct node with the given output dimensions.
  static ValueRef<R> create (Context& c, NodeRefVec&& deps, const Dimension<R>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 2);
    checkNthDependencyIsValue<DepT0>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<DepT1>(typeid (Self), deps, 1);
    // Return 0 if any 0.
    if (std::any_of (deps.begin (), deps.end (), [](const NodeRef& dep) -> bool {
        return dep->hasNumericalProperty (NumericalProperty::ConstantZero);
      }))
    {
      return ConstantZero<R>::create (c, dim);
    }
    // Select node implementation
    bool identityDep0 = deps[0]->hasNumericalProperty (NumericalProperty::ConstantIdentity);
    bool identityDep1 = deps[1]->hasNumericalProperty (NumericalProperty::ConstantIdentity);
    if (identityDep0 && identityDep1)
    {
      // No specific class for Identity
      using namespace numeric;
      return NumericConstant<R>::create (c, identity (dim));
    }
    else if (identityDep0 && !identityDep1)
    {
      return Convert<R, T1>::create (c, {deps[1]}, dim);
    }
    else if (!identityDep0 && identityDep1)
    {
      return Convert<R, T0>::create (c, {deps[0]}, dim);
    }
    else
    {
      return cachedAs<Value<R> >(c, std::make_shared<Self>(std::move (deps), dim));
    }
  }

  MatrixProduct (NodeRefVec&& deps, const Dimension<R>& dim)
    : Value<R>(std::move (deps)), targetDimension_ (dim)
  {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  std::string shape() const override
  {
    return "doubleoctagon";
  }

  std::string color() const override
  {
    return "#9e9e9e";
  }


  std::string description() const override
  {
    return "Matrix Product";
  }

  // MatrixProduct additional arguments = ().
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    return dynamic_cast<const Self*>(&other) != nullptr;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<R>::create (c, targetDimension_);
    }
    const auto& x0 = this->dependency (0);
    const auto& x1 = this->dependency (1);
    auto dx0_prod = Self::create (c, {x0->derive (c, node), x1}, targetDimension_);
    auto dx1_prod = Self::create (c, {x0, x1->derive (c, node)}, targetDimension_);
    return CWiseAdd<R, std::tuple<R, R> >::create (c, {dx0_prod, dx1_prod}, targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), targetDimension_);
  }

private:
  void compute () final
  {
    auto& result = this->accessValueMutable ();
    const auto& x0 = accessValueConstCast<DepT0>(*this->dependency (0));
    const auto& x1 = accessValueConstCast<DepT1>(*this->dependency (1));
    result.noalias () =
      NumericalDependencyTransform<T0>::transform (x0) * NumericalDependencyTransform<T1>::transform (x1);
#ifdef DEBUG
    if ((x1.cols() + x1.rows() < 100 ) &&  (x0.cols() + x0.rows() < 100))
    {
      std::cerr << "=== MatrixProd === " << this << std::endl;
      std::cerr << "x0= "     << x0.rows() << " x " << x0.cols() << std::endl;
      std::cerr << "x0= " << NumericalDependencyTransform<T0>::transform (x0) << std::endl;
      std::cerr << "x1= "     << x1.rows() << " x " << x1.cols() << std::endl;
      std::cerr << "x1= " << NumericalDependencyTransform<T1>::transform (x1) << std::endl;
      std::cerr << "result= " << result << std::endl;
      std::cerr << "=== end MatrixProd === " << this << std::endl << std::endl;
    }
#endif
  }

  Dimension<R> targetDimension_;
};

/*************************************************************************
 * @brief r = n * delta + x.
 * - r: T.
 * - delta: double.
 * - x: T.
 * - n: constant int.
 * - Order of dependencies: (delta, x).
 *
 * Adds n * delta to all values (component wise) of x.
 * Used to generate x +/- delta values for numerical derivation.
 * Node construction should be done with the create static method.
 */
template<typename T> class ShiftDelta : public Value<T>
{
public:
  using Self = ShiftDelta;

  /// Build a new ShiftDelta node with the given output dimensions and shift number.
  static ValueRef<T> create (Context& c, NodeRefVec&& deps, int n, const Dimension<T>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 2);
    checkNthDependencyIsValue<double>(typeid (Self), deps, 0);
    checkNthDependencyIsValue<T>(typeid (Self), deps, 1);
    // Detect if we have a chain of ShiftDelta with the same delta.
    auto& delta = deps[0];
    auto& x = deps[1];
    auto* xAsShiftDelta = dynamic_cast<const ShiftDelta<T>*>(x.get ());
    if (xAsShiftDelta != nullptr && xAsShiftDelta->dependency (0) == delta)
    {
      // Merge with ShiftDelta dependency by summing the n.
      return Self::create (c, NodeRefVec{x->dependencies ()}, n + xAsShiftDelta->getN (), dim);
    }
    // Not a merge, select node implementation.
    if (n == 0 || delta->hasNumericalProperty (NumericalProperty::ConstantZero))
    {
      return convertRef<Value<T> >(x);
    }
    else
    {
      return cachedAs<Value<T> >(c, std::make_shared<Self>(std::move (deps), n, dim));
    }
  }

  ShiftDelta (NodeRefVec&& deps, int n, const Dimension<T>& dim)
    : Value<T>(std::move (deps)), targetDimension_ (dim), n_ (n) {}

  std::string debugInfo () const override
  {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_) +
           " n=" + std::to_string (n_);
  }

  // ShiftDelta additional arguments = (n_).
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    const auto* derived = dynamic_cast<const Self*>(&other);
    return derived != nullptr && n_ == derived->n_;
  }
  std::size_t hashAdditionalArguments () const final
  {
    std::size_t seed = 0;
    combineHash (seed, n_);
    return seed;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<T>::create (c, targetDimension_);
    }
    auto& delta = this->dependency (0);
    auto& x = this->dependency (1);
    return Self::create (c, {delta->derive (c, node), x->derive (c, node)}, n_, targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), n_, targetDimension_);
  }

  int getN () const { return n_; }

private:
  void compute () final
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& delta = accessValueConstCast<double>(*this->dependency (0));
    const auto& x = accessValueConstCast<T>(*this->dependency (1));
    result = n_ * delta + cwise (x);
  }

  Dimension<T> targetDimension_;
  int n_;
};


/*************************************************************************
 * @brief r = (1/delta)^n * sum_i coeffs_i * x_i.
 * - r: T.
 * - delta: double.
 * - x_i: T.
 * - n: constant int.
 * - coeffs_i: constant double.
 * - Order of dependencies: (delta, x_i).
 *
 * Weighted sum of dependencies, multiplied by a double.
 * Used to combine f(x+n*delta) values in numerical derivation.
 * Node construction should be done with the create static method.
 *
 * Lambda represents the 1/delta^n for a nth-order numerical derivation.
 * Note that in the whole class, coeffs[i] is the coefficient of deps[i + 1] !
 */

template<typename T> class CombineDeltaShifted : public Value<T>
{
public:
  using Self = CombineDeltaShifted;

  /// Build a new CombineDeltaShifted node with the given output dimensions, exponent and weights.
  static ValueRef<T> create (Context& c, NodeRefVec&& deps, int n, std::vector<double>&& coeffs,
                             const Dimension<T>& dim)
  {
    // Check dependencies
    checkDependenciesNotNull (typeid (Self), deps);
    checkDependencyVectorSize (typeid (Self), deps, 1 + coeffs.size ());
    checkNthDependencyIsValue<double>(typeid (Self), deps, 0);
    checkDependencyRangeIsValue<T>(typeid (Self), deps, 1, deps.size ());
    //
    auto cleanAndCreateNode = [&c, &dim](NodeRefVec&& deps2, int n2,
                                         std::vector<double>&& coeffs2) -> ValueRef<T> {
                                // Clean deps : remove constant 0, of deps with a 0 factor.
                                for (std::size_t i = 0; i < coeffs2.size ();)
                                {
                                  if (coeffs2[i] == 0. || deps2[i + 1]->hasNumericalProperty (NumericalProperty::ConstantZero))
                                  {
                                    coeffs2.erase (coeffs2.begin () + std::ptrdiff_t (i));
                                    deps2.erase (deps2.begin () + std::ptrdiff_t (i + 1));
                                  }
                                  else
                                  {
                                    ++i;
                                  }
                                }
                                // Final node selection
                                if (coeffs2.empty ())
                                {
                                  return ConstantZero<T>::create (c, dim);
                                }
                                else
                                {
                                  return cachedAs<Value<T> >(
                                    c, std::make_shared<Self>(std::move (deps2), n2, std::move (coeffs2), dim));
                                }
                              };
    // Detect if we can merge this node with its dependencies
    const auto& delta = deps[0];
    auto isSelfWithSameDelta = [&delta](const NodeRef& dep) -> bool {
                                 return dynamic_cast<const Self*>(dep.get ()) != nullptr && dep->dependency (0) == delta;
                               };
    if (!coeffs.empty () && std::all_of (deps.begin () + 1, deps.end (), isSelfWithSameDelta))
    {
      const auto depN = static_cast<const Self&>(*deps[1]).getN ();
      auto useSameNasDep1 = [depN](const NodeRef& dep) {
                              return static_cast<const Self&>(*dep).getN () == depN;
                            };
      if (std::all_of (deps.begin () + 2, deps.end (), useSameNasDep1))
      {
        /* Merge with dependencies because they use the same delta and a common N.
         *
         * V(this) = 1/delta^n * sum_i V(deps[i]) * coeffs[i].
         * V(deps[i]) = 1/delta^depN * sum_j V(depDeps[i][j]) * depCoeffs[i][j].
         * V(this) = 1/delta^(n + depN) * sum_ij V(depDeps[i][j]) * coeffs[i] * depCoeffs[i][j].
         * And simplify the sum by summing coeffs for each unique depDeps[i][j].
         */
        NodeRefVec mergedDeps{delta};
        std::vector<double> mergedCoeffs;
        for (std::size_t i = 0; i < coeffs.size (); ++i)
        {
          const auto& dep = static_cast<const Self&>(*deps[i + 1]);
          const auto& depCoeffs = dep.getCoeffs ();
          for (std::size_t j = 0; j < depCoeffs.size (); ++j)
          {
            const auto& subDep = dep.dependency (j + 1);
            auto it = std::find (mergedDeps.begin () + 1, mergedDeps.end (), subDep);
            if (it != mergedDeps.end ())
            {
              // Found
              const auto subDepIndexInMerged =
                static_cast<std::size_t>(std::distance (mergedDeps.begin () + 1, it));
              mergedCoeffs[subDepIndexInMerged] += coeffs[i] * depCoeffs[j];
            }
            else
            {
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

  CombineDeltaShifted (NodeRefVec&& deps, int n, std::vector<double>&& coeffs, const Dimension<T>& dim)
    : Value<T>(std::move (deps)), targetDimension_ (dim), coeffs_ (std::move (coeffs)), n_ (n)
  {
    assert (this->coeffs_.size () + 1 == this->nbDependencies ());
  }

  std::string description() const override
  {
    std::string s = " CombineDeltaShifted  n=" + std::to_string (n_) + " coeffs={";
    if (!coeffs_.empty ())
    {
      s += std::to_string (coeffs_[0]);
      for (std::size_t i = 1; i < coeffs_.size (); ++i)
      {
        s += ';' + std::to_string (coeffs_[i]);
      }
    }
    s += '}';
    return s;
  }

  std::string debugInfo () const override
  {
    using namespace numeric;
    std::string s = debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    return s;
  }

  // CombineDeltaShifted additional arguments = (n_, coeffs_).
  bool compareAdditionalArguments (const Node_DF& other) const final
  {
    const auto* derived = dynamic_cast<const Self*>(&other);
    return derived != nullptr && n_ == derived->n_ && coeffs_ == derived->coeffs_;
  }
  std::size_t hashAdditionalArguments () const final
  {
    std::size_t seed = 0;
    combineHash (seed, n_);
    for (const auto d : coeffs_)
    {
      combineHash (seed, d);
    }
    return seed;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    if (&node == this)
    {
      return ConstantOne<T>::create (c, targetDimension_);
    }
    // For simplicity, we assume delta is a constant with respect to derivation node.
    auto& delta = this->dependency (0);
    if (isTransitivelyDependentOn (node, *delta))
    {
      // Fail if delta is not constant for node.
      failureDeltaNotDerivable (typeid (Self));
    }
    // Derivation is a simple weighted sums of derivatives.
    const auto nbDeps = this->nbDependencies ();
    NodeRefVec derivedDeps (nbDeps);
    derivedDeps[0] = delta;
    for (std::size_t i = 1; i < nbDeps; ++i)
    {
      derivedDeps[i] = this->dependency (i)->derive (c, node);
    }
    return Self::create (c, std::move (derivedDeps), n_, std::vector<double>{coeffs_}, targetDimension_);
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    return Self::create (c, std::move (deps), n_, std::vector<double>{coeffs_}, targetDimension_);
  }

  const std::vector<double>& getCoeffs () const { return coeffs_; }

  int getN () const { return n_; }

private:
  void compute() override { compute<T>();}

  template<class U>
  typename std::enable_if<!std::is_same<U, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& delta = accessValueConstCast<double>(*this->dependency (0));
    const double lambda = pow (delta, -n_);
    result = zero (targetDimension_);
    for (std::size_t i = 0; i < coeffs_.size (); ++i)
    {
      const auto& x = accessValueConstCast<T>(*this->dependency (1 + i));
      cwise (result) += (lambda * coeffs_[i]) * cwise (x);
    }
  }

  template<class U>
  typename std::enable_if<std::is_same<U, TransitionFunction>::value, void>::type
  compute ()
  {
    using namespace numeric;
    auto& result = this->accessValueMutable ();
    const auto& delta = accessValueConstCast<double>(*this->dependency (0));
    const double lambda = pow (delta, -n_);

    std::vector<const T*> vT;
    const auto& dep = this->dependencies();

    for (size_t i = 1; i < dep.size(); i++)
    {
      vT.push_back(&accessValueConstCast<T>(*dep[i]));
    }

    result = [vT, lambda, this](const VectorLik& x)->VectorLik {
               using namespace numeric;
               VectorLik r = zero (Dimension<VectorLik>(this->targetDimension_.cols, (Eigen::Index)1));

               for (std::size_t i = 0; i < this->coeffs_.size (); ++i)
               {
                 cwise(r) += (lambda * this->coeffs_[i]) * cwise((*vT[i])(x));
               }
               return r;
             };
  }

  Dimension<T> targetDimension_;
  std::vector<double> coeffs_;
  int n_;
};

// Precompiled instantiations of numeric nodes

extern template class CWiseApply<MatrixLik, MatrixLik, TransitionFunction>;

extern template class CWiseAdd<double, std::tuple<double, double> >;
extern template class CWiseAdd<ExtendedFloat, std::tuple<ExtendedFloat, ExtendedFloat> >;
extern template class CWiseAdd<VectorLik, std::tuple<VectorLik, VectorLik> >;
extern template class CWiseAdd<RowLik, std::tuple<RowLik, RowLik> >;
extern template class CWiseAdd<MatrixLik, std::tuple<MatrixLik, MatrixLik> >;
extern template class CWiseAdd<TransitionFunction, std::tuple<TransitionFunction, TransitionFunction> >;

extern template class CWiseAdd<RowLik, MatrixLik>;
extern template class CWiseAdd<VectorLik, MatrixLik>;
extern template class CWiseAdd<DataLik, VectorLik>;
extern template class CWiseAdd<DataLik, RowLik>;

extern template class CWiseAdd<double, ReductionOf<double> >;
extern template class CWiseAdd<ExtendedFloat, ReductionOf<ExtendedFloat> >;
extern template class CWiseAdd<VectorLik, ReductionOf<VectorLik> >;
extern template class CWiseAdd<RowLik, ReductionOf<RowLik> >;
extern template class CWiseAdd<MatrixLik, ReductionOf<MatrixLik> >;
extern template class CWiseAdd<TransitionFunction, ReductionOf<TransitionFunction> >;

extern template class CWiseMean<VectorLik, ReductionOf<VectorLik>, ReductionOf<double> >;
extern template class CWiseMean<RowLik, ReductionOf<RowLik>, ReductionOf<double> >;
extern template class CWiseMean<MatrixLik, ReductionOf<MatrixLik>, ReductionOf<double> >;
extern template class CWiseMean<double, ReductionOf<double>, ReductionOf<double> >;
extern template class CWiseMean<ExtendedFloat, ReductionOf<ExtendedFloat>, ReductionOf<double> >;

extern template class CWiseMean<VectorLik, ReductionOf<VectorLik>, Eigen::VectorXd>;
extern template class CWiseMean<RowLik, ReductionOf<RowLik>,  Eigen::VectorXd>;
extern template class CWiseMean<MatrixLik, ReductionOf<MatrixLik>, Eigen::VectorXd>;
extern template class CWiseMean<VectorLik, ReductionOf<VectorLik>, Eigen::RowVectorXd>;
extern template class CWiseMean<RowLik, ReductionOf<RowLik>, Eigen::RowVectorXd>;
extern template class CWiseMean<MatrixLik, ReductionOf<MatrixLik>, Eigen::RowVectorXd>;

extern template class CWiseSub<double, std::tuple<double, double> >;
extern template class CWiseSub<ExtendedFloat, std::tuple<ExtendedFloat, ExtendedFloat> >;
extern template class CWiseSub<VectorLik, std::tuple<VectorLik, VectorLik> >;
extern template class CWiseSub<RowLik, std::tuple<RowLik, RowLik> >;
extern template class CWiseSub<MatrixLik, std::tuple<MatrixLik, MatrixLik> >;

extern template class CWiseSub<VectorLik, std::tuple<VectorLik, DataLik> >;
extern template class CWiseSub<RowLik, std::tuple<RowLik, DataLik> >;

extern template class CWiseMul<double, std::tuple<double, double> >;
extern template class CWiseMul<double, std::tuple<double, uint> >;
extern template class CWiseMul<ExtendedFloat, std::tuple<ExtendedFloat, ExtendedFloat> >;
extern template class CWiseMul<ExtendedFloat, std::tuple<ExtendedFloat, uint> >;
extern template class CWiseMul<VectorLik, std::tuple<VectorLik, VectorLik> >;
extern template class CWiseMul<RowLik, std::tuple<RowLik, RowLik> >;
extern template class CWiseMul<MatrixLik, std::tuple<MatrixLik, MatrixLik> >;

extern template class CWiseMul<RowLik, std::tuple<RowLik, Eigen::RowVectorXi> >;
extern template class CWiseMul<VectorLik, std::tuple<VectorLik, Eigen::RowVectorXi> >;
extern template class CWiseMul<VectorLik, std::tuple<DataLik, VectorLik> >;
extern template class CWiseMul<RowLik, std::tuple<DataLik, RowLik> >;
extern template class CWiseMul<MatrixLik, std::tuple<DataLik, MatrixLik> >;
extern template class CWiseMul<TransitionFunction, std::tuple<TransitionFunction, TransitionFunction> >;
extern template class CWiseMul<TransitionFunction, std::tuple<double, TransitionFunction> >;

extern template class CWiseMul<double, ReductionOf<double> >;
extern template class CWiseMul<ExtendedFloat, ReductionOf<ExtendedFloat> >;
extern template class CWiseMul<VectorLik, ReductionOf<VectorLik> >;
extern template class CWiseMul<RowLik, ReductionOf<RowLik> >;
extern template class CWiseMul<MatrixLik, ReductionOf<MatrixLik> >;

extern template class CWiseNegate<double>;
extern template class CWiseNegate<ExtendedFloat>;
extern template class CWiseNegate<VectorLik>;
extern template class CWiseNegate<RowLik>;
extern template class CWiseNegate<MatrixLik>;

extern template class CWiseInverse<double>;
extern template class CWiseInverse<ExtendedFloat>;
extern template class CWiseInverse<VectorLik>;
extern template class CWiseInverse<RowLik>;
extern template class CWiseInverse<MatrixLik>;

extern template class CWiseLog<double>;
extern template class CWiseLog<ExtendedFloat>;
extern template class CWiseLog<VectorLik>;
extern template class CWiseLog<RowLik>;
extern template class CWiseLog<MatrixLik>;

extern template class CWiseExp<double>;
extern template class CWiseExp<ExtendedFloat>;
extern template class CWiseExp<VectorLik>;
extern template class CWiseExp<RowLik>;
extern template class CWiseExp<MatrixLik>;

extern template class CWiseConstantPow<double>;
extern template class CWiseConstantPow<ExtendedFloat>;
extern template class CWiseConstantPow<VectorLik>;
extern template class CWiseConstantPow<RowLik>;
extern template class CWiseConstantPow<MatrixLik>;

extern template class ScalarProduct<DataLik, VectorLik, VectorLik>;
extern template class ScalarProduct<DataLik, RowLik, RowLik>;

extern template class LogSumExp<DataLik, VectorLik, Eigen::VectorXd>;
extern template class LogSumExp<DataLik, RowLik, Eigen::RowVectorXd>;

extern template class SumOfLogarithms<VectorLik>;
extern template class SumOfLogarithms<RowLik>;

extern template class MatrixProduct<ExtendedFloatRowVectorXd, Eigen::RowVectorXd, ExtendedFloatMatrixXd>;
extern template class MatrixProduct<Eigen::RowVectorXd, Eigen::RowVectorXd, Eigen::MatrixXd>;
extern template class MatrixProduct<MatrixLik, MatrixLik, Eigen::MatrixXd>;
extern template class MatrixProduct<MatrixLik, MatrixLik, Transposed<Eigen::MatrixXd> >;

extern template class ShiftDelta<double>;
extern template class ShiftDelta<VectorLik>;
extern template class ShiftDelta<RowLik>;
extern template class ShiftDelta<MatrixLik>;

extern template class CombineDeltaShifted<double>;
extern template class CombineDeltaShifted<VectorLik>;
extern template class CombineDeltaShifted<RowLik>;
extern template class CombineDeltaShifted<MatrixLik>;
extern template class CombineDeltaShifted<TransitionFunction>;


/*****************************************************************************
 * Numerical derivation helpers.
 */
enum class NumericalDerivativeType { Disabled, ThreePoints, FivePoints };

/// Configuration for a numerical derivation: what delta to use, and type of derivation.
struct NumericalDerivativeConfiguration
{
  NumericalDerivativeType type{NumericalDerivativeType::Disabled};
  ValueRef<double> delta{};
};

/** @brief Helper used to generate data flow expressions computing a numerical derivative.
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

template<typename NodeT, typename DepT, typename B>
ValueRef<NodeT> generateNumericalDerivative (Context& c, const NumericalDerivativeConfiguration& config,
                                             NodeRef dep,
                                             const Dimension<DepT>& depDim,
                                             Dimension<NodeT>& nodeDim,
                                             B buildNodeWithDep)
{
  if (config.delta == nullptr)
  {
    failureNumericalDerivationNotConfigured ();
  }
  switch (config.type)
  {
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

template<typename NodeT, typename B>
ValueRef<NodeT> generateNumericalDerivative (Context& c,
                                             const NumericalDerivativeConfiguration& config,
                                             NodeRef dep,
                                             const Dimension<NodeT>& nodeDim,
                                             B buildNodeWithDep)
{
  if (config.delta == nullptr)
  {
    failureNumericalDerivationNotConfigured ();
  }
  ConfiguredParameter* param = dynamic_cast<ConfiguredParameter*>(dep.get());
  if (param == nullptr)
    throw Exception("generateNumericalDerivative : dependency should be ConfiguredParameter");

  switch (config.type)
  {
  case NumericalDerivativeType::ThreePoints: {
    // Shift {-1, +1}, coeffs {-0.5, +0.5}
    auto shift_m1 = ShiftParameter::create (c, {dep->dependency(0), config.delta}, *param, -1);
    auto shift_p1 = ShiftParameter::create (c, {dep->dependency(0), config.delta}, *param, 1);
    NodeRefVec combineDeps (3);
    combineDeps[0] = config.delta;
    combineDeps[1] = buildNodeWithDep (std::move (shift_m1));
    combineDeps[2] = buildNodeWithDep (std::move (shift_p1));
    return CombineDeltaShifted<NodeT>::create (c, std::move (combineDeps), 1, {-0.5, 0.5}, nodeDim);
  } break;
  case NumericalDerivativeType::FivePoints: {
    // Shift {-2, -1, +1, +2}, coeffs {1/12, -2/3, 2/3, -1/12}
    auto shift_m2 = ShiftParameter::create (c, {dep->dependency(0), config.delta}, *param, -2);
    auto shift_m1 = ShiftParameter::create (c, {dep->dependency(0), config.delta}, *param, -1);
    auto shift_p1 = ShiftParameter::create (c, {dep->dependency(0), config.delta}, *param, 1);
    auto shift_p2 = ShiftParameter::create (c, {dep->dependency(0), config.delta}, *param, 2);
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
} // namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_DATAFLOWCWISECOMPUTING_H
