//
// File: DataFlowNumeric.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-10-09 00:00:00
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

#include <Bpp/Exceptions.h> // checks
#include <Bpp/NewPhyl/DataFlowInternalTemplates.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/NewPhyl/DataFlowTemplates.h>
#include <Bpp/NewPhyl/Debug.h> // checks
#include <Bpp/NewPhyl/Range.h>
#include <memory>
#include <typeinfo>
#include <utility>

namespace bpp {
namespace DF {
	/* TODO
	 * - opt: aggregate constants
	 * - merge constants by type
	 * - support for general register (reuse code in DF.cpp, declare in DFTinternal).
	 * - opt: keep 0s instead of recreating them
	 * - deduce size (require 1 arg in reductions)
	 */

	// Utils
	namespace {
		// Constant builders
		auto zeroVector (SizeType size) -> decltype (VectorDouble::Zero (size)) {
			return VectorDouble::Zero (size);
		}
		auto zeroMatrix (MatrixDimension dim) -> decltype (MatrixDouble::Zero (dim.rows, dim.cols)) {
			return MatrixDouble::Zero (dim.rows, dim.cols);
		}
		auto onesMatrix (MatrixDimension dim) -> decltype (MatrixDouble::Ones (dim.rows, dim.cols)) {
			return MatrixDouble::Ones (dim.rows, dim.cols);
		}

		// Constant checking predicate (const T & -> bool)
		bool isExactZeroVector (const VectorDouble & v) { return v == zeroVector (dimensions (v)); }
		bool isExactZeroMatrix (const MatrixDouble & m) { return m == zeroMatrix (dimensions (m)); }
		bool isExactOnesMatrix (const MatrixDouble & m) { return m == onesMatrix (dimensions (m)); }
		bool isExactIdentityMatrix (const MatrixDouble & m) {
			return m.rows () == m.cols () && m == MatrixDouble::Identity (m.rows (), m.cols ());
		}

		// Constant node checking predicates (const NodeRef & -> bool)
		auto isConstantZeroDouble = predicateIsConstantValueMatching<double> (0.);
		auto isConstantOneDouble = predicateIsConstantValueMatching<double> (1.);
		auto isConstantZeroVector = predicateIsConstantValueMatching<VectorDouble> (isExactZeroVector);
		auto isConstantZeroMatrix = predicateIsConstantValueMatching<MatrixDouble> (isExactZeroMatrix);
		auto isConstantOnesMatrix = predicateIsConstantValueMatching<MatrixDouble> (isExactOnesMatrix);
		auto isConstantIdentityMatrix =
		    predicateIsConstantValueMatching<MatrixDouble> (isExactIdentityMatrix);

		// Matrix dimension checking TODO drop ?
		void checkDepsHaveRequiredDimension (const std::type_info & inNodeType, const NodeRefVec & deps,
		                                     MatrixDimension dim) {
			for (auto i : index_range (deps)) {
				auto depDim = dimensions (dynamic_cast<const Value<MatrixDouble> &> (*deps[i]));
				if (dim != depDim)
					throw Exception (prettyTypeName (inNodeType) + ": matrix size mismatch for " +
					                 std::to_string (i) + "-th argument: expected " + dim.toString () +
					                 ", got " + depDim.toString ());
			}
		}
		template <typename NodeType>
		void checkDepsHaveRequiredDimension (const NodeType & node, MatrixDimension dim) {
			checkDepsHaveRequiredDimension (typeid (NodeType), node.dependencies (), dim);
		}
	} // namespace

	// Debug info specialisations
	template <> std::string Value<double>::debugInfo () const {
		return std::to_string (this->accessValue ());
	}
	template <> std::string Value<VectorDouble>::debugInfo () const {
		auto & v = this->accessValue ();
		auto s = std::string ("size=") + std::to_string (dimensions (v)) + " props{";
		if (isExactZeroMatrix (v))
			s += '0';
		if (isExactOnesMatrix (v))
			s += '1';
		s += "}";
		return s;
	}
	template <> std::string Value<MatrixDouble>::debugInfo () const {
		auto & m = this->accessValue ();
		auto s = std::string ("dim") + dimensions (m).toString () + " props{";
		if (isExactZeroMatrix (m))
			s += '0';
		if (isExactOnesMatrix (m))
			s += '1';
		if (isExactIdentityMatrix (m))
			s += 'I';
		s += "}";
		return s;
	}

	// Dimensions
	std::string MatrixDimension::toString () const {
		return "(" + std::to_string (rows) + "," + std::to_string (cols) + ")";
	}

	// Parameter specialisations
	template <> NodeRef Parameter<double>::derive (const Node & node) {
		if (&node == static_cast<const Node *> (this)) {
			return ConstantDouble::one;
		} else {
			return ConstantDouble::zero;
		}
	}

	// ConstantDouble
	ConstantDouble::ConstantDouble (double d) : Value<double> (noDependency, d) {
		this->makeValid ();
	}
	void ConstantDouble::compute () { failureComputeWasCalled (typeid (ConstantDouble)); }
	bool ConstantDouble::isConstant () const { return true; }
	NodeRef ConstantDouble::derive (const Node &) { return zero; }
	std::shared_ptr<ConstantDouble> ConstantDouble::zero = std::make_shared<ConstantDouble> (0.);
	std::shared_ptr<ConstantDouble> ConstantDouble::one = std::make_shared<ConstantDouble> (1.);
	std::shared_ptr<ConstantDouble> ConstantDouble::create (double d) {
		if (d == 0.) {
			return zero;
		} else if (d == 1.) {
			return one;
		} else {
			return std::make_shared<ConstantDouble> (d);
		}
	}

	// AddDouble
	AddDouble::AddDouble (NodeRefVec && deps) : Value<double> (std::move (deps)) {
		checkDependencies (*this);
	}
	void AddDouble::compute () {
		callWithValues (*this, [](double & r) { r = 0.; }, [](double & r, double d) { r += d; });
	}
	NodeRef AddDouble::derive (const Node & node) {
		return AddDouble::create (
		    this->dependencies ().map ([&node](const NodeRef & dep) { return dep->derive (node); }));
	}
	ValueRef<double> AddDouble::create (NodeRefVec && deps) {
		checkDependencies<Dependencies> (deps, typeid (AddDouble));
		// Remove '0s' from deps
		removeDependenciesIf (deps, isConstantZeroDouble);
		// Node choice
		if (deps.size () == 1) {
			return convertRef<Value<double>> (deps[0]);
		} else if (deps.size () == 0) {
			return ConstantDouble::zero;
		} else {
			return std::make_shared<AddDouble> (std::move (deps));
		}
	}

	// MulDouble
	MulDouble::MulDouble (NodeRefVec && deps) : Value<double> (std::move (deps)) {
		checkDependencies (*this);
	}
	void MulDouble::compute () {
		callWithValues (*this, [](double & r) { r = 1.; }, [](double & r, double d) { r *= d; });
	}
	NodeRef MulDouble::derive (const Node & node) {
		NodeRefVec addDeps;
		for (auto i : bpp::index_range (this->dependencies ())) {
			NodeRefVec mulDeps = this->dependencies ();
			mulDeps[i] = this->dependencies ()[i]->derive (node);
			addDeps.emplace_back (MulDouble::create (std::move (mulDeps)));
		}
		return AddDouble::create (std::move (addDeps));
	}
	ValueRef<double> MulDouble::create (NodeRefVec && deps) {
		checkDependencies<Dependencies> (deps, typeid (MulDouble));
		// Return 0 if any dep is 0
		if (std::any_of (deps.begin (), deps.end (), isConstantZeroDouble)) {
			return ConstantDouble::zero;
		}
		// Remove any 1s
		removeDependenciesIf (deps, isConstantOneDouble);
		// Node choice
		if (deps.size () == 1) {
			return convertRef<Value<double>> (deps[0]);
		} else if (deps.size () == 0) {
			return ConstantDouble::one;
		} else {
			return std::make_shared<MulDouble> (std::move (deps));
		}
	}

	// ScalarProdDouble
	ScalarProdDouble::ScalarProdDouble (NodeRefVec && deps) : Value<double> (std::move (deps)) {
		checkDependencies (*this);
	}
	void ScalarProdDouble::compute () {
		callWithValues (*this, [](double & r, const VectorDouble & lhs, const VectorDouble & rhs) {
			r = lhs.dot (rhs);
		});
	}
	NodeRef ScalarProdDouble::derive (const Node & node) {
		auto & lhs = this->dependencies ()[0];
		auto & rhs = this->dependencies ()[1];
		auto dLhs = ScalarProdDouble::create (NodeRefVec{lhs->derive (node), rhs});
		auto dRhs = ScalarProdDouble::create (NodeRefVec{lhs, rhs->derive (node)});
		return AddDouble::create (NodeRefVec{std::move (dLhs), std::move (dRhs)});
	}
	ValueRef<double> ScalarProdDouble::create (NodeRefVec && deps) {
		checkDependencies<Dependencies> (deps, typeid (ScalarProdDouble));
		auto & lhs = deps[0];
		auto & rhs = deps[1];
		// Return 0 if any dep is 0
		if (isConstantZeroVector (lhs) || isConstantZeroVector (rhs)) {
			return ConstantDouble::zero;
		}
		return std::make_shared<ScalarProdDouble> (std::move (deps));
	}

	// ConstantVectorDouble
	void ConstantVectorDouble::compute () { failureComputeWasCalled (typeid (ConstantVectorDouble)); }
	NodeRef ConstantVectorDouble::derive (const Node &) { return createZero (dimensions (*this)); }
	std::shared_ptr<ConstantVectorDouble> ConstantVectorDouble::createZero (SizeType size) {
		return create (zeroVector (size));
	}

	// AddVectorDouble
	AddVectorDouble::AddVectorDouble (NodeRefVec && deps, SizeType size)
	    : Value<VectorDouble> (std::move (deps), size) {
		checkDependencies (*this);
	}
	void AddVectorDouble::compute () {
		callWithValues (*this, [](VectorDouble & r) { r.fill (0.); },
		                [](VectorDouble & r, const VectorDouble & m) { r += m; });
	}
	NodeRef AddVectorDouble::derive (const Node & node) {
		return AddVectorDouble::create (
		    this->dependencies ().map ([&node](const NodeRef & dep) { return dep->derive (node); }),
		    dimensions (*this));
	}
	ValueRef<VectorDouble> AddVectorDouble::create (NodeRefVec && deps, SizeType size) {
		checkDependencies<Dependencies> (deps, typeid (AddVectorDouble));
		// Remove Os
		removeDependenciesIf (deps, isConstantZeroVector);
		// Select node impl
		if (deps.size () == 1) {
			return convertRef<Value<VectorDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return ConstantVectorDouble::createZero (size);
		} else {
			return std::make_shared<AddVectorDouble> (std::move (deps), size);
		}
	}

	// ConstantMatrixDouble
	void ConstantMatrixDouble::compute () { failureComputeWasCalled (typeid (ConstantMatrixDouble)); }
	NodeRef ConstantMatrixDouble::derive (const Node &) { return createZero (dimensions (*this)); }
	std::shared_ptr<ConstantMatrixDouble> ConstantMatrixDouble::createZero (MatrixDimension dim) {
		return create (zeroMatrix (dim));
	}

	// AddMatrixDouble
	AddMatrixDouble::AddMatrixDouble (NodeRefVec && deps, MatrixDimension dim)
	    : Value<MatrixDouble> (std::move (deps), dim.rows, dim.cols) {
		checkDependencies (*this);
		checkDepsHaveRequiredDimension (*this, dim);
	}
	void AddMatrixDouble::compute () {
		callWithValues (*this, [](MatrixDouble & r) { r.fill (0.); },
		                [](MatrixDouble & r, const MatrixDouble & m) { r += m; });
	}
	NodeRef AddMatrixDouble::derive (const Node & node) {
		return AddMatrixDouble::create (
		    this->dependencies ().map ([&node](const NodeRef & dep) { return dep->derive (node); }),
		    dimensions (*this));
	}
	ValueRef<MatrixDouble> AddMatrixDouble::create (NodeRefVec && deps, MatrixDimension dim) {
		checkDependencies<Dependencies> (deps, typeid (AddMatrixDouble));
		// Remove Os
		removeDependenciesIf (deps, isConstantZeroMatrix);
		// Select node impl
		if (deps.size () == 1) {
			return convertRef<Value<MatrixDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return ConstantMatrixDouble::createZero (dim);
		} else {
			return std::make_shared<AddMatrixDouble> (std::move (deps), dim);
		}
	}

	// MulMatrixDouble
	MulMatrixDouble::MulMatrixDouble (NodeRefVec && deps, MatrixDimension dim)
	    : Value<MatrixDouble> (std::move (deps), dim.rows, dim.cols) {
		checkDependencies (*this);
	}
	void MulMatrixDouble::compute () {
		callWithValues (*this, [](MatrixDouble & r, const MatrixDouble & lhs,
		                          const MatrixDouble & rhs) { r.noalias () = lhs * rhs; });
	}
	NodeRef MulMatrixDouble::derive (const Node & node) {
		auto dim = dimensions (*this);
		auto & lhs = this->dependencies ()[0];
		auto & rhs = this->dependencies ()[1];
		auto dLhs = MulMatrixDouble::create (NodeRefVec{lhs->derive (node), rhs}, dim);
		auto dRhs = MulMatrixDouble::create (NodeRefVec{lhs, rhs->derive (node)}, dim);
		return AddMatrixDouble::create (NodeRefVec{std::move (dLhs), std::move (dRhs)}, dim);
	}
	ValueRef<MatrixDouble> MulMatrixDouble::create (NodeRefVec && deps, MatrixDimension dim) {
		checkDependencies<Dependencies> (deps, typeid (MulMatrixDouble));
		auto & lhs = deps[0];
		auto & rhs = deps[1];
		// Return O if any matrix is 0
		if (isConstantZeroMatrix (lhs) || isConstantZeroMatrix (rhs)) {
			return ConstantMatrixDouble::createZero (dim);
		}
		// Return the other if one is identity
		if (isConstantIdentityMatrix (lhs)) {
			return convertRef<Value<MatrixDouble>> (rhs);
		}
		if (isConstantIdentityMatrix (rhs)) {
			return convertRef<Value<MatrixDouble>> (lhs);
		}
		return std::make_shared<MulMatrixDouble> (std::move (deps), dim);
	}
	ValueRef<MatrixDouble> MulMatrixDouble::create (ValueRef<MatrixDouble> lhs,
	                                                ValueRef<MatrixDouble> rhs, MatrixDimension dim) {
		return create (NodeRefVec{std::move (lhs), std::move (rhs)}, dim);
	}

	// CWiseMulMatrixDouble
	CWiseMulMatrixDouble::CWiseMulMatrixDouble (NodeRefVec && deps, MatrixDimension dim)
	    : Value<MatrixDouble> (std::move (deps), dim.rows, dim.cols) {
		checkDependencies (*this);
		checkDepsHaveRequiredDimension (*this, dim);
	}
	void CWiseMulMatrixDouble::compute () {
		callWithValues (*this, [](MatrixDouble & r) { r.fill (1.); },
		                [](MatrixDouble & r, const MatrixDouble & m) { r = r.cwiseProduct (m); });
	}
	NodeRef CWiseMulMatrixDouble::derive (const Node & node) {
		auto dim = dimensions (*this);
		NodeRefVec addDeps;
		for (auto i : bpp::index_range (this->dependencies ())) {
			NodeRefVec mulDeps = this->dependencies ();
			mulDeps[i] = this->dependencies ()[i]->derive (node);
			addDeps.emplace_back (CWiseMulMatrixDouble::create (std::move (mulDeps), dim));
		}
		return AddMatrixDouble::create (std::move (addDeps), dim);
	}
	ValueRef<MatrixDouble> CWiseMulMatrixDouble::create (NodeRefVec && deps, MatrixDimension dim) {
		checkDependencies<Dependencies> (deps, typeid (CWiseMulMatrixDouble));
		// Return 0 if any 0 dep
		if (std::any_of (deps.begin (), deps.end (), isConstantZeroMatrix)) {
			return ConstantMatrixDouble::createZero (dim);
		}
		// Remove 1s
		removeDependenciesIf (deps, isConstantOnesMatrix);
		// Select node
		if (deps.size () == 1) {
			return convertRef<Value<MatrixDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return ConstantMatrixDouble::create (onesMatrix (dim));
		} else {
			return std::make_shared<CWiseMulMatrixDouble> (std::move (deps), dim);
		}
	}

	// MulTransposedMatrixVectorDouble
	MulTransposedMatrixVectorDouble::MulTransposedMatrixVectorDouble (NodeRefVec && deps,
	                                                                  SizeType size)
	    : Value<VectorDouble> (std::move (deps), size) {
		checkDependencies (*this);
	}
	void MulTransposedMatrixVectorDouble::compute () {
		callWithValues (*this, [](VectorDouble & r, const MatrixDouble & lhs,
		                          const VectorDouble & rhs) { r.noalias () = lhs.transpose () * rhs; });
	}
	NodeRef MulTransposedMatrixVectorDouble::derive (const Node & node) {
		auto dim = dimensions (*this);
		auto & lhs = this->dependencies ()[0];
		auto & rhs = this->dependencies ()[1];
		auto dLhs = MulTransposedMatrixVectorDouble::create (NodeRefVec{lhs->derive (node), rhs}, dim);
		auto dRhs = MulTransposedMatrixVectorDouble::create (NodeRefVec{lhs, rhs->derive (node)}, dim);
		return AddVectorDouble::create (NodeRefVec{std::move (dLhs), std::move (dRhs)}, dim);
	}
	ValueRef<VectorDouble> MulTransposedMatrixVectorDouble::create (NodeRefVec && deps,
	                                                                SizeType size) {
		checkDependencies<Dependencies> (deps, typeid (MulTransposedMatrixVectorDouble));
		auto & lhs = deps[0];
		auto & rhs = deps[1];
		// Return O if any matrix is 0
		if (isConstantZeroMatrix (lhs) || isConstantZeroVector (rhs)) {
			return ConstantVectorDouble::createZero (size);
		}
		// Return the other if one is identity
		return std::make_shared<MulTransposedMatrixVectorDouble> (std::move (deps), size);
	}
	ValueRef<VectorDouble> MulTransposedMatrixVectorDouble::create (ValueRef<MatrixDouble> lhs,
	                                                                ValueRef<VectorDouble> rhs,
	                                                                SizeType size) {
		return create (NodeRefVec{std::move (lhs), std::move (rhs)}, size);
	}
} // namespace DF
} // namespace bpp
