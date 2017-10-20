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
	 * - deduce size (require 1 arg in reductions) ?
	 * - opt: merge * of *, + of +, ... ?
	 * - common template for *,+ opts ?
	 */

	/******************************** Utils *******************************/
	namespace {
		// Constant builders
		auto zeroVector (SizeType size) -> decltype (VectorDouble::Zero (size)) {
			return VectorDouble::Zero (size);
		}
		auto oneVector (SizeType size) -> decltype (VectorDouble::Ones (size)) {
			return VectorDouble::Ones (size);
		}
		auto zeroMatrix (MatrixDimension dim) -> decltype (MatrixDouble::Zero (dim.rows, dim.cols)) {
			return MatrixDouble::Zero (dim.rows, dim.cols);
		}
		auto oneMatrix (MatrixDimension dim) -> decltype (MatrixDouble::Ones (dim.rows, dim.cols)) {
			return MatrixDouble::Ones (dim.rows, dim.cols);
		}

		// Constant checking predicate (const T & -> bool)
		bool isExactZeroVector (const VectorDouble & v) { return v == zeroVector (dimensions (v)); }
		bool isExactOnesVector (const VectorDouble & v) { return v == oneVector (dimensions (v)); }
		bool isExactZeroMatrix (const MatrixDouble & m) { return m == zeroMatrix (dimensions (m)); }
		bool isExactOnesMatrix (const MatrixDouble & m) { return m == oneMatrix (dimensions (m)); }
		bool isExactIdentityMatrix (const MatrixDouble & m) {
			return m.rows () == m.cols () && m == MatrixDouble::Identity (m.rows (), m.cols ());
		}

		// Constant node checking predicates (const NodeRef & -> bool)
		auto isConstantZeroDouble = predicateIsConstantValueMatching<double> (0.);
		auto isConstantOneDouble = predicateIsConstantValueMatching<double> (1.);
		auto isConstantZeroVector = predicateIsConstantValueMatching<VectorDouble> (isExactZeroVector);
		auto isConstantOnesVector = predicateIsConstantValueMatching<VectorDouble> (isExactOnesVector);
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

	// Dimensions
	std::string MatrixDimension::toString () const {
		return "(" + std::to_string (rows) + "," + std::to_string (cols) + ")";
	}

	/******************************** External specialisations *******************************/

	// Debug info specialisations
	template <> std::string Value<double>::debugInfo () const {
		return std::to_string (this->accessValue ());
	}
	template <> std::string Value<VectorDouble>::debugInfo () const {
		auto & v = this->accessValue ();
		auto s = std::string ("size=") + std::to_string (dimensions (v)) + " props{";
		if (isExactZeroVector (v))
			s += '0';
		if (isExactOnesVector (v))
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

	// Parameter<double> specialisation
	template <> NodeRef Parameter<double>::derive (const Node & node) {
		if (&node == static_cast<const Node *> (this)) {
			return Builder<Constant<double>>::makeOne ();
		} else {
			return Builder<Constant<double>>::makeZero ();
		}
	}

	// Constant<double> specialisation
	template <> NodeRef Constant<double>::derive (const Node &) {
		return Builder<Constant<double>>::makeZero ();
	}
	std::shared_ptr<Constant<double>> Builder<Constant<double>>::make (double d) {
		if (d == 0.) {
			return makeZero ();
		} else if (d == 1.) {
			return makeOne ();
		} else {
			return std::make_shared<Constant<double>> (d);
		}
	}
	std::shared_ptr<Constant<double>> Builder<Constant<double>>::makeZero () {
		static auto zero = std::make_shared<Constant<double>> (0.);
		return zero;
	}
	std::shared_ptr<Constant<double>> Builder<Constant<double>>::makeOne () {
		static auto one = std::make_shared<Constant<double>> (1.);
		return one;
	}

	// Constant<VectorDouble> specialisation
	template <> NodeRef Constant<VectorDouble>::derive (const Node &) {
		return Builder<Constant<VectorDouble>>::makeZero (dimensions (*this));
	}
	std::shared_ptr<Constant<VectorDouble>>
	Builder<Constant<VectorDouble>>::makeZero (SizeType size) {
		return make (zeroVector (size));
	}

	// Constant<MatrixDouble> specialisation
	template <> NodeRef Constant<MatrixDouble>::derive (const Node &) {
		return Builder<Constant<MatrixDouble>>::makeZero (dimensions (*this));
	}
	std::shared_ptr<Constant<MatrixDouble>>
	Builder<Constant<MatrixDouble>>::makeZero (const MatrixDimension & dim) {
		return make (zeroMatrix (dim));
	}
	std::shared_ptr<Constant<MatrixDouble>>
	Builder<Constant<MatrixDouble>>::makeOne (const MatrixDimension & dim) {
		return make (oneMatrix (dim));
	}

	/******************************** New Nodes *******************************/
	// TODO templatize impls

	/******************************** Old nodes *******************************/

	// AddDouble
	AddDouble::AddDouble (NodeRefVec && deps) : Value<double> (std::move (deps)) {
		checkDependencies (*this);
	}
	void AddDouble::compute () {
		callWithValues (*this, [](double & r) { r = 0.; }, [](double & r, double d) { r += d; });
	}
	NodeRef AddDouble::derive (const Node & node) {
		return makeNode<AddDouble> (
		    this->dependencies ().map ([&node](const NodeRef & dep) { return dep->derive (node); }));
	}
	ValueRef<double> Builder<AddDouble>::make (NodeRefVec && deps) {
		checkDependencies<AddDouble> (deps);
		// Remove '0s' from deps
		removeDependenciesIf (deps, isConstantZeroDouble);
		// Node choice
		if (deps.size () == 1) {
			return convertRef<Value<double>> (deps[0]);
		} else if (deps.size () == 0) {
			return makeNode<Constant<double>> (0.);
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
			addDeps.emplace_back (makeNode<MulDouble> (std::move (mulDeps)));
		}
		return makeNode<AddDouble> (std::move (addDeps));
	}
	ValueRef<double> Builder<MulDouble>::make (NodeRefVec && deps) {
		checkDependencies<MulDouble> (deps);
		// Return 0 if any dep is 0
		if (std::any_of (deps.begin (), deps.end (), isConstantZeroDouble)) {
			return makeNode<Constant<double>> (0.);
		}
		// Remove any 1s
		removeDependenciesIf (deps, isConstantOneDouble);
		// Node choice
		if (deps.size () == 1) {
			return convertRef<Value<double>> (deps[0]);
		} else if (deps.size () == 0) {
			return makeNode<Constant<double>> (1.);
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
		auto dLhs = makeNode<ScalarProdDouble> ({lhs->derive (node), rhs});
		auto dRhs = makeNode<ScalarProdDouble> ({lhs, rhs->derive (node)});
		return makeNode<AddDouble> ({std::move (dLhs), std::move (dRhs)});
	}
	ValueRef<double> Builder<ScalarProdDouble>::make (NodeRefVec && deps) {
		checkDependencies<ScalarProdDouble> (deps);
		auto & lhs = deps[0];
		auto & rhs = deps[1];
		if (isConstantZeroVector (lhs) || isConstantZeroVector (rhs)) {
			return makeNode<Constant<double>> (0.);
		} else {
			return std::make_shared<ScalarProdDouble> (std::move (deps));
		}
	}

	// AddVectorDouble
	AddVectorDouble::AddVectorDouble (NodeRefVec && deps, SizeType size)
	    : Value<VectorDouble> (std::move (deps), size) {
		checkDependencies (*this);
	}
	void AddVectorDouble::compute () {
		callWithValues (*this, [](VectorDouble & r) { r.fill (0.); },
		                [](VectorDouble & r, const VectorDouble & v) { r += v; });
	}
	NodeRef AddVectorDouble::derive (const Node & node) {
		return makeNode<AddVectorDouble> (
		    this->dependencies ().map ([&node](const NodeRef & dep) { return dep->derive (node); }),
		    dimensions (*this));
	}
	ValueRef<VectorDouble> Builder<AddVectorDouble>::make (NodeRefVec && deps, SizeType size) {
		checkDependencies<AddVectorDouble> (deps);
		// Remove Os
		removeDependenciesIf (deps, isConstantZeroVector);
		// Select node impl
		if (deps.size () == 1) {
			return convertRef<Value<VectorDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return Builder<Constant<VectorDouble>>::makeZero (size);
		} else {
			return std::make_shared<AddVectorDouble> (std::move (deps), size);
		}
	}

	// CWiseInverseVectorDouble
	CWiseInverseVectorDouble::CWiseInverseVectorDouble (NodeRefVec && deps, SizeType size)
	    : Value<VectorDouble> (std::move (deps), size) {
		checkDependencies (*this);
	}
	void CWiseInverseVectorDouble::compute () {
		callWithValues (*this, [](VectorDouble & r, const VectorDouble & v) { r = v.cwiseInverse (); });
	}
	ValueRef<VectorDouble> Builder<CWiseInverseVectorDouble>::make (NodeRefVec && deps,
	                                                                SizeType size) {
		checkDependencies<CWiseInverseVectorDouble> (deps);
		auto & arg = deps[0];
		if (isConstantOnesVector (arg)) {
			return convertRef<Value<VectorDouble>> (arg);
		} else {
			return std::make_shared<CWiseInverseVectorDouble> (std::move (deps), size);
		}
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
		return makeNode<AddMatrixDouble> (
		    this->dependencies ().map ([&node](const NodeRef & dep) { return dep->derive (node); }),
		    dimensions (*this));
	}
	ValueRef<MatrixDouble> Builder<AddMatrixDouble>::make (NodeRefVec && deps, MatrixDimension dim) {
		checkDependencies<AddMatrixDouble> (deps);
		// Remove Os
		removeDependenciesIf (deps, isConstantZeroMatrix);
		// Select node impl
		if (deps.size () == 1) {
			return convertRef<Value<MatrixDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return Builder<Constant<MatrixDouble>>::makeZero (dim);
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
		auto dLhs = makeNode<MulMatrixDouble> ({lhs->derive (node), rhs}, dim);
		auto dRhs = makeNode<MulMatrixDouble> ({lhs, rhs->derive (node)}, dim);
		return makeNode<AddMatrixDouble> ({std::move (dLhs), std::move (dRhs)}, dim);
	}
	ValueRef<MatrixDouble> Builder<MulMatrixDouble>::make (NodeRefVec && deps, MatrixDimension dim) {
		checkDependencies<MulMatrixDouble> (deps);
		auto & lhs = deps[0];
		auto & rhs = deps[1];
		if (isConstantZeroMatrix (lhs) || isConstantZeroMatrix (rhs)) {
			return Builder<Constant<MatrixDouble>>::makeZero (dim);
		} else if (isConstantIdentityMatrix (lhs)) {
			return convertRef<Value<MatrixDouble>> (rhs);
		} else if (isConstantIdentityMatrix (rhs)) {
			return convertRef<Value<MatrixDouble>> (lhs);
		} else {
			return std::make_shared<MulMatrixDouble> (std::move (deps), dim);
		}
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
			addDeps.emplace_back (makeNode<CWiseMulMatrixDouble> (std::move (mulDeps), dim));
		}
		return makeNode<AddMatrixDouble> (std::move (addDeps), dim);
	}
	ValueRef<MatrixDouble> Builder<CWiseMulMatrixDouble>::make (NodeRefVec && deps,
	                                                            MatrixDimension dim) {
		checkDependencies<CWiseMulMatrixDouble> (deps);
		// Return 0 if any 0 dep
		if (std::any_of (deps.begin (), deps.end (), isConstantZeroMatrix)) {
			return Builder<Constant<MatrixDouble>>::makeZero (dim);
		}
		// Remove 1s
		removeDependenciesIf (deps, isConstantOnesMatrix);
		// Select node
		if (deps.size () == 1) {
			return convertRef<Value<MatrixDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return Builder<Constant<MatrixDouble>>::makeOne (dim);
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
		auto dLhs = makeNode<MulTransposedMatrixVectorDouble> ({lhs->derive (node), rhs}, dim);
		auto dRhs = makeNode<MulTransposedMatrixVectorDouble> ({lhs, rhs->derive (node)}, dim);
		return makeNode<AddVectorDouble> ({std::move (dLhs), std::move (dRhs)}, dim);
	}
	ValueRef<VectorDouble> Builder<MulTransposedMatrixVectorDouble>::make (NodeRefVec && deps,
	                                                                       SizeType size) {
		checkDependencies<MulTransposedMatrixVectorDouble> (deps);
		auto & lhs = deps[0];
		auto & rhs = deps[1];
		if (isConstantZeroMatrix (lhs) || isConstantZeroVector (rhs)) {
			return Builder<Constant<VectorDouble>>::makeZero (size);
		} else if (isConstantIdentityMatrix (lhs)) {
			return convertRef<Value<VectorDouble>> (rhs);
		} else {
			return std::make_shared<MulTransposedMatrixVectorDouble> (std::move (deps), size);
		}
	}

	// MulScalarMatrixDouble
	MulScalarMatrixDouble::MulScalarMatrixDouble (NodeRefVec && deps, MatrixDimension dim)
	    : Value<MatrixDouble> (std::move (deps), dim.rows, dim.cols) {
		checkDependencies (*this);
	}
	void MulScalarMatrixDouble::compute () {
		callWithValues (
		    *this, [](MatrixDouble & r, double d, const MatrixDouble & m) { r.noalias () = d * m; });
	}
	NodeRef MulScalarMatrixDouble::derive (const Node & node) {
		auto dim = dimensions (*this);
		auto & lhs = this->dependencies ()[0];
		auto & rhs = this->dependencies ()[1];
		auto dLhs = makeNode<MulScalarMatrixDouble> ({lhs->derive (node), rhs}, dim);
		auto dRhs = makeNode<MulScalarMatrixDouble> ({lhs, rhs->derive (node)}, dim);
		return makeNode<AddMatrixDouble> ({std::move (dLhs), std::move (dRhs)}, dim);
	}
	ValueRef<MatrixDouble> Builder<MulScalarMatrixDouble>::make (NodeRefVec && deps,
	                                                             MatrixDimension dim) {
		checkDependencies<MulScalarMatrixDouble> (deps);
		auto & lhs = deps[0];
		auto & rhs = deps[1];
		if (isConstantZeroDouble (lhs) || isConstantZeroMatrix (rhs)) {
			return Builder<Constant<MatrixDouble>>::makeZero (dim);
		} else if (isConstantOneDouble (lhs)) {
			return convertRef<Value<MatrixDouble>> (rhs);
		} else {
			return std::make_shared<MulScalarMatrixDouble> (std::move (deps), dim);
		}
	}
} // namespace DF
} // namespace bpp
