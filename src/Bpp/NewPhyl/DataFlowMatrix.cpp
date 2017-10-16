//
// File: DataFlowMatrix.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-10-11
// Last modified: 2017-10-11
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
#include <Bpp/NewPhyl/DataFlowMatrix.h>
#include <Bpp/NewPhyl/DataFlowTemplateUtils.h>
#include <Bpp/NewPhyl/Debug.h> // checks
#include <Bpp/NewPhyl/Range.h> // checks
#include <typeinfo>

namespace bpp {
namespace DF {
	// Dimensions
	std::string MatrixDimension::toString () const {
		return "(" + std::to_string (rows) + "," + std::to_string (cols) + ")";
	}

	namespace {
		auto zeroMatrix (MatrixDimension dim) -> decltype (MatrixDouble::Zero (dim.rows, dim.cols)) {
			return MatrixDouble::Zero (dim.rows, dim.cols);
		}
		auto onesMatrix (MatrixDimension dim) -> decltype (MatrixDouble::Ones (dim.rows, dim.cols)) {
			return MatrixDouble::Ones (dim.rows, dim.cols);
		}

		bool isExactZeroMatrix (const MatrixDouble & m) { return m == zeroMatrix (dimensions (m)); }
		bool isExactOnesMatrix (const MatrixDouble & m) { return m == onesMatrix (dimensions (m)); }
		bool isExactIdentityMatrix (const MatrixDouble & m) {
			return m.rows () == m.cols () && m == MatrixDouble::Identity (m.rows (), m.cols ());
		}

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

	std::string debugInfoFor (const MatrixDouble & m) {
		auto s = std::string ("dim=") + dimensions (m).toString () + " props=";
		if (isExactZeroMatrix (m))
			s += '0';
		if (isExactOnesMatrix (m))
			s += '1';
		if (isExactIdentityMatrix (m))
			s += 'I';
		return s;
	}

	// ConstantMatrixDouble
	void ConstantMatrixDouble::compute () { failureComputeWasCalled (typeid (ConstantMatrixDouble)); }
	std::string ConstantMatrixDouble::debugInfo () const {
		return debugInfoFor (this->accessValue ());
	}
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
	std::string AddMatrixDouble::debugInfo () const { return debugInfoFor (this->accessValue ()); }
	NodeRef AddMatrixDouble::derive (const Node & node) {
		return AddMatrixDouble::create (
		    this->dependencies ().map ([&node](const NodeRef & dep) { return dep->derive (node); }),
		    dimensions (*this));
	}
	ValueRef<MatrixDouble> AddMatrixDouble::create (NodeRefVec && deps, MatrixDimension dim) {
		checkDependencies<Dependencies> (deps, typeid (AddMatrixDouble));
		// Remove Os
		removeDependenciesIf (deps, predicateIsConstantValueMatching<MatrixDouble> (isExactZeroMatrix));
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
		// TODO extract as func & improve
		auto & lhs = accessValueUncheckedCast<MatrixDouble> (*this->dependencies ()[0]);
		auto & rhs = accessValueUncheckedCast<MatrixDouble> (*this->dependencies ()[1]);
		if (!(lhs.cols () == rhs.rows () && lhs.rows () == dim.rows && rhs.cols () == dim.cols))
			throw Exception (prettyTypeName<MulMatrixDouble> () +
			                 ": matrix size mismatch: " + dimensions (lhs).toString () + " * " +
			                 dimensions (rhs).toString () + " = " + dimensions (*this).toString ());
	}
	void MulMatrixDouble::compute () {
		callWithValues (*this, [](MatrixDouble & r, const MatrixDouble & lhs,
		                          const MatrixDouble & rhs) { r.noalias () = lhs * rhs; });
	}
	std::string MulMatrixDouble::debugInfo () const { return debugInfoFor (this->accessValue ()); }
	NodeRef MulMatrixDouble::derive (const Node & node) {
		auto dim = dimensions (*this);
		auto & lhs = this->dependencies ()[0];
		auto & rhs = this->dependencies ()[1];
		auto lhsDerived = MulMatrixDouble::create (NodeRefVec{lhs->derive (node), rhs}, dim);
		auto rhsDerived = MulMatrixDouble::create (NodeRefVec{lhs, rhs->derive (node)}, dim);
		return AddMatrixDouble::create (NodeRefVec{std::move (lhsDerived), std::move (rhsDerived)},
		                                dim);
	}
	ValueRef<MatrixDouble> MulMatrixDouble::create (NodeRefVec && deps, MatrixDimension dim) {
		checkDependencies<Dependencies> (deps, typeid (MulMatrixDouble));
		auto isConstantZeroMatrix = predicateIsConstantValueMatching<MatrixDouble> (isExactZeroMatrix);
		auto isConstantIdentityMatrix =
		    predicateIsConstantValueMatching<MatrixDouble> (isExactIdentityMatrix);
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
	std::string CWiseMulMatrixDouble::debugInfo () const {
		return debugInfoFor (this->accessValue ());
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
		if (std::any_of (deps.begin (), deps.end (),
		                 predicateIsConstantValueMatching<MatrixDouble> (isExactZeroMatrix))) {
			return ConstantMatrixDouble::createZero (dim);
		}
		// Remove 1s
		removeDependenciesIf (deps, predicateIsConstantValueMatching<MatrixDouble> (isExactOnesMatrix));
		// Select node
		if (deps.size () == 1) {
			return convertRef<Value<MatrixDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return ConstantMatrixDouble::create (onesMatrix (dim));
		} else {
			return std::make_shared<CWiseMulMatrixDouble> (std::move (deps), dim);
		}
	}
} // namespace DF
} // namespace bpp
