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

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowMatrix.h>
#include <Bpp/NewPhyl/DataFlowTemplateUtils.h>
#include <Eigen/Core>
#include <typeinfo>

namespace bpp {
namespace DF {
	// Explicit template instantiation TODO test usefulness
	// template class Value<Eigen::MatrixXd>;

	namespace {
		bool isExactZeroMatrix (const MatrixDouble & m) {
			return m == MatrixDouble::Zero (m.rows (), m.cols ());
		}
		auto zeroMatrix (SizeType nbRows, SizeType nbCols)
		    -> decltype (MatrixDouble::Zero (nbRows, nbCols)) {
			return MatrixDouble::Zero (nbRows, nbCols);
		}
		auto zeroMatrix (const MatrixDouble & m) -> decltype (zeroMatrix (m.rows (), m.cols ())) {
			return zeroMatrix (m.rows (), m.cols ());
		}
	} // namespace

	// ConstantMatrixDouble
	void ConstantMatrixDouble::compute () { failureComputeWasCalled (typeid (ConstantMatrixDouble)); }
	NodeRef ConstantMatrixDouble::derive (const Node &) {
		return ConstantMatrixDouble::create (zeroMatrix (this->accessValue ()));
	}

	// AddMatrixDouble
	AddMatrixDouble::AddMatrixDouble (NodeRefVec && deps, SizeType nbRows, SizeType nbCols)
	    : Value<MatrixDouble> (std::move (deps), nbRows, nbCols) {
		checkDependencies (*this);
    // TODO check dimensions
	}
	void AddMatrixDouble::compute () {
		callWithValues (*this, [](MatrixDouble & r) { r.fill (0.); },
		                [](MatrixDouble & r, const MatrixDouble & m) { r += m; });
	}
	NodeRef AddMatrixDouble::derive (const Node & node) {
		const auto & m = this->accessValue ();
		return AddMatrixDouble::create (
		    this->dependencies ().map ([&node](const NodeRef & dep) { return dep->derive (node); }),
		    m.rows (), m.cols ());
	}
	ValueRef<MatrixDouble> AddMatrixDouble::create (NodeRefVec && deps, SizeType nbRows,
	                                                SizeType nbCols) {
		// Remove Os
		removeDependenciesIf (deps, predicateIsConstantValueMatching<MatrixDouble> (isExactZeroMatrix));
		// Select node impl
		if (deps.size () == 1) {
			return convertRef<Value<MatrixDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return ConstantMatrixDouble::create (zeroMatrix (nbRows, nbCols));
		} else {
			return std::make_shared<AddMatrixDouble> (std::move (deps), nbRows, nbCols);
		}
	}

	// MulMatrixDouble
	MulMatrixDouble::MulMatrixDouble (NodeRefVec && deps, SizeType nbRows, SizeType nbCols)
	    : Value<MatrixDouble> (std::move (deps), nbRows, nbCols) {
		checkDependencies (*this);
    // TODO check dimensions
	}
	void MulMatrixDouble::compute () {
		callWithValues (*this, [](MatrixDouble & r, const MatrixDouble & lhs,
		                          const MatrixDouble & rhs) { r = lhs * rhs; });
	}
	NodeRef MulMatrixDouble::derive (const Node & node) {
		auto rows = this->accessValue ().rows ();
		auto cols = this->accessValue ().cols ();
		auto & lhs = this->dependencies ()[0];
		auto & rhs = this->dependencies ()[1];
		auto lhsDerived = MulMatrixDouble::create (NodeRefVec{lhs->derive (node), rhs}, rows, cols);
		auto rhsDerived = MulMatrixDouble::create (NodeRefVec{lhs, rhs->derive (node)}, rows, cols);
		return AddMatrixDouble::create (NodeRefVec{std::move (lhsDerived), std::move (rhsDerived)},
		                                rows, cols);
	}
	ValueRef<MatrixDouble> MulMatrixDouble::create (NodeRefVec && deps, SizeType nbRows,
	                                                SizeType nbCols) {
		// Return O if any matrix is 0
		if (std::any_of (deps.begin (), deps.end (),
		                 predicateIsConstantValueMatching<MatrixDouble> (isExactZeroMatrix))) {
			return ConstantMatrixDouble::create (zeroMatrix (nbRows, nbCols));
		} else {
			return std::make_shared<MulMatrixDouble> (std::move (deps), nbRows, nbCols);
		}
    // TODO improve
    // only 2 args -> check deps beforehand
    // if one of them is 0, return 0
    // if one of them is identity, return other
    // default : return mul
	}

} // namespace DF
} // namespace bpp
