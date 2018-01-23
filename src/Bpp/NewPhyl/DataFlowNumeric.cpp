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
#include <Bpp/NewPhyl/DataFlowInternal.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/NewPhyl/Debug.h> // checks
#include <Bpp/NewPhyl/IntegerRange.h>
#include <Bpp/NewPhyl/LinearAlgebra.h>
#include <algorithm>
#include <memory>
#include <typeinfo>
#include <utility>

namespace bpp {
/* TODO
 * - opt: aggregate constants
 * - merge constants by type
 * - support for general register (reuse code in DF.cpp, declare in DFTinternal).
 * - opt: keep 0s instead of recreating them
 * - deduce size (require 1 arg in reductions) ?
 * - opt: merge * of *, + of +, ... ?
 * - common template for *,+ opts ?
 */

// TODO fwd decl: cannot convert sp<Value<T>> to sp<Node>: provide make(Value<T>, ...) ?

namespace DF {
	/******************************** Utils *******************************/
	namespace {
		// Constant node checking predicates (const NodeRef & -> bool)
		auto isConstantZeroDouble = predicateIsConstantValueMatching<double> (isExactZero);
		auto isConstantOneDouble = predicateIsConstantValueMatching<double> (isExactOne);
		auto isConstantZeroVector = predicateIsConstantValueMatching<VectorDouble> (isExactZero);
		auto isConstantOnesVector = predicateIsConstantValueMatching<VectorDouble> (isExactOne);
		auto isConstantZeroMatrix = predicateIsConstantValueMatching<MatrixDouble> (isExactZero);
		auto isConstantOnesMatrix = predicateIsConstantValueMatching<MatrixDouble> (isExactOne);
		auto isConstantIdentityMatrix =
		    predicateIsConstantValueMatching<MatrixDouble> (isExactIdentity);

	} // namespace

	bool derivableIfAllDepsAre (const Node & toDerive, const Node & node) {
		auto & deps = toDerive.dependencies ();
		return std::all_of (deps.begin (), deps.end (),
		                    [&node](const NodeRef & dep) { return dep->isDerivable (node); });
	}

	/******************************** New Nodes *******************************/
	// TODO templatize impls

	/******************************** Old nodes *******************************/

	// AddDouble
	class AddDouble : public Value<double> {
	public:
		using Dependencies = ReductionOfValue<double>;

		AddDouble (NodeRefVec && deps) : Value<double> (std::move (deps)) {}
		NodeRef derive (const Node & node) final {
			return makeNode<AddDouble> (mapToVector (
			    this->dependencies (), [&node](const NodeRef & dep) { return dep->derive (node); }));
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<AddDouble> (std::move (deps));
		}

	private:
		void compute () final {
			callWithValues (*this, [](double & r) { r = 0.; }, [](double & r, double d) { r += d; });
		}
	};
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
	class MulDouble : public Value<double> {
	public:
		using Dependencies = ReductionOfValue<double>;

		MulDouble (NodeRefVec && deps) : Value<double> (std::move (deps)) {}
		NodeRef derive (const Node & node) final {
			NodeRefVec addDeps;
			for (auto i : bpp::range (this->nbDependencies ())) {
				NodeRefVec mulDeps = this->dependencies ();
				mulDeps[i] = this->dependency (i)->derive (node);
				addDeps.emplace_back (makeNode<MulDouble> (std::move (mulDeps)));
			}
			return makeNode<AddDouble> (std::move (addDeps));
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<MulDouble> (std::move (deps));
		}

	private:
		void compute () final {
			callWithValues (*this, [](double & r) { r = 1.; }, [](double & r, double d) { r *= d; });
		}
	};
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

	// NegDouble
	class NegDouble : public Value<double> {
	public:
		using Dependencies = FunctionOfValues<double>;

		NegDouble (NodeRefVec && deps) : Value<double> (std::move (deps)) {}
		NodeRef derive (const Node & node) final {
			return makeNode<NegDouble> ({this->dependency (0)->derive (node)});
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<NegDouble> (std::move (deps));
		}

	private:
		void compute () final {
			callWithValues (*this, [](double & r, double d) { r = -d; });
		}
	};
	ValueRef<double> Builder<NegDouble>::make (NodeRefVec && deps) {
		checkDependencies<NegDouble> (deps);
		auto & dep = deps[0];
		if (dep->isConstant ()) {
			return makeNode<Constant<double>> (-accessValidValueConstCast<double> (dep));
		} else {
			return std::make_shared<NegDouble> (std::move (deps));
		}
	}

	// ScalarProdDouble
	class ScalarProdDouble : public Value<double> {
	public:
		using Dependencies = FunctionOfValues<VectorDouble, VectorDouble>;

		ScalarProdDouble (NodeRefVec && deps) : Value<double> (std::move (deps)) {}
		NodeRef derive (const Node & node) final {
			auto & lhs = this->dependency (0);
			auto & rhs = this->dependency (1);
			auto dLhs = makeNode<ScalarProdDouble> ({lhs->derive (node), rhs});
			auto dRhs = makeNode<ScalarProdDouble> ({lhs, rhs->derive (node)});
			return makeNode<AddDouble> ({std::move (dLhs), std::move (dRhs)});
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<ScalarProdDouble> (std::move (deps));
		}

	private:
		void compute () final {
			callWithValues (*this, [](double & r, const VectorDouble & lhs, const VectorDouble & rhs) {
				r = lhs.dot (rhs);
			});
		}
	};
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
	class AddVectorDouble : public Value<VectorDouble> {
	public:
		using Dependencies = ReductionOfValue<VectorDouble>;

		AddVectorDouble (NodeRefVec && deps, const Dimension<VectorDouble> & dim)
		    : Value<VectorDouble> (std::move (deps), dim.size) {}
		NodeRef derive (const Node & node) final {
			return makeNode<AddVectorDouble> (
			    mapToVector (this->dependencies (),
			                 [&node](const NodeRef & dep) { return dep->derive (node); }),
			    dimensions (*this));
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<AddVectorDouble> (std::move (deps), dimensions (*this));
		}

	private:
		void compute () final {
			callWithValues (*this, [](VectorDouble & r) { r.setZero (); },
			                [](VectorDouble & r, const VectorDouble & v) { r += v; });
		}
	};
	ValueRef<VectorDouble> Builder<AddVectorDouble>::make (NodeRefVec && deps,
	                                                       const Dimension<VectorDouble> & dim) {
		checkDependencies<AddVectorDouble> (deps);
		// Remove Os
		removeDependenciesIf (deps, isConstantZeroVector);
		// Select node impl
		if (deps.size () == 1) {
			return convertRef<Value<VectorDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return Builder<Constant<VectorDouble>>::makeZero (dim);
		} else {
			return std::make_shared<AddVectorDouble> (std::move (deps), dim);
		}
	}

	// CWiseMulVectorDouble
	class CWiseMulVectorDouble : public Value<VectorDouble> {
	public:
		using Dependencies = ReductionOfValue<VectorDouble>;

		CWiseMulVectorDouble (NodeRefVec && deps, const Dimension<VectorDouble> & dim)
		    : Value<VectorDouble> (std::move (deps), dim.size) {}
		NodeRef derive (const Node & node) final {
			auto dim = dimensions (*this);
			NodeRefVec addDeps;
			for (auto i : bpp::range (this->nbDependencies ())) {
				NodeRefVec mulDeps = this->dependencies ();
				mulDeps[i] = this->dependency (i)->derive (node);
				addDeps.emplace_back (makeNode<CWiseMulVectorDouble> (std::move (mulDeps), dim));
			}
			return makeNode<AddVectorDouble> (std::move (addDeps), dim);
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<CWiseMulVectorDouble> (std::move (deps), dimensions (*this));
		}

	private:
		void compute () final {
			callWithValues (*this, [](VectorDouble & r) { r.setOnes (); },
			                [](VectorDouble & r, const VectorDouble & v) { r = r.cwiseProduct (v); });
		}
	};
	ValueRef<VectorDouble> Builder<CWiseMulVectorDouble>::make (NodeRefVec && deps,
	                                                            const Dimension<VectorDouble> & dim) {
		checkDependencies<CWiseMulVectorDouble> (deps);
		// Return 0 if any 0 dep
		if (std::any_of (deps.begin (), deps.end (), isConstantZeroVector)) {
			return Builder<Constant<VectorDouble>>::makeZero (dim);
		}
		// Remove 1s
		removeDependenciesIf (deps, isConstantOnesVector);
		// Select node
		if (deps.size () == 1) {
			return convertRef<Value<VectorDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return Builder<Constant<VectorDouble>>::makeOne (dim);
		} else {
			return std::make_shared<CWiseMulVectorDouble> (std::move (deps), dim);
		}
	}

	// CWiseNegVectorDouble
	class CWiseNegVectorDouble : public Value<VectorDouble> {
	public:
		using Dependencies = FunctionOfValues<VectorDouble>;

		CWiseNegVectorDouble (NodeRefVec && deps, const Dimension<VectorDouble> & dim)
		    : Value<VectorDouble> (std::move (deps), dim.size) {}
		NodeRef derive (const Node & node) final {
			return makeNode<CWiseNegVectorDouble> ({this->dependency (0)->derive (node)},
			                                       dimensions (*this));
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<CWiseNegVectorDouble> (std::move (deps), dimensions (*this));
		}

	private:
		void compute () final {
			callWithValues (*this, [](VectorDouble & r, const VectorDouble & v) { r.noalias () = -v; });
		}
	};
	ValueRef<VectorDouble> Builder<CWiseNegVectorDouble>::make (NodeRefVec && deps,
	                                                            const Dimension<VectorDouble> & dim) {
		checkDependencies<CWiseNegVectorDouble> (deps);
		auto & dep = deps[0];
		if (dep->isConstant ()) {
			return makeNode<Constant<VectorDouble>> (-accessValidValueConstCast<VectorDouble> (dep));
		} else {
			return std::make_shared<CWiseNegVectorDouble> (std::move (deps), dim);
		}
	}

	// CWiseInverseVectorDouble
	class CWiseInverseVectorDouble : public Value<VectorDouble> {
	public:
		using Dependencies = FunctionOfValues<VectorDouble>;

		CWiseInverseVectorDouble (NodeRefVec && deps, const Dimension<VectorDouble> & dim)
		    : Value<VectorDouble> (std::move (deps), dim.size) {}
		NodeRef derive (const Node & node) final {
			auto dim = dimensions (*this);
			auto & arg = this->dependency (0);
			return makeNode<CWiseNegVectorDouble> (
			    {makeNode<CWiseConstantPowVectorDouble> ({arg}, dim, -2.), arg->derive (node)}, dim);
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<CWiseInverseVectorDouble> (std::move (deps), dimensions (*this));
		}

	private:
		void compute () final {
			callWithValues (*this, [](VectorDouble & r, const VectorDouble & v) {
				r.noalias () = v.cwiseInverse ();
			});
		}
	};
	ValueRef<VectorDouble>
	Builder<CWiseInverseVectorDouble>::make (NodeRefVec && deps,
	                                         const Dimension<VectorDouble> & dim) {
		checkDependencies<CWiseInverseVectorDouble> (deps);
		auto & arg = deps[0];
		if (isConstantOnesVector (arg)) {
			return convertRef<Value<VectorDouble>> (arg);
		} else {
			return std::make_shared<CWiseInverseVectorDouble> (std::move (deps), dim);
		}
	}

	// CWiseConstantPowVectorDouble
	class CWiseConstantPowVectorDouble : public Value<VectorDouble> {
	public:
		using Dependencies = FunctionOfValues<VectorDouble>;

		CWiseConstantPowVectorDouble (NodeRefVec && deps, const Dimension<VectorDouble> & dim,
		                              double exp)
		    : Value<VectorDouble> (std::move (deps), dim.size), exp_ (exp) {}
		NodeRef derive (const Node & node) final {
			auto dim = dimensions (*this);
			auto & arg = this->dependency (0);
			auto powDerivative = makeNode<CWiseMulScalarVectorDouble> (
			    {makeNode<Constant<double>> (exp_),
			     makeNode<CWiseConstantPowVectorDouble> ({arg}, dim, exp_ - 1.)},
			    dim);
			return makeNode<CWiseMulVectorDouble> ({std::move (powDerivative), arg->derive (node)}, dim);
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<CWiseConstantPowVectorDouble> (std::move (deps), dimensions (*this), exp_);
		}

	private:
		void compute () final {
			callWithValues (*this, [this](VectorDouble & r, const VectorDouble & v) {
				r.noalias () = v.array ().pow (exp_).matrix ();
			});
		}
		double exp_;
	};
	ValueRef<VectorDouble>
	Builder<CWiseConstantPowVectorDouble>::make (NodeRefVec && deps,
	                                             const Dimension<VectorDouble> & dim, double exp) {
		checkDependencies<CWiseConstantPowVectorDouble> (deps);
		auto & arg = deps[0];
		if (exp == 1.) {
			return convertRef<Value<VectorDouble>> (arg);
		} else if (exp == 0.) {
			return Builder<Constant<VectorDouble>>::makeOne (dim);
		} else if (exp == -1.) {
			return makeNode<CWiseInverseVectorDouble> (std::move (deps), dim);
		} else if (isConstantOnesVector (arg)) {
			return convertRef<Value<VectorDouble>> (arg);
		} else {
			return std::make_shared<CWiseConstantPowVectorDouble> (std::move (deps), dim, exp);
		}
	}

	// AddMatrixDouble
	class AddMatrixDouble : public Value<MatrixDouble> {
	public:
		using Dependencies = ReductionOfValue<MatrixDouble>;

		AddMatrixDouble (NodeRefVec && deps, const Dimension<MatrixDouble> & dim)
		    : Value<MatrixDouble> (std::move (deps), dim.rows, dim.cols) {}
		NodeRef derive (const Node & node) final {
			return makeNode<AddMatrixDouble> (
			    mapToVector (this->dependencies (),
			                 [&node](const NodeRef & dep) { return dep->derive (node); }),
			    dimensions (*this));
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<AddMatrixDouble> (std::move (deps), dimensions (*this));
		}

	private:
		void compute () final {
			callWithValues (*this, [](MatrixDouble & r) { r.setZero (); },
			                [](MatrixDouble & r, const MatrixDouble & m) { r += m; });
		}
	};
	ValueRef<MatrixDouble> Builder<AddMatrixDouble>::make (NodeRefVec && deps,
	                                                       const Dimension<MatrixDouble> & dim) {
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
	class MulMatrixDouble : public Value<MatrixDouble> {
	public:
		using Dependencies = FunctionOfValues<MatrixDouble, MatrixDouble>;

		MulMatrixDouble (NodeRefVec && deps, const Dimension<MatrixDouble> & dim)
		    : Value<MatrixDouble> (std::move (deps), dim.rows, dim.cols) {}
		NodeRef derive (const Node & node) final {
			auto dim = dimensions (*this);
			auto & lhs = this->dependency (0);
			auto & rhs = this->dependency (1);
			auto dLhs = makeNode<MulMatrixDouble> ({lhs->derive (node), rhs}, dim);
			auto dRhs = makeNode<MulMatrixDouble> ({lhs, rhs->derive (node)}, dim);
			return makeNode<AddMatrixDouble> ({std::move (dLhs), std::move (dRhs)}, dim);
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<MulMatrixDouble> (std::move (deps), dimensions (*this));
		}

	private:
		void compute () final {
			callWithValues (*this, [](MatrixDouble & r, const MatrixDouble & lhs,
			                          const MatrixDouble & rhs) { r.noalias () = lhs * rhs; });
		}
	};
	ValueRef<MatrixDouble> Builder<MulMatrixDouble>::make (NodeRefVec && deps,
	                                                       const Dimension<MatrixDouble> & dim) {
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
	class CWiseMulMatrixDouble : public Value<MatrixDouble> {
	public:
		using Dependencies = ReductionOfValue<MatrixDouble>;

		CWiseMulMatrixDouble (NodeRefVec && deps, const Dimension<MatrixDouble> & dim)
		    : Value<MatrixDouble> (std::move (deps), dim.rows, dim.cols) {}
		NodeRef derive (const Node & node) final {
			auto dim = dimensions (*this);
			NodeRefVec addDeps;
			for (auto i : bpp::range (this->nbDependencies ())) {
				NodeRefVec mulDeps = this->dependencies ();
				mulDeps[i] = this->dependency (i)->derive (node);
				addDeps.emplace_back (makeNode<CWiseMulMatrixDouble> (std::move (mulDeps), dim));
			}
			return makeNode<AddMatrixDouble> (std::move (addDeps), dim);
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<CWiseMulMatrixDouble> (std::move (deps), dimensions (*this));
		}

	private:
		void compute () final {
			callWithValues (*this, [](MatrixDouble & r) { r.setOnes (); },
			                [](MatrixDouble & r, const MatrixDouble & m) { r = r.cwiseProduct (m); });
		}
	};
	ValueRef<MatrixDouble> Builder<CWiseMulMatrixDouble>::make (NodeRefVec && deps,
	                                                            const Dimension<MatrixDouble> & dim) {
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
	class MulTransposedMatrixVectorDouble : public Value<VectorDouble> {
	public:
		using Dependencies = FunctionOfValues<MatrixDouble, VectorDouble>;

		MulTransposedMatrixVectorDouble (NodeRefVec && deps, const Dimension<VectorDouble> & dim)
		    : Value<VectorDouble> (std::move (deps), dim.size) {}
		NodeRef derive (const Node & node) final {
			auto dim = dimensions (*this);
			auto & lhs = this->dependency (0);
			auto & rhs = this->dependency (1);
			auto dLhs = makeNode<MulTransposedMatrixVectorDouble> ({lhs->derive (node), rhs}, dim);
			auto dRhs = makeNode<MulTransposedMatrixVectorDouble> ({lhs, rhs->derive (node)}, dim);
			return makeNode<AddVectorDouble> ({std::move (dLhs), std::move (dRhs)}, dim);
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<MulTransposedMatrixVectorDouble> (std::move (deps), dimensions (*this));
		}

	private:
		void compute () final {
			callWithValues (*this,
			                [](VectorDouble & r, const MatrixDouble & lhs, const VectorDouble & rhs) {
				                r.noalias () = lhs.transpose () * rhs;
			                });
		}
	};
	ValueRef<VectorDouble>
	Builder<MulTransposedMatrixVectorDouble>::make (NodeRefVec && deps,
	                                                const Dimension<VectorDouble> & dim) {
		checkDependencies<MulTransposedMatrixVectorDouble> (deps);
		auto & lhs = deps[0];
		auto & rhs = deps[1];
		if (isConstantZeroMatrix (lhs) || isConstantZeroVector (rhs)) {
			return Builder<Constant<VectorDouble>>::makeZero (dim);
		} else if (isConstantIdentityMatrix (lhs)) {
			return convertRef<Value<VectorDouble>> (rhs);
		} else {
			return std::make_shared<MulTransposedMatrixVectorDouble> (std::move (deps), dim);
		}
	}

	// CWiseMulScalarVectorDouble
	class CWiseMulScalarVectorDouble : public Value<VectorDouble> {
	public:
		using Dependencies = FunctionOfValues<double, VectorDouble>;

		CWiseMulScalarVectorDouble (NodeRefVec && deps, const Dimension<VectorDouble> & dim)
		    : Value<VectorDouble> (std::move (deps), dim.size) {}
		NodeRef derive (const Node & node) final {
			auto dim = dimensions (*this);
			auto & lhs = this->dependency (0);
			auto & rhs = this->dependency (1);
			auto dLhs = makeNode<CWiseMulScalarVectorDouble> ({lhs->derive (node), rhs}, dim);
			auto dRhs = makeNode<CWiseMulScalarVectorDouble> ({lhs, rhs->derive (node)}, dim);
			return makeNode<AddVectorDouble> ({std::move (dLhs), std::move (dRhs)}, dim);
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<CWiseMulScalarVectorDouble> (std::move (deps), dimensions (*this));
		}

	private:
		void compute () final {
			callWithValues (
			    *this, [](VectorDouble & r, double d, const VectorDouble & v) { r.noalias () = d * v; });
		}
	};
	ValueRef<VectorDouble>
	Builder<CWiseMulScalarVectorDouble>::make (NodeRefVec && deps,
	                                           const Dimension<VectorDouble> & dim) {
		checkDependencies<CWiseMulScalarVectorDouble> (deps);
		auto & lhs = deps[0];
		auto & rhs = deps[1];
		if (isConstantZeroDouble (lhs) || isConstantZeroVector (rhs)) {
			return Builder<Constant<VectorDouble>>::makeZero (dim);
		} else if (isConstantOneDouble (lhs)) {
			return convertRef<Value<VectorDouble>> (rhs);
		} else {
			return std::make_shared<CWiseMulScalarVectorDouble> (std::move (deps), dim);
		}
	}

	// CWiseMulScalarMatrixDouble
	class CWiseMulScalarMatrixDouble : public Value<MatrixDouble> {
	public:
		using Dependencies = FunctionOfValues<double, MatrixDouble>;

		CWiseMulScalarMatrixDouble (NodeRefVec && deps, const Dimension<MatrixDouble> & dim)
		    : Value<MatrixDouble> (std::move (deps), dim.rows, dim.cols) {}
		NodeRef derive (const Node & node) final {
			auto dim = dimensions (*this);
			auto & lhs = this->dependency (0);
			auto & rhs = this->dependency (1);
			auto dLhs = makeNode<CWiseMulScalarMatrixDouble> ({lhs->derive (node), rhs}, dim);
			auto dRhs = makeNode<CWiseMulScalarMatrixDouble> ({lhs, rhs->derive (node)}, dim);
			return makeNode<AddMatrixDouble> ({std::move (dLhs), std::move (dRhs)}, dim);
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<CWiseMulScalarMatrixDouble> (std::move (deps), dimensions (*this));
		}

	private:
		void compute () final {
			callWithValues (
			    *this, [](MatrixDouble & r, double d, const MatrixDouble & m) { r.noalias () = d * m; });
		}
	};
	ValueRef<MatrixDouble>
	Builder<CWiseMulScalarMatrixDouble>::make (NodeRefVec && deps,
	                                           const Dimension<MatrixDouble> & dim) {
		checkDependencies<CWiseMulScalarMatrixDouble> (deps);
		auto & lhs = deps[0];
		auto & rhs = deps[1];
		if (isConstantZeroDouble (lhs) || isConstantZeroMatrix (rhs)) {
			return Builder<Constant<MatrixDouble>>::makeZero (dim);
		} else if (isConstantOneDouble (lhs)) {
			return convertRef<Value<MatrixDouble>> (rhs);
		} else {
			return std::make_shared<CWiseMulScalarMatrixDouble> (std::move (deps), dim);
		}
	}
} // namespace DF
} // namespace bpp
