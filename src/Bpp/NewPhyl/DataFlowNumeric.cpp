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

#include <Bpp/NewPhyl/DataFlowInternal.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/NewPhyl/DataFlowTemplates.h>
#include <Bpp/NewPhyl/Debug.h> // checks
#include <Bpp/NewPhyl/IntegerRange.h>
#include <Bpp/NewPhyl/LinearAlgebra.h>
#include <Bpp/NewPhyl/LinearAlgebraUtils.h>
#include <algorithm>
#include <memory>
#include <string>
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
		auto isConstantZeroNode = [](const NodeRef & n) { return n->isConstantZero (); };
		auto isConstantOneNode = [](const NodeRef & n) { return n->isConstantOne (); };

		auto isConstantIdentityMatrix = [](const NodeRef & n) {
			// FIXME assume Constant<MatrixDouble>, not lazily evaluated
			return n->isConstant () &&
			       isExactIdentity (nodeValueCast<MatrixDouble> (*n).accessValueConst ());
		};

		/* Generate a string describing useful numeric props of Eigen vector/matrix.
		 * Used in Value<T>::debugInfo specialisations.
		 */
		template <typename T> std::string numericProps (const T & t) {
			std::string s{"props{"};
			auto zero = linearAlgebraZeroValue (dimensions (t));
			// Property for all elements
			if (t.isZero (0.))
				s += "[0]";
			if (t.isOnes (0.))
				s += "[1]";
			if (t.rows () == t.cols () && t.isIdentity (0.))
				s += "[I]";
			// Property on any element
			if (t.array ().isNaN ().any ())
				s += "N";
			if (t.array ().isInf ().any ())
				s += "i";
			if ((t.array () == zero.array ()).any ())
				s += "0";
			if ((t.array () > zero.array ()).any ())
				s += "+";
			if ((t.array () < zero.array ()).any ())
				s += "-";
			s += "}";
			return s;
		}
	} // namespace

	bool derivableIfAllDepsAre (const Node & toDerive, const Node & node) {
		auto & deps = toDerive.dependencies ();
		return std::all_of (deps.begin (), deps.end (),
		                    [&node](const NodeRef & dep) { return dep->isDerivable (node); });
	}

	/**************************************************************************
	 * Specialisations of Value/Constant/Mutable for numeric types.
	 */

	// Debug info of Value<T>
	template <> std::string Value<double>::debugInfo () const {
		return std::to_string (this->accessValueConst ());
	}
	template <> std::string Value<VectorDouble>::debugInfo () const {
		using std::to_string;
		auto & v = this->accessValueConst ();
		return "targetDim=" + to_string (this->getTargetDimension ()) +
		       " dim=" + to_string (dimensions (v)) + " " + numericProps (v);
	}
	template <> std::string Value<MatrixDouble>::debugInfo () const {
		using std::to_string;
		auto & m = this->accessValueConst ();
		return "targetDim" + to_string (this->getTargetDimension ()) + " dim" +
		       to_string (dimensions (m)) + " " + numericProps (m);
	}

	// Mutable<double> specialisation
	template <> NodeRef Mutable<double>::derive (const Node & node) {
		if (&node == static_cast<const Node *> (this)) {
			return makeNode<ConstantOne<double>> ();
		} else {
			return makeNode<ConstantZero<double>> ();
		}
	}
	template <> bool Mutable<double>::isDerivable (const Node &) const { return true; }

	// Constant<T> specialisation
	template <> NodeRef Constant<double>::derive (const Node &) {
		return makeNode<ConstantZero<double>> ();
	}
	template <> bool Constant<double>::isDerivable (const Node &) const { return true; }

	template <> NodeRef Constant<VectorDouble>::derive (const Node &) {
		return makeNode<ConstantZero<VectorDouble>> (dimensions (this->accessValueConst ()));
	}
	template <> bool Constant<VectorDouble>::isDerivable (const Node &) const { return true; }

	template <> NodeRef Constant<MatrixDouble>::derive (const Node &) {
		return makeNode<ConstantZero<MatrixDouble>> (dimensions (this->accessValueConst ()));
	}
	template <> bool Constant<MatrixDouble>::isDerivable (const Node &) const { return true; }

	/***********************************************************************
	 * New style numeric nodes.
	 * Nodes are defined as template classes and template construction functions.
	 * This reduce code duplication and simplify changes.
	 * These class templates are only visible in this file.
	 * Builder<NodeType>::make functions are provided for explicit types.
	 */

	// ConstantZero
	template <typename T> class ConstantZero : public Value<T> {
	public:
		explicit ConstantZero (const Dimension<T> & dim) : Value<T> (noDependency) {
			this->setTargetDimension (dim);
		}

		bool isConstant () const final { return true; }
		bool isConstantZero () const final { return true; }

		NodeRef derive (const Node &) final {
			return makeNode<ConstantZero<T>> (this->getTargetDimension ());
		}
		bool isDerivable (const Node &) const final { return true; }

	private:
		void compute () final {
			this->accessValueMutable () = linearAlgebraZeroValue (this->getTargetDimension ());
		}
	};

	ValueRef<double> Builder<ConstantZero<double>>::make (const Dimension<double> & dim) {
		return std::make_shared<ConstantZero<double>> (dim);
	}
	ValueRef<VectorDouble>
	Builder<ConstantZero<VectorDouble>>::make (const Dimension<VectorDouble> & dim) {
		return std::make_shared<ConstantZero<VectorDouble>> (dim);
	}
	ValueRef<MatrixDouble>
	Builder<ConstantZero<MatrixDouble>>::make (const Dimension<MatrixDouble> & dim) {
		return std::make_shared<ConstantZero<MatrixDouble>> (dim);
	}

	// ConstantOne
	template <typename T> class ConstantOne : public Value<T> {
	public:
		explicit ConstantOne (const Dimension<T> & dim) : Value<T> (noDependency) {
			this->setTargetDimension (dim);
		}

		bool isConstant () const final { return true; }
		bool isConstantOne () const final { return true; }

		NodeRef derive (const Node &) final {
			return makeNode<ConstantZero<T>> (this->getTargetDimension ());
		}
		bool isDerivable (const Node &) const final { return true; }

	private:
		void compute () final {
			this->accessValueMutable () = linearAlgebraOneValue (this->getTargetDimension ());
		}
	};

	ValueRef<double> Builder<ConstantOne<double>>::make (const Dimension<double> & dim) {
		return std::make_shared<ConstantOne<double>> (dim);
	}
	ValueRef<VectorDouble>
	Builder<ConstantOne<VectorDouble>>::make (const Dimension<VectorDouble> & dim) {
		return std::make_shared<ConstantOne<VectorDouble>> (dim);
	}
	ValueRef<MatrixDouble>
	Builder<ConstantOne<MatrixDouble>>::make (const Dimension<MatrixDouble> & dim) {
		return std::make_shared<ConstantOne<MatrixDouble>> (dim);
	}

#if 0
	// TODO
	template <typename Result, typename Dependencies> class CWiseAdd : public Value<Result> {
	public:
		CWiseAdd (NodeRefVec && deps, const Dimension<Result> & targetDim)
		    : Value<Result> (std::move (deps)) {
			this->setTargetDimension (targetDim);
		}

		NodeRef derive (const Node & node) final {
			return makeNode<CWiseAdd<Result, Dependencies>> (
			    mapToVector (this->dependencies (),
			                 [&node](const NodeRef & dep) { return dep->derive (node); }),
			    this->getTargetDimension ());
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }

		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<CWiseAdd<Result, Dependencies>> (std::move (deps),
			                                                 this->getTargetDimension ());
		}

	private:
		void compute () final {} // FIXME tried to use partial spec, fails :/
	};

	ValueRef<double>
	Builder<CWiseAdd<double, ReductionOfValue<double>>>::make (NodeRefVec && deps,
	                                                           const Dimension<double> & dim) {
		return std::make_shared<CWiseAdd<double, ReductionOfValue<double>>> (std::move (deps), dim);
	}
	ValueRef<double>
	Builder<CWiseAdd<double, TupleOfValues<double, double>>>::make (NodeRefVec && deps,
	                                                                const Dimension<double> & dim) {
		return std::make_shared<CWiseAdd<double, TupleOfValues<double, double>>> (std::move (deps),
		                                                                          dim);
	}
#endif

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
		removeDependenciesIf (deps, isConstantZeroNode);
		// Node choice
		if (deps.size () == 1) {
			return convertRef<Value<double>> (deps[0]);
		} else if (deps.size () == 0) {
			return makeNode<ConstantZero<double>> ();
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
		if (std::any_of (deps.begin (), deps.end (), isConstantZeroNode)) {
			return makeNode<ConstantZero<double>> ();
		}
		// Remove any 1s
		removeDependenciesIf (deps, isConstantOneNode);
		// Node choice
		if (deps.size () == 1) {
			return convertRef<Value<double>> (deps[0]);
		} else if (deps.size () == 0) {
			return makeNode<ConstantOne<double>> ();
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
		if (isConstantZeroNode (lhs) || isConstantZeroNode (rhs)) {
			return makeNode<ConstantZero<double>> ();
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
		removeDependenciesIf (deps, isConstantZeroNode);
		// Select node impl
		if (deps.size () == 1) {
			return convertRef<Value<VectorDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return makeNode<ConstantZero<VectorDouble>> (dim);
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
		if (std::any_of (deps.begin (), deps.end (), isConstantZeroNode)) {
			return makeNode<ConstantZero<VectorDouble>> (dim);
		}
		// Remove 1s
		removeDependenciesIf (deps, isConstantOneNode);
		// Select node
		if (deps.size () == 1) {
			return convertRef<Value<VectorDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return makeNode<ConstantOne<VectorDouble>> (dim);
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
		if (isConstantOneNode (arg)) {
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
			return makeNode<ConstantOne<VectorDouble>> (dim);
		} else if (exp == -1.) {
			return makeNode<CWiseInverseVectorDouble> (std::move (deps), dim);
		} else if (isConstantOneNode (arg)) {
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
		removeDependenciesIf (deps, isConstantZeroNode);
		// Select node impl
		if (deps.size () == 1) {
			return convertRef<Value<MatrixDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return makeNode<ConstantZero<MatrixDouble>> (dim);
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
		if (isConstantZeroNode (lhs) || isConstantZeroNode (rhs)) {
			return makeNode<ConstantZero<MatrixDouble>> (dim);
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
		if (std::any_of (deps.begin (), deps.end (), isConstantZeroNode)) {
			return makeNode<ConstantZero<MatrixDouble>> (dim);
		}
		// Remove 1s
		removeDependenciesIf (deps, isConstantOneNode);
		// Select node
		if (deps.size () == 1) {
			return convertRef<Value<MatrixDouble>> (deps[0]);
		} else if (deps.size () == 0) {
			return makeNode<ConstantOne<MatrixDouble>> (dim);
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
		if (isConstantZeroNode (lhs) || isConstantZeroNode (rhs)) {
			return makeNode<ConstantZero<VectorDouble>> (dim);
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
		if (isConstantZeroNode (lhs) || isConstantZeroNode (rhs)) {
			return makeNode<ConstantZero<VectorDouble>> (dim);
		} else if (isConstantOneNode (lhs)) {
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
		if (isConstantZeroNode (lhs) || isConstantZeroNode (rhs)) {
			return makeNode<ConstantZero<MatrixDouble>> (dim);
		} else if (isConstantOneNode (lhs)) {
			return convertRef<Value<MatrixDouble>> (rhs);
		} else {
			return std::make_shared<CWiseMulScalarMatrixDouble> (std::move (deps), dim);
		}
	}
} // namespace DF
} // namespace bpp
