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
#include <Bpp/NewPhyl/Utils.h>
#include <algorithm>
#include <memory>
#include <string>
#include <typeinfo>
#include <utility>

namespace bpp {
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
			auto zero = linearAlgebraZeroValue (dimension (t));
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
		return "targetDim=" + to_string (targetDimension (v)) + " dim=" + to_string (dimension (v)) +
		       " " + numericProps (v);
	}
	template <> std::string Value<MatrixDouble>::debugInfo () const {
		using std::to_string;
		auto & m = this->accessValueConst ();
		return "targetDim" + to_string (targetDimension (m)) + " dim" + to_string (dimension (m)) +
		       " " + numericProps (m);
	}

	// MulMatrixDouble
	class MulMatrixDouble : public Value<MatrixDouble> {
	public:
		using Dependencies = TupleOfValues<MatrixDouble, MatrixDouble>;

		MulMatrixDouble (NodeRefVec && deps, const Dimension<MatrixDouble> & dim)
		    : Value<MatrixDouble> (std::move (deps)) {
			setTargetDimension (this->accessValueMutable (), dim);
		}
		NodeRef derive (const Node & node) final {
			auto dim = targetDimension (this->accessValueConst ());
			auto & lhs = this->dependency (0);
			auto & rhs = this->dependency (1);
			auto dLhs = makeNode<MulMatrixDouble> ({lhs->derive (node), rhs}, dim);
			auto dRhs = makeNode<MulMatrixDouble> ({lhs, rhs->derive (node)}, dim);
			return makeNode<AddMatrixDouble> ({std::move (dLhs), std::move (dRhs)}, dim);
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<MulMatrixDouble> (std::move (deps),
			                                  targetDimension (this->accessValueConst ()));
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

	// MulTransposedMatrixVectorDouble
	class MulTransposedMatrixVectorDouble : public Value<VectorDouble> {
	public:
		using Dependencies = TupleOfValues<MatrixDouble, VectorDouble>;

		MulTransposedMatrixVectorDouble (NodeRefVec && deps, const Dimension<VectorDouble> & dim)
		    : Value<VectorDouble> (std::move (deps)) {
			setTargetDimension (this->accessValueMutable (), dim);
		}
		NodeRef derive (const Node & node) final {
			auto dim = targetDimension (this->accessValueConst ());
			auto & lhs = this->dependency (0);
			auto & rhs = this->dependency (1);
			auto dLhs = makeNode<MulTransposedMatrixVectorDouble> ({lhs->derive (node), rhs}, dim);
			auto dRhs = makeNode<MulTransposedMatrixVectorDouble> ({lhs, rhs->derive (node)}, dim);
			return makeNode<AddVectorDouble> ({std::move (dLhs), std::move (dRhs)}, dim);
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<MulTransposedMatrixVectorDouble> (
			    std::move (deps), targetDimension (this->accessValueConst ()));
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

} // namespace DF
} // namespace bpp
