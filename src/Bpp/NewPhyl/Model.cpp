//
// File: Model.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-29
// Last modified: 2017-05-29
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

#include <Bpp/Exceptions.h>
#include <Bpp/NewPhyl/DataFlowInternal.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/NewPhyl/IntegerRange.h>
#include <Bpp/NewPhyl/LinearAlgebra.h>
#include <Bpp/NewPhyl/Model.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <algorithm>
#include <cassert>

namespace bpp {
namespace DF {

	/* Model values derivation:
	 *
	 * We only have analytical derivation for transition matrix with respect to brlen.
	 * A model is defined as "derivable(x)" if its parameters do not depend on "x".
	 *
	 * With a "derivable(x)" model:
	 * - equilibrium frequencies are a constant (derive to 0s).
	 * - transition matrices are only dependent on brlen (which may depend on "x", constant if not).
	 *
	 * derive(x) methods for model value compute nodes assume the model is derivable(x).
	 * derive(x) will not fail is not the case, but the derivative will be wrong.
	 * This is checked by an assert in debug mode.
	 * derive(x) for these nodes should not be called if isDerivable(x) is false.
	 * FIXME check anyway at derive(x), throw exception ?
	 *
	 * A non derivable(x) model is a model whose parameters depend on "x".
	 * It can be derived numerically.
	 */

	// Compute node functions

	class TransitionMatrixFromModelBrlenDerivative : public Value<MatrixDouble> {
	public:
		using Dependencies = TupleOfValues<const TransitionModel *, double>;

		TransitionMatrixFromModelBrlenDerivative (NodeRefVec && deps,
		                                          const TransitionMatrixDimension & dim)
		    : Value<MatrixDouble> (std::move (deps)) {
			setTargetDimension (this->accessValueMutable (), dim);
		}
		std::string debugInfo () const final {
			return Value<MatrixDouble>::debugInfo () + " " +
			       to_string (TransitionMatrixDimension (targetDimension (this->accessValueConst ())));
		}
		NodeRef derive (const Node & node) final {
			assert (isDerivable (node));
			auto dim = TransitionMatrixDimension (targetDimension (this->accessValueConst ()));
			auto & modelNode = this->dependency (0);
			auto & brlenNode = this->dependency (1);
			auto d2TransMat_dBrlen2 =
			    makeNode<TransitionMatrixFromModelBrlenSecondDerivative> ({modelNode, brlenNode}, dim);
			return makeNode<CWiseMulScalarMatrixDouble> (
			    {brlenNode->derive (node), std::move (d2TransMat_dBrlen2)}, dim);
		}
		bool isDerivable (const Node & node) const final { return derivableIfAllDepsAre (*this, node); }
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<TransitionMatrixFromModelBrlenDerivative> (
			    std::move (deps), targetDimension (this->accessValueConst ()));
		}

	private:
		void compute () final {
			callWithValues (*this, [](MatrixDouble & matrix, const TransitionModel * model,
			                          double brlen) { bppToEigen (model->getdPij_dt (brlen), matrix); });
		}
	};
	ValueRef<MatrixDouble>
	Builder<TransitionMatrixFromModelBrlenDerivative>::make (NodeRefVec && deps,
	                                                         const TransitionMatrixDimension & dim) {
		checkDependencies<TransitionMatrixFromModelBrlenDerivative> (deps);
		return std::make_shared<TransitionMatrixFromModelBrlenDerivative> (std::move (deps), dim);
	}

	class TransitionMatrixFromModelBrlenSecondDerivative : public Value<MatrixDouble> {
	public:
		using Dependencies = TupleOfValues<const TransitionModel *, double>;

		TransitionMatrixFromModelBrlenSecondDerivative (NodeRefVec && deps,
		                                                const TransitionMatrixDimension & dim)
		    : Value<MatrixDouble> (std::move (deps)) {
			setTargetDimension (this->accessValueMutable (), dim);
		}
		std::string debugInfo () const final {
			return Value<MatrixDouble>::debugInfo () + " " +
			       to_string (TransitionMatrixDimension (targetDimension (this->accessValueConst ())));
		}
		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<TransitionMatrixFromModelBrlenSecondDerivative> (
			    std::move (deps), targetDimension (this->accessValueConst ()));
		}

	private:
		void compute () final {
			callWithValues (*this,
			                [](MatrixDouble & matrix, const TransitionModel * model, double brlen) {
				                bppToEigen (model->getd2Pij_dt2 (brlen), matrix);
			                });
		}
	};
	ValueRef<MatrixDouble> Builder<TransitionMatrixFromModelBrlenSecondDerivative>::make (
	    NodeRefVec && deps, const TransitionMatrixDimension & dim) {
		checkDependencies<TransitionMatrixFromModelBrlenSecondDerivative> (deps);
		return std::make_shared<TransitionMatrixFromModelBrlenSecondDerivative> (std::move (deps), dim);
	}
} // namespace DF
} // namespace bpp
