//
// File: NumericalDerivation.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-12-19
// Last modified: 2017-12-19
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

#include <Bpp/NewPhyl/DataFlowInternalTemplates.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/LinearAlgebra.h>
#include <Bpp/NewPhyl/NumericalDerivation.h>

namespace bpp {
namespace DF {
  /* NumericalDerivationShiftDelta.
   *
   * To avoid code duplication, the class is implemented as a template (private to this file).
   * The template is instantianted with explicit types due to the factory functions.
   */
	template <typename T> class NumericalDerivationShiftDelta : public Value<T> {
	public:
		using Dependencies = FunctionOfValues<double, T>;

		NumericalDerivationShiftDelta (NodeRefVec && deps, int n, const Dimension<T> & targetDim)
		    : Value<T> (std::move (deps)), n_ (n) {
			this->setTargetDimension (targetDim);
		}

		std::string description () const final {
			return std::to_string (n_) + " * delta * " + prettyTypeName<T> ();
		}

		// Never derive the delta side (not part of computation !)
		NodeRef derive (const Node & node) final {
			auto & delta = this->dependency (0);
			auto & x = this->dependency (1);
			return makeNode<NumericalDerivationShiftDelta<T>> ({delta, x->derive (node)}, n_,
			                                                   this->getTargetDimension ());
		}
		bool isDerivable (const Node & node) const final {
			auto & delta = this->dependency (0);
			auto & x = this->dependency (1);
			return x->isDerivable (node) && !delta->isTransitivelyDependentOn (node);
		}

		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<NumericalDerivationShiftDelta<T>> (std::move (deps), n_,
			                                                   this->getTargetDimension ());
		}

		int getShiftFactor () const noexcept { return n_; }

	private:
		int n_;

		void compute () final {
			callWithValues (*this, [this](T & result, const double & delta, const T & arg) {
				result = linearAlgebraValueFilledWith (this->getTargetDimension (), this->n_ * delta) + arg;
			});
		}
	};

	template <typename T>
	ValueRef<T> makeNumericalDerivationShiftDelta (NodeRefVec && deps, int n,
	                                               const Dimension<T> & targetDim) {
		checkDependencies<NumericalDerivationShiftDelta<T>> (deps);
		auto & delta = deps[0];
		auto & x = deps[1];
		auto * xAsShiftDelta = dynamic_cast<const NumericalDerivationShiftDelta<T> *> (x.get ());
		if (xAsShiftDelta != nullptr && xAsShiftDelta->dependency (0) == delta) {
			// Merge chains of ShiftDelta nodes recursively
			// n*delta* (m*delta*x) -> (n+m)*delta*x
			return makeNode<NumericalDerivationShiftDelta<T>> (
			    NodeRefVec{x->dependencies ()}, n + xAsShiftDelta->getShiftFactor (), targetDim);
		} else {
			// Final node choice (remove ShiftDelta node if n is 0)
			if (n == 0) {
				return convertRef<Value<T>> (x);
			} else {
				return std::make_shared<NumericalDerivationShiftDelta<T>> (std::move (deps), n, targetDim);
			}
		}
	}

	ValueRef<double>
	Builder<NumericalDerivationShiftDelta<double>>::make (NodeRefVec && deps, int n,
	                                                      const Dimension<double> & targetDim) {
		return makeNumericalDerivationShiftDelta<double> (std::move (deps), n, targetDim);
	}
	ValueRef<double> Builder<NumericalDerivationShiftDelta<double>>::make (NodeRefVec && deps,
	                                                                       int n) {
		return make (std::move (deps), n, Dimension<double>{});
	}

	ValueRef<VectorDouble> Builder<NumericalDerivationShiftDelta<VectorDouble>>::make (
	    NodeRefVec && deps, int n, const Dimension<VectorDouble> & targetDim) {
		return makeNumericalDerivationShiftDelta<VectorDouble> (std::move (deps), n, targetDim);
	}

	ValueRef<MatrixDouble> Builder<NumericalDerivationShiftDelta<MatrixDouble>>::make (
	    NodeRefVec && deps, int n, const Dimension<MatrixDouble> & targetDim) {
		return makeNumericalDerivationShiftDelta<MatrixDouble> (std::move (deps), n, targetDim);
	}

  // NumericalDerivationCombineShifted

} // namespace DF
} // namespace bpp
