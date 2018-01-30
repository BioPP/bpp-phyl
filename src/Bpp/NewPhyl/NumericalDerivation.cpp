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

#include <Bpp/NewPhyl/DataFlowInternal.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/IntegerRange.h>
#include <Bpp/NewPhyl/LinearAlgebraUtils.h>
#include <Bpp/NewPhyl/NumericalDerivation.h>
#include <algorithm>

namespace bpp {
namespace DF {
	/* Data flow nodes in this file are implemented using the 'new' pattern.
	 * The class is forward declared as template (header file).
	 * Templated by the type, to have one code for double / VectorDouble / MatrixDouble.
	 * It is implemented as an actual template internally (this file).
	 * Then the factory function is implemented as a template as well.
	 * Finally, explicit specialisations of Builder<T>::make are instanciated here.
	 */

	// NumericalDerivationShiftDelta.
	template <typename T> class NumericalDerivationShiftDelta : public Value<T> {
	public:
		using Dependencies = TupleOfValues<double, T>;

		NumericalDerivationShiftDelta (NodeRefVec && deps, int n, const Dimension<T> & targetDim)
		    : Value<T> (std::move (deps)), n_ (n) {
			this->setTargetDimension (targetDim);
		}

		std::string description () const final { return std::to_string (n_) + " * delta * x"; }

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
				result = linearAlgebraMakeValueWith (this->getTargetDimension (), this->n_ * delta) + arg;
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
			// n * delta * (m * delta * x) -> (n + m) * delta * x
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

	ValueRef<VectorDouble> Builder<NumericalDerivationShiftDelta<VectorDouble>>::make (
	    NodeRefVec && deps, int n, const Dimension<VectorDouble> & targetDim) {
		return makeNumericalDerivationShiftDelta<VectorDouble> (std::move (deps), n, targetDim);
	}

	ValueRef<MatrixDouble> Builder<NumericalDerivationShiftDelta<MatrixDouble>>::make (
	    NodeRefVec && deps, int n, const Dimension<MatrixDouble> & targetDim) {
		return makeNumericalDerivationShiftDelta<MatrixDouble> (std::move (deps), n, targetDim);
	}

	/* NumericalDerivationCombineShifted
	 *
	 * Combine shifted has uncommon dependencies: double followed by ArrayOfValues<T>.
	 * This is not yet representable in the dependency pattern type tag system.
	 * Checks and compute are done by manually combining primitives.
	 */
	template <typename T> class NumericalDerivationCombineShifted : public Value<T> {
	public:
		NumericalDerivationCombineShifted (NodeRefVec && deps, const Vector<double> & coeffs,
		                                   const Dimension<T> & targetDim)
		    : Value<T> (std::move (deps)), coeffs_ (coeffs) {
			this->setTargetDimension (targetDim);
		}

		std::string description () const final {
			auto r = std::string ("lambda * sum (deps[i] * {");
			if (!coeffs_.empty ()) {
				for (auto i : range (coeffs_.size () - 1))
					r += std::to_string (coeffs_[i]) + " ";
				r += std::to_string (coeffs_.back ());
			}
			r += "}[i])";
			return r;
		}

		// Never derive the lambda side (not part of computation !)
		NodeRef derive (const Node & node) final {
			NodeRefVec deps (this->nbDependencies ());
			deps[0] = this->dependency (0); // lambda
			for (auto i : range (SizeType (1), this->nbDependencies ()))
				deps[i] = this->dependency (i)->derive (node);
			return makeNode<NumericalDerivationCombineShifted<T>> (std::move (deps), coeffs_,
			                                                       this->getTargetDimension ());
		}
		bool isDerivable (const Node & node) const final {
			const auto & deps = this->dependencies ();
			auto & lambda = deps[0];
			return std::all_of (deps.begin () + 1, deps.end (),
			                    [&node](const NodeRef & ref) { return ref->isDerivable (node); }) &&
			       !lambda->isTransitivelyDependentOn (node);
		}

		NodeRef rebuild (NodeRefVec && deps) const final {
			return makeNode<NumericalDerivationCombineShifted<T>> (std::move (deps), coeffs_,
			                                                       this->getTargetDimension ());
		}

		const Vector<double> & getCoeffs () const noexcept { return coeffs_; }

	private:
		Vector<double> coeffs_;

		void compute () final {
			const auto & deps = this->dependencies ();
			double lambda = accessValidValueConstCast<double> (deps[0]);
			T & value = this->accessValueMutable ();

			value = linearAlgebraZeroValue (this->getTargetDimension ());
			for (auto i : range (coeffs_.size ())) {
				value += linearAlgebraMakeValueWith (this->getTargetDimension (), coeffs_[i] * lambda) *
				         accessValidValueConstCast<T> (deps[i + 1]);
			}
		}
	};

	template <typename T>
	ValueRef<T> makeNumericalDerivationCombineShifted (NodeRefVec && deps,
	                                                   const Vector<double> & coeffs,
	                                                   const Dimension<T> & targetDim) {
		// CombineShifted dependencies are manually checked using primitives.
		auto && typeId = typeid (NumericalDerivationCombineShifted<T>);
		checkDependencyVectorSize (typeId, deps, 1 + coeffs.size ());
		checkDependenciesNotNull (typeId, deps);
		checkNthDependencyIsValue<double> (typeId, deps, 0); // lambda
		checkDependencyPatternImpl (typeId, deps, 1, ArrayOfValues<T>{coeffs.size ()});

		auto & lambda = deps[0];

		// Merge 2 layers of NumericalDerivationCombineShifted if they have the same lambda.
		auto canMergeThisDep = [&lambda](const NodeRef & dep) {
			return isNodeType<NumericalDerivationCombineShifted<T>> (*dep) &&
			       dep->dependency (0) == lambda;
		};
		if (std::all_of (deps.begin () + 1, deps.end (), canMergeThisDep)) {
			/* Shortening:
			 * L = lambda.
			 * V = Current node value.
			 * DC[i] = coeffs of deps[i]. In these, deps[i] ignores lambda (deps[0]).
			 * DD[i] = deps of deps[i].
			 *
			 * V = L * sum_i (deps[i] * coeffs[i])
			 * After the test, we decompose deps:
			 * V = L * sum_i (L * sum_j (DD[i][j] * DC[i][j]) * coeffs[i]
			 * V = L^2 * sum_ij (DD[i][j] * DC[i][j] * coeffs[i])
			 * Finally, we merge identical dependencies by summing their coefficients.
			 */

			// Create L^2
			NodeRefVec mergedDeps;
			mergedDeps.emplace_back (makeNode<MulDouble> ({lambda, lambda}));

			// Add unique dependencies with merged coefficients
			Vector<double> mergedCoeffs;
			for (auto i : range (coeffs.size ())) {
				auto & dep = nodeCast<NumericalDerivationCombineShifted<T>> (*deps[i + 1]);
				auto & c = coeffs[i];
				for (auto j : range (dep.getCoeffs ().size ())) {
					auto & subDep = dep.dependency (j + 1);
					auto & subCoeff = dep.getCoeffs ()[j];

					auto it = std::find (mergedDeps.begin () + 1, mergedDeps.end (), subDep);
					if (it != mergedDeps.end ()) {
						// If found, sum to merged coefficient
						auto mergedDepIndexIgnoringLambda =
						    static_cast<SizeType> (std::distance (mergedDeps.begin () + 1, it));
						mergedCoeffs[mergedDepIndexIgnoringLambda] += c * subCoeff;
					} else {
						// Not found, add dep and new coefficient
						mergedDeps.emplace_back (subDep);
						mergedCoeffs.emplace_back (c * subCoeff);
					}
				}
			}

			return makeNode<NumericalDerivationCombineShifted<T>> (std::move (mergedDeps), mergedCoeffs,
			                                                       targetDim);
		} else {
			// Non merging case
			return std::make_shared<NumericalDerivationCombineShifted<T>> (std::move (deps), coeffs,
			                                                               targetDim);
		}
	}

	ValueRef<double> Builder<NumericalDerivationCombineShifted<double>>::make (
	    NodeRefVec && deps, const Vector<double> & coeffs, const Dimension<double> & targetDim) {
		return makeNumericalDerivationCombineShifted<double> (std::move (deps), coeffs, targetDim);
	}
	ValueRef<VectorDouble> Builder<NumericalDerivationCombineShifted<VectorDouble>>::make (
	    NodeRefVec && deps, const Vector<double> & coeffs,
	    const Dimension<VectorDouble> & targetDim) {
		return makeNumericalDerivationCombineShifted<VectorDouble> (std::move (deps), coeffs,
		                                                            targetDim);
	}
	ValueRef<MatrixDouble> Builder<NumericalDerivationCombineShifted<MatrixDouble>>::make (
	    NodeRefVec && deps, const Vector<double> & coeffs,
	    const Dimension<MatrixDouble> & targetDim) {
		return makeNumericalDerivationCombineShifted<MatrixDouble> (std::move (deps), coeffs,
		                                                            targetDim);
	}

} // namespace DF
} // namespace bpp
