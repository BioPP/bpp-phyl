//
// File: DataFlowNumericalDerivation.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-12-19 00:00:00
// Last modified: 2018-01-31
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
#include <Bpp/NewPhyl/DataFlowNumericalDerivation.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/IntegerRange.h>
#include <Bpp/NewPhyl/LinearAlgebraUtils.h>
#include <algorithm>

namespace bpp {
namespace DF {
	template <typename T>
	ValueRef<T> makeNumericalDerivationCombineShifted (NodeRefVec && deps,
	                                                   const std::vector<double> & coeffs,
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
			std::vector<double> mergedCoeffs;
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
						    static_cast<std::size_t> (std::distance (mergedDeps.begin () + 1, it));
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
} // namespace DF
} // namespace bpp
