//
// File: Likelihood.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-03
// Last modified: 2017-05-03
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

#pragma once
#ifndef BPP_NEWPHYL_LIKELIHOOD_H
#define BPP_NEWPHYL_LIKELIHOOD_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/NewPhyl/Signed.h>
#include <string>

namespace bpp {
class Sequence;

namespace Phyl {
	namespace DF {
		using namespace bpp::DF;
	}

	/* Likelihood probabilities (final or intermediate) are stored in a matrix.
	 * Frequencies for site k are stored in column k.
	 * TODO accessors ?
	 */
	using LikelihoodData = DF::MatrixDouble;

	// defines a MatrixDimension compatible struct.
	struct LikelihoodDataDimension : public DF::MatrixDimension {
		constexpr LikelihoodDataDimension (SizeType nbSitesArg, SizeType nbStatesArg) noexcept
		    : DF::MatrixDimension (nbStatesArg, nbSitesArg) {}
		constexpr LikelihoodDataDimension (const DF::MatrixDimension & matDim) noexcept
		    : DF::MatrixDimension (matDim) {}

		SizeType nbStates () const { return rows; }
		SizeType nbSites () const { return cols; }
		std::string toString () const;
	};

	/* TODO: <list>
	 * - wrapper classes.
	 * - take LikelihoodDataDimension as dims.
	 * - override debugInfo
	 * - externalize Likelihood as a numeric op
	 */

	namespace DF {
		struct ConditionalLikelihoodFromSequence : public Value<MatrixDouble> {
			// (sequence) -> MatrixDouble
			using Dependencies = FunctionOfValues<const Sequence *>;
			ConditionalLikelihoodFromSequence (NodeRefVec && deps, LikelihoodDataDimension dim);
			void compute () override final;
			std::string debugInfo () const override final;
			NodeRef derive (const Node &) override final;
		};

		// vec<fwdLik> -> condLik
		using ConditionalLikelihoodFromChildrens = CWiseMulMatrixDouble;

		// (transitionMatrix, condLik) -> fwdLik
		using ForwardLikelihoodFromChild = MulMatrixDouble;

		// (condLik, equFreqs) -> likBySiteVector
		using Likelihood = MulTransposedMatrixVectorDouble;

		struct TotalLogLikelihood : public Value<double> {
			// likelihood by site -> total log likelihood
			using Dependencies = FunctionOfValues<VectorDouble>;
			TotalLogLikelihood (NodeRefVec && deps);
			void compute () override final;
			NodeRef derive (const Node &) override final;
		};
	} // namespace DF
} // namespace Phyl
} // namespace bpp

#endif // BPP_NEWPHYL_LIKELIHOOD_H
