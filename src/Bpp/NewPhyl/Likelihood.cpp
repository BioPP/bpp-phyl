//
// File: Likelihood.cpp
// Authors:
// Created: 2017-06-06
// Last modified: 2017-06-06
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
#include <Bpp/NewPhyl/DataFlowInternalTemplates.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/ExtendedFloat.h>
#include <Bpp/NewPhyl/Likelihood.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/Seq/Sequence.h>
#include <cmath>

namespace bpp {
namespace Phyl {
	std::string LikelihoodDataDimension::toString () const {
		return "(sites=" + std::to_string (nbSites ()) + ",states=" + std::to_string (nbStates ()) +
		       ")";
	}

	namespace DF {
		ConditionalLikelihoodFromSequence::ConditionalLikelihoodFromSequence (
		    NodeRefVec && deps, LikelihoodDataDimension dim)
		    : Value<MatrixDouble> (std::move (deps), dim.rows, dim.cols) {
			checkDependencies (*this);
		}
		void ConditionalLikelihoodFromSequence::compute () {
			callWithValues (*this, [](MatrixDouble & condLikBySite, const Sequence * sequence) {
				// Check sizes
				if (sequence == nullptr)
					throw Exception (prettyTypeName<ConditionalLikelihoodFromSequence> () +
					                 ": null sequence");
				auto matDim = LikelihoodDataDimension (dimensions (condLikBySite));
				auto seqDim = LikelihoodDataDimension (static_cast<SizeType> (sequence->size ()),
				                                       sequence->getAlphabet ()->getSize ());
				if (matDim != seqDim)
					throw Exception (prettyTypeName<ConditionalLikelihoodFromSequence> () +
					                 ": size mismatch: sequence " + seqDim.toString () + " -> " +
					                 matDim.toString ());
				// Put 1s at the right places, 0s elsewhere
				condLikBySite.fill (0.);
				for (auto siteIndex : range (condLikBySite.cols ())) {
					auto siteValue =
					    static_cast<IndexType> (sequence->getValue (static_cast<std::size_t> (siteIndex)));
					condLikBySite (siteValue, siteIndex) = 1.;
				}
			});
		}
		std::string ConditionalLikelihoodFromSequence::debugInfo () const {
			// TODO add a common LikelihoodDataValue class that overrides this
			auto dim = LikelihoodDataDimension (dimensions (*this));
			return Value<MatrixDouble>::debugInfo () + " " + dim.toString ();
		}
		NodeRef ConditionalLikelihoodFromSequence::derive (const Node &) {
			// Sequence is a constant with respect to all parameters.
			return Builder<Constant<MatrixDouble>>::makeZero (dimensions (*this));
		}
		std::shared_ptr<ConditionalLikelihoodFromSequence>
		ConditionalLikelihoodFromSequence::create (NodeRefVec && deps, LikelihoodDataDimension dim) {
			return std::make_shared<ConditionalLikelihoodFromSequence> (std::move (deps), dim);
		}
		std::shared_ptr<ConditionalLikelihoodFromSequence>
		ConditionalLikelihoodFromSequence::create (ValueRef<const Sequence *> sequence,
		                                           LikelihoodDataDimension dim) {
			return create (NodeRefVec{std::move (sequence)}, dim);
		}

		TotalLogLikelihood::TotalLogLikelihood (NodeRefVec && deps) : Value<double> (std::move (deps)) {
			checkDependencies (*this);
		}
		void TotalLogLikelihood::compute () {
			callWithValues (*this, [](double & logLik, const VectorDouble & likelihood) {
				auto lik = likelihood
				               .unaryExpr ([](double d) {
					               ExtendedFloat ef{d};
					               ef.normalize_small ();
					               return ef;
				               })
				               .redux ([](const ExtendedFloat & lhs, const ExtendedFloat & rhs) {
					               auto r = denorm_mul (lhs, rhs);
					               r.normalize_small ();
					               return r;
				               });
				logLik = log (lik);
			});
		}
		NodeRef TotalLogLikelihood::derive (const Node & node) {
			// TODO improve this case
			auto likelihoodVector = convertRef<Value<VectorDouble>> (this->dependencies ()[0]);
			return ScalarProdDouble::create (
			    NodeRefVec{likelihoodVector->derive (node),
			               CWiseInverseVectorDouble::create (NodeRefVec{likelihoodVector},
			                                                 dimensions (*likelihoodVector))

			    });
		}
		std::shared_ptr<TotalLogLikelihood> TotalLogLikelihood::create (NodeRefVec && deps) {
			return std::make_shared<TotalLogLikelihood> (std::move (deps));
		}
		std::shared_ptr<TotalLogLikelihood>
		TotalLogLikelihood::create (ValueRef<VectorDouble> likelihood) {
			return create (NodeRefVec{std::move (likelihood)});
		}
	} // namespace DF
} // namespace Phyl
} // namespace bpp
