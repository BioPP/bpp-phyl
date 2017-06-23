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

#include <Bpp/NewPhyl/Likelihood.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/Seq/Sequence.h>
#include <cassert>
#include <cmath>

namespace bpp {
namespace Phyl {
	void ComputeConditionalLikelihoodFromDataOp::compute (LikelihoodVectorBySite & condLikBySite,
	                                                      const Sequence * sequence) {
		assert (sequence != nullptr);
		assert (condLikBySite.size () == sequence->size ());
		for (auto siteIndex : index_range (condLikBySite)) {
			auto & lik = condLikBySite[siteIndex];
			assert (lik.size () == sequence->getAlphabet ()->getSize ());
			lik.fill (0.);
			auto siteValue = static_cast<int> (sequence->getValue (siteIndex));
			lik[siteValue] = 1.;
		}
	}

	void ComputeConditionalLikelihoodFromChildrensOp::reset (LikelihoodVectorBySite & condLikBySite) {
		for (auto & lik : condLikBySite)
			lik.fill (1.);
	}
	void ComputeConditionalLikelihoodFromChildrensOp::reduce (
	    LikelihoodVectorBySite & condLikBySite, const LikelihoodVectorBySite & fwdLikBySite) {
		assert (condLikBySite.size () == fwdLikBySite.size ());
		for (auto siteIndex : index_range (condLikBySite))
			condLikBySite[siteIndex] = condLikBySite[siteIndex].cwiseProduct (fwdLikBySite[siteIndex]);
	}

	void ComputeForwardLikelihoodOp::compute (LikelihoodVectorBySite & fwdLikBySite,
	                                          const LikelihoodVectorBySite & condLikBySite,
	                                          const TransitionMatrix & transitionMatrix) {
		assert (fwdLikBySite.size () == condLikBySite.size ());
		for (auto siteIndex : index_range (fwdLikBySite))
			fwdLikBySite[siteIndex].noalias () = transitionMatrix * condLikBySite[siteIndex];
	}

	void ComputeLogLikelihoodOp::compute (double & logLik,
	                                      const LikelihoodVectorBySite & condLikBySite,
	                                      const FrequencyVector & equilibriumFreqs) {
		logLik = 0.;
		for (auto & condLik : condLikBySite)
			logLik += std::log (condLik.dot (equilibriumFreqs));
	}
}
}
