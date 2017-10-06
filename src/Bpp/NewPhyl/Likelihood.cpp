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

#include <Bpp/NewPhyl/DataFlowTemplateUtils.h>
#include <Bpp/NewPhyl/ExtendedFloat.h>
#include <Bpp/NewPhyl/Likelihood.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/Seq/Sequence.h>
#include <cassert>
#include <cmath>

namespace bpp {
namespace Phyl {
	ComputeConditionalLikelihoodFromDataNode::ComputeConditionalLikelihoodFromDataNode (
	    DF::NodeRefVec && deps, SizeType nbSites, SizeType nbStates)
	    : DF::Value<LikelihoodVectorBySite> (std::move (deps), nbSites, nbStates) {
		DF::checkDependencies (*this);
	}
	void ComputeConditionalLikelihoodFromDataNode::compute () {
		DF::callWithValues (
		    *this, [](LikelihoodVectorBySite & condLikBySite, const Sequence * sequence) {
			    assert (sequence != nullptr);
			    assert (condLikBySite.size () == sequence->size ());
			    for (auto siteIndex : index_range (condLikBySite)) {
				    LikelihoodVectorBySite::reference lik = condLikBySite[siteIndex];
				    assert (lik.size () == sequence->getAlphabet ()->getSize ());
				    lik.fill (0.);
				    auto siteValue =
				        static_cast<IndexType> (sequence->getValue (static_cast<std::size_t> (siteIndex)));
				    lik[siteValue] = 1.;
			    }
			  });
	}
	std::string ComputeConditionalLikelihoodFromDataNode::description () const {
		return "CondLikFromData";
	}

	ComputeConditionalLikelihoodFromChildrensNode::ComputeConditionalLikelihoodFromChildrensNode (
	    DF::NodeRefVec && deps, SizeType nbSites, SizeType nbStates)
	    : DF::Value<LikelihoodVectorBySite> (std::move (deps), nbSites, nbStates) {
		DF::checkDependencies (*this);
	}
	void ComputeConditionalLikelihoodFromChildrensNode::compute () {
		DF::callWithValues (
		    *this, [](LikelihoodVectorBySite & condLikBySite) { condLikBySite.asMatrix ().fill (1.); },
		    [](LikelihoodVectorBySite & condLikBySite, const LikelihoodVectorBySite & fwdLikBySite) {
			    condLikBySite.asMatrix () =
			        condLikBySite.asMatrix ().cwiseProduct (fwdLikBySite.asMatrix ());
			  });
	}
	std::string ComputeConditionalLikelihoodFromChildrensNode::description () const {
		return "CondLikFromChildrens";
	}

	ComputeForwardLikelihoodNode::ComputeForwardLikelihoodNode (DF::NodeRefVec && deps,
	                                                            SizeType nbSites, SizeType nbStates)
	    : DF::Value<LikelihoodVectorBySite> (std::move (deps), nbSites, nbStates) {
		DF::checkDependencies (*this);
	}
	void ComputeForwardLikelihoodNode::compute () {
		DF::callWithValues (*this, [](LikelihoodVectorBySite & fwdLikBySite,
		                              const LikelihoodVectorBySite & condLikBySite,
		                              const TransitionMatrix & transitionMatrix) {
			fwdLikBySite.asMatrix ().noalias () = transitionMatrix * condLikBySite.asMatrix ();
		});
	}
	std::string ComputeForwardLikelihoodNode::description () const { return "FwdLik"; }

	ComputeLogLikelihoodNode::ComputeLogLikelihoodNode (DF::NodeRefVec && deps)
	    : DF::Value<double> (std::move (deps)) {
		DF::checkDependencies (*this);
	}
	void ComputeLogLikelihoodNode::compute () {
		DF::callWithValues (*this, [](double & logLik, const LikelihoodVectorBySite & condLikBySite,
		                              const FrequencyVector & equilibriumFreqs) {
			auto lik = (equilibriumFreqs.transpose () * condLikBySite.asMatrix ())
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
	std::string ComputeLogLikelihoodNode::description () const { return "LogLikFromCondLik"; }
}
}
