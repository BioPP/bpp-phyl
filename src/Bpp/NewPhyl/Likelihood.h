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

#include <Bpp/NewPhyl/Dataflow.h>
#include <Eigen/Dense>
#include <cassert>
#include <vector>

namespace bpp {

using LikelihoodVector = Eigen::VectorXd;
using LikelihoodVectorBySite = std::vector<LikelihoodVector>;

// TODO have special types for Nucleotide (N=4) ?
using NucleotideLikelihoodVector = Eigen::Vector4d;

namespace DF {

	class ConditionalLikelihood : public Value<LikelihoodVectorBySite>::Impl {
	public:
		using ArgumentType = Value<LikelihoodVectorBySite>;

		ConditionalLikelihood (std::size_t nbCharacters, std::size_t nbSites)
		    : Value<LikelihoodVectorBySite>::Impl (nbSites, LikelihoodVector (nbCharacters)) {}

		void addDependency (ArgumentType forwardLikelihood) {
			this->appendDependency (Node (forwardLikelihood));
		}

	private:
		void compute () override final {
			auto & deps = this->dependencyNodes_;
			assert (!deps.empty ());
      // Store vect of reference wrappers to deps values
      // Process liks site by site
			this->value_ = static_cast<ArgumentType::Ref> (deps.front ().getImpl ()).getValue ();
			for (auto it = deps.begin () + 1; it != deps.end (); ++it)
				;
		}
	};
}
}

#endif // BPP_NEWPHYL_LIKELIHOOD_H
