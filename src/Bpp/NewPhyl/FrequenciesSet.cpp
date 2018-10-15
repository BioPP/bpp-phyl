//
// File: FrequenciesSet.cpp
// Authors: Laurent Guéguen
// Created: jeudi 11 octobre 2018, à 07h 10
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
#include <Bpp/NewPhyl/FrequenciesSet.h>
#include <Bpp/Phyl/Model/FrequenciesSet/FrequenciesSet.h>

using namespace std;

namespace bpp {
  namespace dataflow {
    // FrequenciesSet node

    std::shared_ptr<ConfiguredFrequenciesSet> ConfiguredFrequenciesSet::create (Context & c, NodeRefVec && deps,
                                                                                std::unique_ptr<FrequenciesSet> && freqset) {
      if (!freqset) {
        throw Exception ("ConfiguredFrequenciesSet(): nullptr freqset");
      }
      // Check dependencies
      const auto nbParameters = freqset->getParameters ().size ();
      checkDependenciesNotNull (typeid (Self), deps);
      checkDependencyVectorSize (typeid (Self), deps, nbParameters);
      checkDependencyRangeIsValue<double> (typeid (Self), deps, 0, nbParameters);
      return cachedAs<Self> (c, std::make_shared<Self> (std::move (deps), std::move (freqset)));
    }

    ConfiguredFrequenciesSet::ConfiguredFrequenciesSet (NodeRefVec && deps, std::unique_ptr<FrequenciesSet> && freqset)
      : ConfiguredParametrizable(std::move(deps), std::move(freqset)), freqset_ (dynamic_cast<const FrequenciesSet*>(getValue())) {}
    
    ConfiguredFrequenciesSet::~ConfiguredFrequenciesSet () = default;
    
    std::string ConfiguredFrequenciesSet::description () const { return "FreqSet(" + freqset_->getName () + ")"; }
    
    std::string ConfiguredFrequenciesSet::debugInfo () const {
      return "nbState=" + std::to_string (freqset_->getAlphabet ()->getSize ());
    }


    NodeRef ConfiguredFrequenciesSet::recreate (Context & c, NodeRefVec && deps) {
      auto m = Self::create (c, std::move (deps), std::unique_ptr<FrequenciesSet>{freqset_->clone ()});
      m->config = this->config; // Duplicate derivation config
      return m;
    }
    
  } // namespace dataflow
} // namespace bpp
