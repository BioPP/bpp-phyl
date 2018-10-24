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
      : Value<const FrequenciesSet*> (std::move (deps), freqset.get ()), freqset_(std::move(freqset)) {}

    ConfiguredFrequenciesSet::~ConfiguredFrequenciesSet () = default;
    
    std::string ConfiguredFrequenciesSet::description () const { return "FreqSet(" + freqset_->getName () + ")"; }
    
    std::string ConfiguredFrequenciesSet::debugInfo () const {
      return "nbState=" + std::to_string (freqset_->getAlphabet ()->getSize ());
    }

    const std::string & ConfiguredFrequenciesSet::getParameterName (std::size_t index) {
      return freqset_->getParameters ()[index].getName ();
    }
    
    std::size_t ConfiguredFrequenciesSet::getParameterIndex (const std::string & name) {
      return static_cast<std::size_t> (freqset_->getParameters ().whichParameterHasName (name));
    }

    // FrequenciesSet node additional arguments = (type of bpp::FrequenciesSet).
    // Everything else is determined by the node dependencies.

    bool ConfiguredFrequenciesSet::compareAdditionalArguments (const Node & other) const {
      const auto * derived = dynamic_cast<const Self *> (&other);
      if (derived == nullptr) {
        return false;
      } else {
        const auto & thisFS = *freqset_;
        const auto & otherFS = *derived->freqset_;
        return typeid (thisFS) == typeid (otherFS);
      }
    }
    
    std::size_t ConfiguredFrequenciesSet::hashAdditionalArguments () const {
      const auto & bppFS = *freqset_;
      return typeid (bppFS).hash_code ();
    }


    NodeRef ConfiguredFrequenciesSet::recreate (Context & c, NodeRefVec && deps) {
      auto m = Self::create (c, std::move (deps), std::unique_ptr<FrequenciesSet>{freqset_->clone ()});
      m->config = this->config; // Duplicate derivation config
      return m;
    }

    void ConfiguredFrequenciesSet::compute () {
      // Update each internal freqseq bpp::Parameter with the dependency
      auto & parameters = freqset_->getParameters ();
      const auto nbParameters = this->nbDependencies ();
      for (std::size_t i = 0; i < nbParameters; ++i) {
        auto & v = accessValueConstCast<double> (*this->dependency (i));
        auto & p = parameters[i];
        if (p.getValue () != v) {
          // TODO improve bpp::Parametrizable interface to change values by index.
          freqset_->setParameterValue (freqset_->getParameterNameWithoutNamespace (p.getName ()), v);
        }
      }
    }


    
  } // namespace dataflow
} // namespace bpp
