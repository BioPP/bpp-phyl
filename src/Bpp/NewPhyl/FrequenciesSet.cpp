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
#include <Bpp/NewPhyl/Parametrizable.h>
#include <Bpp/Phyl/Model/FrequenciesSet/FrequenciesSet.h>

using namespace std;

namespace bpp {
  namespace dataflow {
    // FrequenciesSet node

    ConfiguredFrequenciesSet::ConfiguredFrequenciesSet (const Context& context, NodeRefVec && deps, std::unique_ptr<FrequenciesSet> && freqset)
      : Value<const FrequenciesSet*> (std::move (deps), freqset.get ()), context_(context), freqset_(std::move(freqset)) {}

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
      auto m = ConfiguredParametrizable::createConfigured<Target, Self> (c, std::move (deps), std::unique_ptr<Target>{freqset_->clone ()});
      m->config = this->config; // Duplicate derivation config
      return m;
    }

    void ConfiguredFrequenciesSet::compute () {
      // Update each internal freqseq bpp::Parameter with the dependency
      auto & parameters = freqset_->getParameters ();
      const auto nbParameters = this->nbDependencies ();
      for (std::size_t i = 0; i < nbParameters; ++i) {
        auto v = accessValueConstCast<Parameter*> (*this->dependency (i))->getValue();
        auto & p = parameters[i];
        if (p.getValue () != v) {
          // TODO improve bpp::Parametrizable interface to change values by index.
          freqset_->setParameterValue (freqset_->getParameterNameWithoutNamespace (p.getName ()), v);
        }
      }
    }

    // FrequenciesFromFrequenciesSet

    FrequenciesFromFrequenciesSet::FrequenciesFromFrequenciesSet (
      NodeRefVec && deps, const Dimension<Eigen::RowVectorXd> & dim)
      : Value<Eigen::RowVectorXd> (std::move (deps)), targetDimension_ (dim) {}

    std::string FrequenciesFromFrequenciesSet::debugInfo () const {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    // FrequenciesFromFrequenciesSet additional arguments = ().
    bool FrequenciesFromFrequenciesSet::compareAdditionalArguments (const Node & other) const {
      return dynamic_cast<const Self *> (&other) != nullptr;
    }

    NodeRef FrequenciesFromFrequenciesSet::derive (Context & c, const Node & node) {
      // d(equFreqs)/dn = sum_i d(equFreqs)/dx_i * dx_i/dn (x_i = freqset parameters)
      auto freqSetDep = this->dependency (0);
      auto & freqset = static_cast<ConfiguredFrequenciesSet &> (*freqSetDep);
      auto buildFWithNewFreqSet = [this, &c](NodeRef && newFreqSet) {
        return ConfiguredParametrizable::createVector<ConfiguredFrequenciesSet, Self> (c, {std::move (newFreqSet)}, targetDimension_);
      };
      
      NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<ConfiguredFrequenciesSet, T > (
        c, freqset, node, targetDimension_, buildFWithNewFreqSet);
      return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
    }

    NodeRef FrequenciesFromFrequenciesSet::recreate (Context & c, NodeRefVec && deps) {
      return ConfiguredParametrizable::createVector<ConfiguredFrequenciesSet, Self> (c, std::move (deps), targetDimension_);
    }

    void FrequenciesFromFrequenciesSet::compute () {
      const auto * freqset = accessValueConstCast<const FrequenciesSet *> (*this->dependency (0));
      const auto & freqsFromFS = freqset->getFrequencies ();
      auto & r = this->accessValueMutable ();
      r = Eigen::Map<const T> (freqsFromFS.data(), static_cast<Eigen::Index> (freqsFromFS.size ()));
    }


    
  } // namespace dataflow
} // namespace bpp
