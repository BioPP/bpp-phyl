//
// File: FrequencySet.cpp
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
#include <Bpp/Phyl/NewLikelihood/DataFlow/FrequencySet.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/Parametrizable.h>

using namespace std;

namespace bpp {
    // FrequencySet node

    ConfiguredFrequencySet::ConfiguredFrequencySet (const Context& context, NodeRefVec && deps, std::unique_ptr<FrequencySet> && freqset)
      : Value<const FrequencySet*> (std::move (deps), freqset.get ()), AbstractParametrizable(freqset->getNamespace()), context_(context), freqset_(std::move(freqset))
    {
      for (const auto& dep:dependencies())
      {
        const auto& param=std::dynamic_pointer_cast<ConfiguredParameter>(dep);
        shareParameter_(param);
      }
    }

    ConfiguredFrequencySet::~ConfiguredFrequencySet () = default;
    
    std::string ConfiguredFrequencySet::description () const { return "FreqSet(" + freqset_->getName () + ")"; }
    
    std::string ConfiguredFrequencySet::debugInfo () const {
      return "nbState=" + std::to_string (freqset_->getAlphabet ()->getSize ());
    }

    // FrequencySet node additional arguments = (type of bpp::FrequencySet).
    // Everything else is determined by the node dependencies.

    bool ConfiguredFrequencySet::compareAdditionalArguments (const Node_DF & other) const {
      const auto * derived = dynamic_cast<const Self *> (&other);
      if (derived == nullptr) {
        return false;
      } else {
        const auto & thisFS = *freqset_;
        const auto & otherFS = *derived->freqset_;
        return typeid (thisFS) == typeid (otherFS);
      }
    }
    
    std::size_t ConfiguredFrequencySet::hashAdditionalArguments () const {
      const auto & bppFS = *freqset_;
      return typeid (bppFS).hash_code ();
    }


    NodeRef ConfiguredFrequencySet::recreate (Context & c, NodeRefVec && deps) {
      auto m = ConfiguredParametrizable::createConfigured<Target, Self> (c, std::move (deps), std::unique_ptr<Target>{freqset_->clone ()});
      m->config = this->config; // Duplicate derivation config
      return m;
    }

    // FrequenciesFromFrequencySet

    FrequenciesFromFrequencySet::FrequenciesFromFrequencySet (
      NodeRefVec && deps, const Dimension<Eigen::RowVectorXd> & dim)
      : Value<Eigen::RowVectorXd> (std::move (deps)), targetDimension_ (dim) {}

    std::string FrequenciesFromFrequencySet::debugInfo () const {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    // FrequenciesFromFrequencySet additional arguments = ().
    bool FrequenciesFromFrequencySet::compareAdditionalArguments (const Node_DF & other) const {
      return dynamic_cast<const Self *> (&other) != nullptr;
    }

    NodeRef FrequenciesFromFrequencySet::derive (Context & c, const Node_DF & node) {
      // d(equFreqs)/dn = sum_i d(equFreqs)/dx_i * dx_i/dn (x_i = freqset parameters)
      auto freqSetDep = this->dependency (0);
      auto & freqset = static_cast<ConfiguredFrequencySet &> (*freqSetDep);
      auto buildFWithNewFreqSet = [this, &c](NodeRef && newFreqSet) {
        return ConfiguredParametrizable::createVector<ConfiguredFrequencySet, Self> (c, {std::move (newFreqSet)}, targetDimension_);
      };
      
      NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<ConfiguredFrequencySet, T > (
        c, freqset, node, targetDimension_, buildFWithNewFreqSet);
      return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
    }

    NodeRef FrequenciesFromFrequencySet::recreate (Context & c, NodeRefVec && deps) {
      return ConfiguredParametrizable::createVector<ConfiguredFrequencySet, Self> (c, std::move (deps), targetDimension_);
    }

    void FrequenciesFromFrequencySet::compute () {
      const auto * freqset = accessValueConstCast<const FrequencySet *> (*this->dependency (0));
      const auto & freqsFromFS = freqset->getFrequencies ();
      auto & r = this->accessValueMutable ();
      r = Eigen::Map<const T> (freqsFromFS.data(), static_cast<Eigen::Index> (freqsFromFS.size ()));
    }

  
} // namespace bpp
