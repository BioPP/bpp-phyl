//
// File: Simplex.cpp
// Authors: Laurent Guéguen
// Created: vendredi 3 avril 2020, à 08h 32
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
#include <Bpp/Phyl/NewLikelihood/DataFlow/Simplex_DF.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/Parametrizable.h>

using namespace std;

namespace bpp {
  // Simplex node

  ConfiguredSimplex::ConfiguredSimplex (const Context& context, NodeRefVec && deps, std::unique_ptr<Simplex> && simplex)
    : Value<const Simplex*> (std::move (deps), simplex.get ()), AbstractParametrizable(simplex->getNamespace()), context_(context), simplex_(std::move(simplex))
  {
    for (const auto& dep:dependencies())
    {
      const auto& param=std::dynamic_pointer_cast<ConfiguredParameter>(dep);
      shareParameter_(param);
    }
  }

  ConfiguredSimplex::~ConfiguredSimplex () = default;
    
  std::string ConfiguredSimplex::description () const { return "Simplex"; }
    
  std::string ConfiguredSimplex::debugInfo () const {
    return "nbState=" + std::to_string (simplex_->dimension ());
  }

  // Simplex node additional arguments = (type of bpp::Simplex).
  // Everything else is determined by the node dependencies.

  bool ConfiguredSimplex::compareAdditionalArguments (const Node_DF & other) const {
    const auto * derived = dynamic_cast<const Self *> (&other);
    if (derived == nullptr) {
      return false;
    } else {
      const auto & thisFS = *simplex_;
      const auto & otherFS = *derived->simplex_;
      return typeid (thisFS) == typeid (otherFS);
    }
  }
    
  std::size_t ConfiguredSimplex::hashAdditionalArguments () const {
    const auto & bppFS = *simplex_;
    return typeid (bppFS).hash_code ();
  }

  NodeRef ConfiguredSimplex::recreate (Context & c, NodeRefVec && deps) {
    auto m = ConfiguredParametrizable::createConfigured<Target, Self> (c, std::move (deps), std::unique_ptr<Target>(dynamic_cast<Target*>(simplex_->clone ())));
    m->config = this->config; // Duplicate derivation config
    return m;
  }

/***********************************************************/
  
  // FrequenciesFromSimplex

  FrequenciesFromSimplex::FrequenciesFromSimplex (
    NodeRefVec && deps, const Dimension<Eigen::RowVectorXd> & dim)
    : Value<Eigen::RowVectorXd> (std::move (deps)), targetDimension_ (dim) {}

  std::string FrequenciesFromSimplex::debugInfo () const {
    using namespace numeric;
    return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
  }

  // FrequenciesFromSimplex additional arguments = ().
  bool FrequenciesFromSimplex::compareAdditionalArguments (const Node_DF & other) const {
    return dynamic_cast<const Self *> (&other) != nullptr;
  }

  NodeRef FrequenciesFromSimplex::derive (Context & c, const Node_DF & node) {
    // d(equFreqs)/dn = sum_i d(equFreqs)/dx_i * dx_i/dn (x_i = simplex parameters)
    auto simplexDep = this->dependency (0);
    auto & simplex = static_cast<ConfiguredSimplex &> (*simplexDep);
    auto buildFWithNewSimplex = [this, &c](NodeRef && newSimplex) {
      return ConfiguredParametrizable::createRowVector<ConfiguredSimplex, Self> (c, {std::move (newSimplex)}, targetDimension_);
    };
      
    NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<ConfiguredSimplex, T > (
      c, simplex, node, targetDimension_, buildFWithNewSimplex);
    
    return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
  }

  NodeRef FrequenciesFromSimplex::recreate (Context & c, NodeRefVec && deps) {
    return ConfiguredParametrizable::createRowVector<ConfiguredSimplex, Self> (c, std::move (deps), targetDimension_);
  }

  void FrequenciesFromSimplex::compute () {
    const auto * simplex = accessValueConstCast<const Simplex *> (*this->dependency (0));
    const auto & freqsFromFS = simplex->getFrequencies ();
    auto & r = this->accessValueMutable ();
    r = Eigen::Map<const T> (freqsFromFS.data(), static_cast<Eigen::Index> (freqsFromFS.size ()));
  }

  
} // namespace bpp
