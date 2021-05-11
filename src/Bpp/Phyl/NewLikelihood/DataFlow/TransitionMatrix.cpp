//
// File: TransitionMatrix.cpp
// Authors:Laurent Guéguen
// Created: vendredi 3 juillet 2020, à 17h 54
// Last modified: vendredi 3 juillet 2020, à 17h 54
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
#include <Bpp/Phyl/NewLikelihood/DataFlow/TransitionMatrix.h>

#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowCWiseComputing.h>

#include <Bpp/Phyl/NewLikelihood/DataFlow/Parametrizable.h>


using namespace std;
using namespace bpp;

// TransitionMatrix node

ConfiguredTransitionMatrix::ConfiguredTransitionMatrix (Context& context, NodeRefVec && deps, std::unique_ptr<HmmTransitionMatrix> && hmm)
  : Value<const HmmTransitionMatrix*> (std::move (deps), hmm.get ()), AbstractParametrizable(hmm->getNamespace())// , context_(context)
  , hmm_(std::move(hmm))
{
  for (const auto& dep:dependencies())
    shareParameter_(std::dynamic_pointer_cast<ConfiguredParameter>(dep));
}

ConfiguredTransitionMatrix::~ConfiguredTransitionMatrix () = default;

std::string ConfiguredTransitionMatrix::description () const { return "TransitionMatrix(HMM)"; }

std::string ConfiguredTransitionMatrix::debugInfo () const {
  return "nbState=" + std::to_string (hmm_->getNumberOfStates ());
}

// TransitionMatrix node additional arguments = (type of BranchTransitionMatrix).
// Everything else is determined by the node dependencies.
bool ConfiguredTransitionMatrix::compareAdditionalArguments (const Node_DF & other) const {
  const auto * derived = dynamic_cast<const Self *> (&other);
  if (derived == nullptr) {
    return false;
  } else {
    const auto & thisHMM = *hmm_;
    const auto & otherHMM = *derived->hmm_;
    return typeid (thisHMM) == typeid (otherHMM);
  }
}
    
std::size_t ConfiguredTransitionMatrix::hashAdditionalArguments () const {
  const auto & bppTransitionMatrix = *hmm_;
  return typeid (bppTransitionMatrix).hash_code ();
}

NodeRef ConfiguredTransitionMatrix::recreate (Context & c, NodeRefVec && deps) {
  auto m = ConfiguredParametrizable::createConfigured<Target, Self> (c, std::move (deps), std::unique_ptr<Target>(dynamic_cast<Target*>(hmm_->clone ())));
  m->config = this->config; // Duplicate derivation config
  return m;
}
    
// EquilibriumFrequenciesFromTransitionMatrix

EquilibriumFrequenciesFromTransitionMatrix::EquilibriumFrequenciesFromTransitionMatrix (
  NodeRefVec && deps, const Dimension<T> & dim)
  : Value<T> (std::move (deps)), targetDimension_ (dim) {}

std::string EquilibriumFrequenciesFromTransitionMatrix::debugInfo () const {
  using namespace numeric;
  return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
}

// EquilibriumFrequenciesFromTransitionMatrix additional arguments = ().
bool EquilibriumFrequenciesFromTransitionMatrix::compareAdditionalArguments (const Node_DF & other) const {
  return dynamic_cast<const Self *> (&other) != nullptr;
}

NodeRef EquilibriumFrequenciesFromTransitionMatrix::derive (Context & c, const Node_DF & node) {
  // d(equFreqs)/dn = sum_i d(equFreqs)/dx_i * dx_i/dn (x_i = model parameters)
  auto hmmDep = this->dependency (0);
  auto & hmm = static_cast<Dep &> (*hmmDep);
  auto buildFWithNewTransitionMatrix = [this, &c](NodeRef && newTransitionMatrix) {
    return ConfiguredParametrizable::createVector<Dep, Self> (c, {std::move (newTransitionMatrix)}, targetDimension_);
  };
      
  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T> (
    c, hmm, node, targetDimension_, buildFWithNewTransitionMatrix);
  
  return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
}

NodeRef EquilibriumFrequenciesFromTransitionMatrix::recreate (Context & c, NodeRefVec && deps) {
  return ConfiguredParametrizable::createVector<Dep, Self> (c, std::move (deps), targetDimension_);
}

void EquilibriumFrequenciesFromTransitionMatrix::compute () {

  const auto phmm = dynamic_cast<const HmmTransitionMatrix*> (accessValueConstCast<const HmmTransitionMatrix *> (*this->dependency (0)));

  const Vdouble& freqs = phmm->getEquilibriumFrequencies();
      
  auto & r = this->accessValueMutable ();
  r = Eigen::Map<const T> (freqs.data (), static_cast<Eigen::Index> (freqs.size ()));
}

// TransitionMatrixFromTransitionMatrix

TransitionMatrixFromTransitionMatrix::TransitionMatrixFromTransitionMatrix (NodeRefVec && deps,
                                                                            const Dimension<Eigen::MatrixXd> & dim)
  : Value<Eigen::MatrixXd> (std::move (deps)), targetDimension_ (dim) {}

std::string TransitionMatrixFromTransitionMatrix::debugInfo () const {
  using namespace numeric;
  return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
}

// TransitionMatrixFromTransitionMatrix additional arguments = ().
bool TransitionMatrixFromTransitionMatrix::compareAdditionalArguments (const Node_DF & other) const {
  return dynamic_cast<const Self *> (&other) != nullptr;
}

NodeRef TransitionMatrixFromTransitionMatrix::derive (Context & c, const Node_DF & node) {
  // dtm/dn = sum_i dtm/dx_i * dx_i/dn + dtm/dbrlen * dbrlen/dn (x_i = model parameters).
  auto hmmDep = this->dependency (0);

  // TransitionMatrix part
  auto & hmm = static_cast<Dep &> (*hmmDep);
  auto buildFWithNewTransitionMatrix = [this, &c](NodeRef && newTransitionMatrix) {
    return ConfiguredParametrizable::createMatrix<Dep, Self> (c, {std::move (newTransitionMatrix)}, targetDimension_);
  };
  
  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T> (
    c, hmm, node, targetDimension_, buildFWithNewTransitionMatrix);
  
  return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
}

NodeRef TransitionMatrixFromTransitionMatrix::recreate (Context & c, NodeRefVec && deps) {
  return ConfiguredParametrizable::createMatrix<Dep, Self> (c, std::move (deps), targetDimension_);
}

void TransitionMatrixFromTransitionMatrix::compute () {
  const auto phmm = dynamic_cast<const HmmTransitionMatrix*> (accessValueConstCast<const HmmTransitionMatrix *> (*this->dependency (0)));

  auto & r = this->accessValueMutable ();
  copyBppToEigen (phmm->getPij(), r);
}
