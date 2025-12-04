// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Exceptions.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowCWiseComputing.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Parametrizable.h>
#include <Bpp/Phyl/Likelihood/DataFlow/TransitionMatrix.h>


using namespace std;
using namespace bpp;

// TransitionMatrix node

ConfiguredTransitionMatrix::ConfiguredTransitionMatrix (Context& context, NodeRefVec&& deps, std::unique_ptr<HmmTransitionMatrix>&& hmm):
  Value<const HmmTransitionMatrix*>(std::move (deps), hmm.get ()), AbstractParametrizable(hmm->getNamespace()),
  config(),
  hmm_(std::move(hmm))
{
  for (const auto& dep:dependencies())
  {
    shareParameter_(std::dynamic_pointer_cast<ConfiguredParameter>(dep));
  }
}

ConfiguredTransitionMatrix::~ConfiguredTransitionMatrix () = default;

std::string ConfiguredTransitionMatrix::description () const { return "TransitionMatrix(HMM)"; }

std::string ConfiguredTransitionMatrix::debugInfo () const
{
  return "nbState=" + std::to_string (hmm_->getNumberOfStates ());
}

// TransitionMatrix node additional arguments = (type of BranchTransitionMatrix).
// Everything else is determined by the node dependencies.
bool ConfiguredTransitionMatrix::compareAdditionalArguments (const Node_DF& other) const
{
  const auto* derived = dynamic_cast<const Self*>(&other);
  if (derived == nullptr)
  {
    return false;
  }
  else
  {
    const auto& thisHMM = *hmm_;
    const auto& otherHMM = *derived->hmm_;
    return typeid (thisHMM) == typeid (otherHMM);
  }
}

std::size_t ConfiguredTransitionMatrix::hashAdditionalArguments () const
{
  const auto& bppTransitionMatrix = *hmm_;
  return typeid (bppTransitionMatrix).hash_code ();
}

NodeRef ConfiguredTransitionMatrix::recreate (Context& c, NodeRefVec&& deps)
{
  auto m = ConfiguredParametrizable::createConfigured<Target, Self>(c, std::move (deps), std::unique_ptr<Target>(dynamic_cast<Target*>(hmm_->clone ())));
  m->config = this->config; // Duplicate derivation config
  return m;
}

// EquilibriumFrequenciesFromTransitionMatrix

EquilibriumFrequenciesFromTransitionMatrix::EquilibriumFrequenciesFromTransitionMatrix (
    NodeRefVec&& deps, const Dimension<T>& dim)
  : Value<T>(std::move (deps)), targetDimension_ (dim) {}

std::string EquilibriumFrequenciesFromTransitionMatrix::debugInfo () const
{
  using namespace numeric;
  return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
}

// EquilibriumFrequenciesFromTransitionMatrix additional arguments = ().
bool EquilibriumFrequenciesFromTransitionMatrix::compareAdditionalArguments (const Node_DF& other) const
{
  return dynamic_cast<const Self*>(&other) != nullptr;
}

NodeRef EquilibriumFrequenciesFromTransitionMatrix::derive (Context& c, const Node_DF& node)
{
  // d(equFreqs)/dn = sum_i d(equFreqs)/dx_i * dx_i/dn (x_i = model parameters)
  auto hmmDep = this->dependency (0);
  auto& hmm = static_cast<Dep&>(*hmmDep);
  auto buildFWithNewTransitionMatrix = [this, &c](NodeRef&& newTransitionMatrix) {
        return ConfiguredParametrizable::createVector<Dep, Self>(c, {std::move (newTransitionMatrix)}, targetDimension_);
      };

  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T>(
        c, hmm, node, targetDimension_, buildFWithNewTransitionMatrix);

  return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
}

NodeRef EquilibriumFrequenciesFromTransitionMatrix::recreate (Context& c, NodeRefVec&& deps)
{
  return ConfiguredParametrizable::createVector<Dep, Self>(c, std::move (deps), targetDimension_);
}

void EquilibriumFrequenciesFromTransitionMatrix::compute ()
{
  const auto phmm = dynamic_cast<const HmmTransitionMatrix*>(accessValueConstCast<const HmmTransitionMatrix*>(*this->dependency (0)));

  const Vdouble& freqs = phmm->getEquilibriumFrequencies();

  auto& r = this->accessValueMutable ();
  r = Eigen::Map<const T>(freqs.data (), static_cast<Eigen::Index>(freqs.size ()));
}

// TransitionMatrixFromTransitionMatrix

TransitionMatrixFromTransitionMatrix::TransitionMatrixFromTransitionMatrix (NodeRefVec&& deps,
    const Dimension<Eigen::MatrixXd>& dim)
  : Value<Eigen::MatrixXd>(std::move (deps)), targetDimension_ (dim) {}

std::string TransitionMatrixFromTransitionMatrix::debugInfo () const
{
  using namespace numeric;
  return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
}

// TransitionMatrixFromTransitionMatrix additional arguments = ().
bool TransitionMatrixFromTransitionMatrix::compareAdditionalArguments (const Node_DF& other) const
{
  return dynamic_cast<const Self*>(&other) != nullptr;
}

NodeRef TransitionMatrixFromTransitionMatrix::derive (Context& c, const Node_DF& node)
{
  // dtm/dn = sum_i dtm/dx_i * dx_i/dn + dtm/dbrlen * dbrlen/dn (x_i = model parameters).
  auto hmmDep = this->dependency (0);

  // TransitionMatrix part
  auto& hmm = static_cast<Dep&>(*hmmDep);
  auto buildFWithNewTransitionMatrix = [this, &c](NodeRef&& newTransitionMatrix) {
        return ConfiguredParametrizable::createMatrix<Dep, Self>(c, {std::move (newTransitionMatrix)}, targetDimension_);
      };

  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T>(
        c, hmm, node, targetDimension_, buildFWithNewTransitionMatrix);

  return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
}

NodeRef TransitionMatrixFromTransitionMatrix::recreate (Context& c, NodeRefVec&& deps)
{
  return ConfiguredParametrizable::createMatrix<Dep, Self>(c, std::move (deps), targetDimension_);
}

void TransitionMatrixFromTransitionMatrix::compute ()
{
  const auto phmm = dynamic_cast<const HmmTransitionMatrix*>(accessValueConstCast<const HmmTransitionMatrix*>(*this->dependency (0)));

  auto& r = this->accessValueMutable ();
  copyBppToEigen (phmm->getPij(), r);
}
