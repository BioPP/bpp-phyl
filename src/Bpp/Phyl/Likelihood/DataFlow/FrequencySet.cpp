// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Exceptions.h>
#include <Bpp/Phyl/Likelihood/DataFlow/FrequencySet.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Parametrizable.h>


using namespace std;

namespace bpp
{
// FrequencySet node

ConfiguredFrequencySet::ConfiguredFrequencySet (const Context& context, NodeRefVec&& deps, std::unique_ptr<FrequencySetInterface>&& freqset) :
  Value<const FrequencySetInterface*>(std::move(deps), freqset.get()),
  AbstractParametrizable(freqset->getNamespace()),
  config(),
  freqset_(std::move(freqset))
{
  for (const auto& dep:dependencies())
  {
    const auto& param = std::dynamic_pointer_cast<ConfiguredParameter>(dep);
    shareParameter_(param);
  }
}

ConfiguredFrequencySet::~ConfiguredFrequencySet () = default;

std::string ConfiguredFrequencySet::description () const { return "FreqSet(" + freqset_->getName () + ")"; }

std::string ConfiguredFrequencySet::debugInfo () const
{
  return "nbState=" + std::to_string (freqset_->getAlphabet ()->getSize ());
}

// FrequencySet node additional arguments = (type of bpp::FrequencySet).
// Everything else is determined by the node dependencies.

bool ConfiguredFrequencySet::compareAdditionalArguments (const Node_DF& other) const
{
  const auto* derived = dynamic_cast<const Self*>(&other);
  if (derived == nullptr)
  {
    return false;
  }
  else
  {
    const auto& thisFS = *freqset_;
    const auto& otherFS = *derived->freqset_;
    return typeid (thisFS) == typeid (otherFS);
  }
}

std::size_t ConfiguredFrequencySet::hashAdditionalArguments () const
{
  const auto& bppFS = *freqset_;
  return typeid (bppFS).hash_code ();
}


NodeRef ConfiguredFrequencySet::recreate (Context& c, NodeRefVec&& deps)
{
  auto m = ConfiguredParametrizable::createConfigured<Target, Self>(c, std::move (deps), std::unique_ptr<Target>{dynamic_cast<Target*>(freqset_->clone ())});
  m->config = this->config; // Duplicate derivation config
  return m;
}

// FrequenciesFromFrequencySet

FrequenciesFromFrequencySet::FrequenciesFromFrequencySet (
    NodeRefVec&& deps, const Dimension<Eigen::RowVectorXd>& dim)
  : Value<Eigen::RowVectorXd>(std::move (deps)), targetDimension_ (dim) {}

std::string FrequenciesFromFrequencySet::debugInfo () const
{
  using namespace numeric;
  return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
}

// FrequenciesFromFrequencySet additional arguments = ().
bool FrequenciesFromFrequencySet::compareAdditionalArguments (const Node_DF& other) const
{
  return dynamic_cast<const Self*>(&other) != nullptr;
}

NodeRef FrequenciesFromFrequencySet::derive (Context& c, const Node_DF& node)
{
  // d(equFreqs)/dn = sum_i d(equFreqs)/dx_i * dx_i/dn (x_i = freqset parameters)
  auto freqSetDep = this->dependency (0);
  auto& freqset = static_cast<ConfiguredFrequencySet&>(*freqSetDep);
  auto buildFWithNewFreqSet = [this, &c](NodeRef&& newFreqSet) {
        return ConfiguredParametrizable::createRowVector<ConfiguredFrequencySet, Self>(c, {std::move (newFreqSet)}, targetDimension_);
      };

  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<ConfiguredFrequencySet, T >(
        c, freqset, node, targetDimension_, buildFWithNewFreqSet);
  return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
}

NodeRef FrequenciesFromFrequencySet::recreate (Context& c, NodeRefVec&& deps)
{
  return ConfiguredParametrizable::createRowVector<ConfiguredFrequencySet, Self>(c, std::move (deps), targetDimension_);
}

void FrequenciesFromFrequencySet::compute ()
{
  const auto* freqset = accessValueConstCast<const FrequencySetInterface*>(*this->dependency (0));
  const auto& freqsFromFS = freqset->getFrequencies ();
  auto& r = this->accessValueMutable ();
  r = Eigen::Map<const T>(freqsFromFS.data(), static_cast<Eigen::Index>(freqsFromFS.size ()));
}
} // namespace bpp
