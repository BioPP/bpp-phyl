// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Exceptions.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Parametrizable.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Simplex_DF.h>


using namespace std;

namespace bpp
{
// Simplex node

ConfiguredSimplex::ConfiguredSimplex (const Context& context, NodeRefVec&& deps, std::unique_ptr<Simplex>&& simplex)
  : Value<const Simplex*>(std::move (deps), simplex.get ()), AbstractParametrizable(simplex->getNamespace())// , context_(context)
  , simplex_(std::move(simplex))
{
  for (const auto& dep:dependencies())
  {
    const auto& param = std::dynamic_pointer_cast<ConfiguredParameter>(dep);
    shareParameter_(param);
  }
}

ConfiguredSimplex::~ConfiguredSimplex () = default;

std::string ConfiguredSimplex::description () const { return "Simplex"; }

std::string ConfiguredSimplex::debugInfo () const
{
  return "nbState=" + std::to_string (simplex_->dimension ());
}

// Simplex node additional arguments = (type of bpp::Simplex).
// Everything else is determined by the node dependencies.

bool ConfiguredSimplex::compareAdditionalArguments (const Node_DF& other) const
{
  const auto* derived = dynamic_cast<const Self*>(&other);
  if (derived == nullptr)
  {
    return false;
  }
  else
  {
    const auto& thisFS = *simplex_;
    const auto& otherFS = *derived->simplex_;
    return typeid (thisFS) == typeid (otherFS);
  }
}

std::size_t ConfiguredSimplex::hashAdditionalArguments () const
{
  const auto& bppFS = *simplex_;
  return typeid (bppFS).hash_code ();
}

NodeRef ConfiguredSimplex::recreate (Context& c, NodeRefVec&& deps)
{
  auto m = ConfiguredParametrizable::createConfigured<Target, Self>(c, std::move (deps), std::unique_ptr<Target>(dynamic_cast<Target*>(simplex_->clone ())));
  m->config = this->config; // Duplicate derivation config
  return m;
}

/***********************************************************/

// FrequenciesFromSimplex

FrequenciesFromSimplex::FrequenciesFromSimplex (
  NodeRefVec&& deps, const Dimension<Eigen::RowVectorXd>& dim)
  : Value<Eigen::RowVectorXd>(std::move (deps)), targetDimension_ (dim) {}

std::string FrequenciesFromSimplex::debugInfo () const
{
  using namespace numeric;
  return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
}

// FrequenciesFromSimplex additional arguments = ().
bool FrequenciesFromSimplex::compareAdditionalArguments (const Node_DF& other) const
{
  return dynamic_cast<const Self*>(&other) != nullptr;
}

NodeRef FrequenciesFromSimplex::derive (Context& c, const Node_DF& node)
{
  // d(equFreqs)/dn = sum_i d(equFreqs)/dx_i * dx_i/dn (x_i = simplex parameters)
  auto simplexDep = this->dependency (0);
  auto& simplex = static_cast<ConfiguredSimplex&>(*simplexDep);
  auto buildFWithNewSimplex = [this, &c](NodeRef&& newSimplex) {
                                return ConfiguredParametrizable::createRowVector<ConfiguredSimplex, Self>(c, {std::move (newSimplex)}, targetDimension_);
                              };

  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<ConfiguredSimplex, T >(
    c, simplex, node, targetDimension_, buildFWithNewSimplex);

  return CWiseAdd<T, ReductionOf<T> >::create (c, std::move (derivativeSumDeps), targetDimension_);
}

NodeRef FrequenciesFromSimplex::recreate (Context& c, NodeRefVec&& deps)
{
  return ConfiguredParametrizable::createRowVector<ConfiguredSimplex, Self>(c, std::move (deps), targetDimension_);
}

void FrequenciesFromSimplex::compute ()
{
  const auto* simplex = accessValueConstCast<const Simplex*>(*this->dependency (0));
  const auto& freqsFromFS = simplex->getFrequencies ();
  auto& r = this->accessValueMutable ();
  r = Eigen::Map<const T>(freqsFromFS.data(), static_cast<Eigen::Index>(freqsFromFS.size ()));
}
} // namespace bpp
