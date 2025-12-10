// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Exceptions.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowCWiseComputing.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DiscreteDistribution.h>


using namespace std;

namespace bpp
{
ConfiguredDistribution::ConfiguredDistribution (Context& context, NodeRefVec&& deps, std::unique_ptr<DiscreteDistributionInterface>&& distrib) :
  Value<const DiscreteDistributionInterface*>(std::move(deps), distrib.get()),
  AbstractParametrizable(distrib->getNamespace()),
  config(),
  distrib_(std::move(distrib))
{
  for (const auto& dep:dependencies())
  {
    const auto& param = std::dynamic_pointer_cast<ConfiguredParameter>(dep);
    shareParameter_(param);
  }
}

ConfiguredDistribution::~ConfiguredDistribution () = default;

std::string ConfiguredDistribution::description () const { return "Distribution(" + distrib_->getName () + ")"; }

std::string ConfiguredDistribution::debugInfo () const
{
  return "nbClass=" + std::to_string (distrib_->getNumberOfCategories ());
}

// Model node additional arguments = (type of bpp::BranchModel).
// Everything else is determined by the node dependencies.
bool ConfiguredDistribution::compareAdditionalArguments (const Node_DF& other) const
{
  const auto* derived = dynamic_cast<const Self*>(&other);
  if (derived == nullptr)
  {
    return false;
  }
  else
  {
    const auto& thisDistrib = *distrib_;
    const auto& otherDistrib = *derived->distrib_;
    return typeid (thisDistrib) == typeid (otherDistrib);
  }
}

std::size_t ConfiguredDistribution::hashAdditionalArguments () const
{
  const auto& bppDistrib = *distrib_;
  return typeid (bppDistrib).hash_code ();
}

NodeRef ConfiguredDistribution::recreate(Context& c, NodeRefVec&& deps)
{
  auto m = ConfiguredParametrizable::createConfigured<Target, Self>(c, std::move (deps), std::unique_ptr<DiscreteDistributionInterface>{dynamic_cast<DiscreteDistributionInterface*>(distrib_->clone ())});
  m->config = this->config; // Duplicate derivation config
  return m;
}

//////////////////////////////////////////////
// ProbabilitiesFromDiscreteDistribution

ProbabilitiesFromDiscreteDistribution::ProbabilitiesFromDiscreteDistribution(
    NodeRefVec&& deps, const Dimension<Eigen::RowVectorXd>& dim)
  : Value<Eigen::RowVectorXd>(std::move (deps)), nbClass_ (dim) {}


std::string ProbabilitiesFromDiscreteDistribution::debugInfo () const
{
  using namespace numeric;
  return debug (this->accessValueConst ()) + " nbClass=" + to_string (nbClass_);
}

// ProbabilitiesFromDiscreteDistribution additional arguments = ().
bool ProbabilitiesFromDiscreteDistribution::compareAdditionalArguments(const Node_DF& other) const
{
  return dynamic_cast<const Self*>(&other) != nullptr;
}

NodeRef ProbabilitiesFromDiscreteDistribution::derive(Context& c, const Node_DF& node)
{
  // d(Prob)/dn = sum_i d(Prob)/dx_i * dx_i/dn (x_i = distrib parameters)
  auto distribDep = this->dependency (0);
  auto& distrib = static_cast<Dep&>(*distribDep);
  auto buildPWithNewDistrib = [this, &c](NodeRef&& newDistrib) {
        return ConfiguredParametrizable::createRowVector<Dep, Self>(c, {std::move(newDistrib)}, nbClass_);
      };

  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<ConfiguredDistribution, T >(
        c, distrib, node, nbClass_, buildPWithNewDistrib);
  return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), nbClass_);
}

NodeRef ProbabilitiesFromDiscreteDistribution::recreate (Context& c, NodeRefVec&& deps)
{
  return ConfiguredParametrizable::createRowVector<Dep, Self>(c, std::move (deps), nbClass_);
}

void ProbabilitiesFromDiscreteDistribution::compute()
{
  const auto* distrib = accessValueConstCast<const DiscreteDistributionInterface*>(*this->dependency (0));
  const auto& probasFromDistrib = distrib->getProbabilities();
  auto& r = this->accessValueMutable ();
  r = Eigen::Map<const T>(probasFromDistrib.data(), static_cast<Eigen::Index>(probasFromDistrib.size ()));
}


std::shared_ptr<ProbabilitiesFromDiscreteDistribution> ProbabilitiesFromDiscreteDistribution::create(Context& c, NodeRefVec&& deps)
{
  checkDependenciesNotNull (typeid (Self), deps);
  checkDependencyVectorSize (typeid (Self), deps, 1);
  checkNthDependencyIs<ConfiguredDistribution>(typeid (Self), deps, 0);
  size_t nbCat = accessValueConstCast<DiscreteDistributionInterface*>(*deps[0])->getNumberOfCategories();
  return cachedAs<ProbabilitiesFromDiscreteDistribution>(c, std::make_shared<ProbabilitiesFromDiscreteDistribution>(std::move(deps), RowVectorDimension(Eigen::Index(nbCat))));
}


////////////////////////////////////////////////////
// ProbabilityFromDiscreteDistribution

ProbabilityFromDiscreteDistribution::ProbabilityFromDiscreteDistribution (
    NodeRefVec&& deps, unsigned int nCat)
  : Value<double>(std::move (deps)), nCat_ (nCat) {}


std::string ProbabilityFromDiscreteDistribution::debugInfo () const
{
  using namespace numeric;
  return "proba=" + TextTools::toString(accessValueConst()) + ":nCat=" + TextTools::toString(nCat_);
}

// ProbabilityFromDiscreteDistribution additional arguments = ().
bool ProbabilityFromDiscreteDistribution::compareAdditionalArguments (const Node_DF& other) const
{
  const auto* derived = dynamic_cast<const Self*>(&other);
  return derived != nullptr && nCat_ == derived->nCat_;
}

std::shared_ptr<ProbabilityFromDiscreteDistribution> ProbabilityFromDiscreteDistribution::create (Context& c, NodeRefVec&& deps, unsigned int nCat)
{
  checkDependenciesNotNull (typeid (Self), deps);
  checkDependencyVectorSize (typeid (Self), deps, 1);
  checkNthDependencyIs<ConfiguredDistribution>(typeid (Self), deps, 0);
  return cachedAs<ProbabilityFromDiscreteDistribution>(c, std::make_shared<ProbabilityFromDiscreteDistribution>(std::move (deps), nCat));
}

NodeRef ProbabilityFromDiscreteDistribution::derive (Context& c, const Node_DF& node)
{
  // d(Prob)/dn = sum_i d(Prob)/dx_i * dx_i/dn (x_i = distrib parameters)
  auto distribDep = this->dependency (0);
  auto& distrib = static_cast<Dep&>(*distribDep);
  auto buildPWithNewDistrib = [this, &c](NodeRef&& newDistrib) {
        return this->create (c, {std::move (newDistrib)}, nCat_);
      };

  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T >(
        c, distrib, node, 1, buildPWithNewDistrib);
  return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), 1);
}

NodeRef ProbabilityFromDiscreteDistribution::recreate (Context& c, NodeRefVec&& deps)
{
  return ProbabilityFromDiscreteDistribution::create (c, {std::move (deps)}, nCat_);
}

void ProbabilityFromDiscreteDistribution::compute ()
{
  const auto* distrib = accessValueConstCast<const DiscreteDistributionInterface*>(*this->dependency (0));
  this->accessValueMutable () = distrib->getProbability((size_t)nCat_);
}

////////////////////////////////////////////////////
// CategoryFromDiscreteDistribution

CategoryFromDiscreteDistribution::CategoryFromDiscreteDistribution (
    NodeRefVec&& deps, unsigned int nCat)
  : Value<double>(std::move (deps)), nCat_ (nCat) {}


std::string CategoryFromDiscreteDistribution::debugInfo () const
{
  using namespace numeric;
  return "rate=" + TextTools::toString(accessValueConst()) + ":nCat=" + TextTools::toString(nCat_);
}

// CategoryFromDiscreteDistribution additional arguments = ().
bool CategoryFromDiscreteDistribution::compareAdditionalArguments (const Node_DF& other) const
{
  const auto* derived = dynamic_cast<const Self*>(&other);
  return derived != nullptr && nCat_ == derived->nCat_;
}

std::shared_ptr<CategoryFromDiscreteDistribution> CategoryFromDiscreteDistribution::create (Context& c, NodeRefVec&& deps, unsigned int nCat)
{
  checkDependenciesNotNull (typeid (Self), deps);
  checkDependencyVectorSize (typeid (Self), deps, 1);
  checkNthDependencyIs<ConfiguredDistribution>(typeid (Self), deps, 0);
  return cachedAs<CategoryFromDiscreteDistribution>(c, std::make_shared<CategoryFromDiscreteDistribution>(std::move (deps), nCat));
}

NodeRef CategoryFromDiscreteDistribution::derive (Context& c, const Node_DF& node)
{
  // d(Prob)/dn = sum_i d(Prob)/dx_i * dx_i/dn (x_i = distrib parameters)
  auto distribDep = this->dependency (0);
  auto& distrib = static_cast<Dep&>(*distribDep);
  auto buildPWithNewDistrib = [this, &c](NodeRef&& newDistrib) {
        return this->create (c, {std::move (newDistrib)}, nCat_);
      };

  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T >(
        c, distrib, node, 1, buildPWithNewDistrib);
  return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), 1);
}

NodeRef CategoryFromDiscreteDistribution::recreate (Context& c, NodeRefVec&& deps)
{
  return CategoryFromDiscreteDistribution::create (c, {std::move (deps)}, nCat_);
}

void CategoryFromDiscreteDistribution::compute ()
{
  const auto* distrib = accessValueConstCast<const DiscreteDistributionInterface*>(*this->dependency (0));
  double categoryFromDistrib = distrib->getCategory(nCat_);
  this->accessValueMutable () = categoryFromDistrib;
}
} // namespace bpp
