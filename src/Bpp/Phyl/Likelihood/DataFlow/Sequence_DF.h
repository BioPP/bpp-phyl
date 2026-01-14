// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_SEQUENCE_DF_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_SEQUENCE_DF_H

#include <Bpp/Exceptions.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowNumeric.h>
#include <functional>
#include <unordered_map>

#include "Definitions.h"

namespace bpp
{
/** 
 * @brief Data flow node representing a Sequence as a
 * Value<Eigen::MatrixXd> with a name.
 */
class Sequence_DF : public Value<MatrixLik>
{
private:
  std::string name_;

public:
  using Self = Sequence_DF;
  using T = MatrixLik;

  static ValueRef<T> create (Context& c, T&& value, const std::string& name)
  {
    return cachedAs<Self>(c, std::make_shared<Self>(std::move(value), name));
  }

  Sequence_DF (T&& value, const std::string& name) :
    Value<T>(NodeRefVec{}, std::move(value)),
    name_(name)
  {
    this->makeValid (); // Always valid
  }

  const std::string& getName() const
  {
    return name_;
  }

  std::string debugInfo () const final
  {
    using namespace numeric;
    return debug(this->accessValueConst ());
  }

  std::string description () const final
  {
    return Node_DF::description() + " " + name_;
  }

  bool compareAdditionalArguments (const Node_DF& other) const override
  {
    const auto* derived = dynamic_cast<const Self*>(&other);
    return derived != nullptr && name_ == derived->name_ && this->accessValueConst () == derived->accessValueConst ();
  }

  std::string color () const override
  {
    return "grey";
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    checkRecreateWithoutDependencies (typeid (Self), deps);
    return this->shared_from_this ();
  }

  std::size_t hashAdditionalArguments () const override
  {
    using namespace numeric;
    size_t seed = hash (this->accessValueConst ());
    combineHash<std::string>(seed, name_);
    return seed;
  }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized" // Remove STL warning
  NodeRef derive (Context& c, const Node_DF& node) final
  {
    const auto dim = Dimension<T>(this->accessValueConst ());
    if (&node == this)
    {
      return ConstantOne<T>::create (c, dim);
    }
    return ConstantZero<T>::create (c, dim);
  }
#pragma GCC diagnostic pop

private:
  void compute () final
  {
    // Constant is valid from construction
    failureComputeWasCalled (typeid (*this));
  }
};
} // namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_SEQUENCE_DF_H
