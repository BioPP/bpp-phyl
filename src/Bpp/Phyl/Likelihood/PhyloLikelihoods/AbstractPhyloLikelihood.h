// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_ABSTRACTPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_ABSTRACTPHYLOLIKELIHOOD_H


#include "../DataFlow/DataFlowNumeric.h"
#include "../DataFlow/LikelihoodCalculation.h"
#include "../DataFlow/Parametrizable.h"
#include "PhyloLikelihood.h"

namespace bpp
{
class LikelihoodCalculation;

class AbstractPhyloLikelihood :
  public virtual PhyloLikelihoodInterface
{
public:
  // Cache generated nodes representing derivatives, to avoid recreating them every time.
  // Using the mutable keyword because the table must be changed even in const methods.
  struct StringPairHash
  {
    std::size_t operator()(const std::pair<std::string, std::string>& p) const
    {
      std::hash<std::string> strHash{};
      return strHash (p.first) ^ (strHash (p.second) << 1);
    }
  };

protected:
  Context& context_;

  /**
   * @brief the value
   */
  mutable DataLik minusLogLik_;

  /**
   * @brief For Dataflow computing
   */
  mutable std::unordered_map<std::string, ValueRef<DataLik> > firstOrderDerivativeNodes_;

  mutable std::unordered_map<std::pair<std::string, std::string>, ValueRef<DataLik>,
                             StringPairHash>
  secondOrderDerivativeNodes_;

public:
  AbstractPhyloLikelihood(Context& context) :
    context_(context),
    minusLogLik_(0)
  {}

  AbstractPhyloLikelihood(const AbstractPhyloLikelihood& apl) :
    context_(apl.context_),
    minusLogLik_(apl.minusLogLik_)
  {
    shareParameters(apl.getParameters());
  }

  AbstractPhyloLikelihood& operator=(const AbstractPhyloLikelihood& apl)
  {
    context_ = apl.context_;
    minusLogLik_ = apl.minusLogLik_;
    shareParameters(apl.getParameters());
    return *this;
  }

  virtual ~AbstractPhyloLikelihood() {}

  const Context& context() const override { return context_; }
  
  Context& context() override { return context_; }

  /**
   * @brief Sets the computeLikelihoods_ to true.
   *
   */
  virtual bool isInitialized() const override
  {
    return false;
  }

public:

  /**
   * @brief Share Parameters, that are DF_parameters
   */
  void shareParameters(const ParameterList& variableNodes)
  {
    this->getParameters_().shareParameters(variableNodes);
  }

  void setParameters(const ParameterList& parameters) override
  {
    setParametersValues(parameters);
  }

  ValueRef<DataLik> getLikelihoodNode() const override
  {
    return getLikelihoodCalculation()->getLikelihoodNode();
  }


  // bpp::Function

  /**
   * @brief Tell if derivatives must be computed: for Function
   * inheritance.
   *
   */
  virtual void enableFirstOrderDerivatives(bool yn)  override {}
  virtual void enableSecondOrderDerivatives(bool yn)  override {}
  bool enableFirstOrderDerivatives() const override { return true; }
  bool enableSecondOrderDerivatives() const override { return true; }

  /*
   * @brief return the value, ie -loglikelihood
   *
   * !!! check on computeLikelihoods_ is not done here.
   *
   */
  double getValue() const override
  {
    minusLogLik_ = -getLikelihoodNode()->targetValue();
    return convert(minusLogLik_);
  }

  // bpp::DerivableFirstOrder
  double getFirstOrderDerivative (const std::string& variable) const override
  {
    using namespace std; // for isinf
    using namespace numeric; // for isinf
    return convert(-firstOrderDerivativeNode (variable)->targetValue ());
  }

  // Get nodes of derivatives directly
  ValueRef<DataLik> firstOrderDerivativeNode (const std::string& variable) const
  {
    const auto it = firstOrderDerivativeNodes_.find (variable);
    if (it != firstOrderDerivativeNodes_.end ())
    {
      return it->second;
    }
    else
    {
      auto node = getLikelihoodNode()->deriveAsValue (context_, accessVariableNode (variable));
      firstOrderDerivativeNodes_.emplace (variable, node);
      return node;
    }
  }

  // bpp::DerivableSecondOrder
  double getSecondOrderDerivative (const std::string& variable) const override
  {
    return getSecondOrderDerivative (variable, variable);
  }

  double getSecondOrderDerivative (const std::string& variable1,
                                   const std::string& variable2) const override
  {
    using namespace std; // for isinf
    using namespace numeric; // for isinf
    return convert(-secondOrderDerivativeNode (variable1, variable2)->targetValue ());
  }

  ValueRef<DataLik> secondOrderDerivativeNode (const std::string& variable1,
                                               const std::string& variable2) const
  {
    const auto key = std::make_pair (variable1, variable2);
    const auto it = secondOrderDerivativeNodes_.find (key);
    if (it != secondOrderDerivativeNodes_.end ())
    {
      return it->second;
    }
    else
    {
      // Reuse firstOrderDerivative() to generate the first derivative with caching
      auto node =
        firstOrderDerivativeNode (variable1)->deriveAsValue (context_, accessVariableNode (variable2));
      secondOrderDerivativeNodes_.emplace (key, node);
      return node;
    }
  }

protected:
  static Node_DF& accessVariableNode (const Parameter& param)
  {
    return *dynamic_cast<const ConfiguredParameter&>(param).dependency(0);
  }

  Node_DF& accessVariableNode (const std::string& name) const
  {
    return accessVariableNode (parameter(name));
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_ABSTRACTPHYLOLIKELIHOOD_H
