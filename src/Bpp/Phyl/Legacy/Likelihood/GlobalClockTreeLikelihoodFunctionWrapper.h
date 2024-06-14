// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LEGACY_LIKELIHOOD_GLOBALCLOCKTREELIKELIHOODFUNCTIONWRAPPER_H
#define BPP_PHYL_LEGACY_LIKELIHOOD_GLOBALCLOCKTREELIKELIHOODFUNCTIONWRAPPER_H


#include "TreeLikelihood.h"

namespace bpp
{
class GlobalClockTreeLikelihoodFunctionWrapper :
  public virtual SecondOrderDerivable,
  public AbstractParametrizable
{
private:
  std::shared_ptr<TreeLikelihoodInterface> tl_;

public:
  GlobalClockTreeLikelihoodFunctionWrapper(std::shared_ptr<TreeLikelihoodInterface> tl) :
    AbstractParametrizable(""),
    tl_(tl)
  {
    initParameters_();
  }

  GlobalClockTreeLikelihoodFunctionWrapper(const GlobalClockTreeLikelihoodFunctionWrapper& gctlfw) :
    AbstractParametrizable(gctlfw), tl_(gctlfw.tl_)
  {}

  GlobalClockTreeLikelihoodFunctionWrapper& operator=(const GlobalClockTreeLikelihoodFunctionWrapper& gctlfw)
  {
    AbstractParametrizable::operator=(gctlfw);
    tl_ = gctlfw.tl_;
    return *this;
  }

  GlobalClockTreeLikelihoodFunctionWrapper* clone() const { return new GlobalClockTreeLikelihoodFunctionWrapper(*this); }

public:
  void setParameters(const ParameterList& pl)
  {
    // For now we go the hard way and recompute everything:
    matchParametersValues(pl);
  }

  double getValue() const { return tl_->getValue(); }

  void fireParameterChanged(const bpp::ParameterList& pl);

  void enableSecondOrderDerivatives(bool yn) { tl_->enableSecondOrderDerivatives(yn); }
  bool enableSecondOrderDerivatives() const { return tl_->enableSecondOrderDerivatives(); }
  void enableFirstOrderDerivatives(bool yn) { tl_->enableFirstOrderDerivatives(yn); }
  bool enableFirstOrderDerivatives() const { return tl_->enableFirstOrderDerivatives(); }
  double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const { return tl_->getSecondOrderDerivative(variable1, variable2); }
  double getSecondOrderDerivative(const std::string& variable) const { return tl_->getSecondOrderDerivative(variable); }
  double getFirstOrderDerivative(const std::string& variable) const { return tl_->getFirstOrderDerivative(variable); }

  ParameterList getHeightParameters() const;

private:
  void initParameters_();
  void computeBranchLengthsFromHeights_(const Node* node, double height, ParameterList& brlenPl);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_LIKELIHOOD_GLOBALCLOCKTREELIKELIHOODFUNCTIONWRAPPER_H
