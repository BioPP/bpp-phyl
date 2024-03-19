// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_ONEJUMPSUBSTITUTIONCOUNT_H
#define BPP_PHYL_MAPPING_ONEJUMPSUBSTITUTIONCOUNT_H


#include "../Model/SubstitutionModel.h"
#include "SubstitutionCount.h"

namespace bpp
{
/**
 * @brief Computes the probability that at least one jump occured on a branch, given the initial and final state.
 *
 * This probability is defined as
 * @f[
 * p_{x,y}(l) = \left\{\begin{array}{ll}1 & \mathrm{if} x \neq y \\ 1 - \exp{\left(Q \cdot t\right)}_{x,y} & \mathrm{otherwise.}\end{array}\right.
 * @f]
 *
 * @author Julien Dutheil
 */
class OneJumpSubstitutionCount :
  public AbstractSubstitutionCount
{
private:
  std::shared_ptr<const SubstitutionModelInterface> model_;
  mutable RowMatrix<double> tmp_;

public:
  OneJumpSubstitutionCount(std::shared_ptr<const SubstitutionModelInterface> model) :
    AbstractSubstitutionCount(std::make_shared<const TotalSubstitutionRegister>(model->getStateMap())),
    model_(model), tmp_() {}

  OneJumpSubstitutionCount(std::shared_ptr<const StateMapInterface> statemap) :
    AbstractSubstitutionCount(std::make_shared<const TotalSubstitutionRegister>(statemap)),
    model_(nullptr), tmp_() {}

  OneJumpSubstitutionCount(const OneJumpSubstitutionCount& ojsc) :
    AbstractSubstitutionCount(ojsc),
    model_(ojsc.model_), tmp_(ojsc.tmp_) {}

  OneJumpSubstitutionCount& operator=(const OneJumpSubstitutionCount& ojsc)
  {
    AbstractSubstitutionCount::operator=(ojsc);
    model_    = ojsc.model_;
    tmp_      = ojsc.tmp_;
    return *this;
  }

  virtual ~OneJumpSubstitutionCount() {}

  virtual OneJumpSubstitutionCount* clone() const override { return new OneJumpSubstitutionCount(*this); }

public:
  double getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type = 1) const override
  {
    if (!model_)
      throw Exception("OneJumpSubstitutionCount::getNumberOfSubstitutions: model not defined.");

    if (finalState != initialState) return 1.;
    else return 1. - model_->Pij_t(initialState, finalState, length);
  }

  std::unique_ptr< Matrix<double>> getAllNumbersOfSubstitutions(double length, size_t type = 1) const override;

  void storeAllNumbersOfSubstitutions(double length, size_t type, Eigen::MatrixXd& mat) const override;

  std::vector<double> getNumberOfSubstitutionsPerType(size_t initialState, size_t finalState, double length) const override
  {
    std::vector<double> v(0);
    v[0] = getNumberOfSubstitutions(initialState, finalState, length, 0);
    return v;
  }

  void setSubstitutionModel(std::shared_ptr<const SubstitutionModelInterface> model) override { model_ = model; }

  /*
   *@param reg pointer to a SubstitutionRegister
   *
   */
  void setSubstitutionRegister(std::shared_ptr<const SubstitutionRegisterInterface> reg) override
  {
    throw Exception("OneJumpSubstitutionCount::setSubstitutionRegister. This SubstitutionsCount only works with a TotalSubstitutionRegister.");
  }

private:
  void substitutionRegisterHasChanged() override {}
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_ONEJUMPSUBSTITUTIONCOUNT_H
