// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_LAPLACESUBSTITUTIONCOUNT_H
#define BPP_PHYL_MAPPING_LAPLACESUBSTITUTIONCOUNT_H


#include "../Model/SubstitutionModel.h"
#include "SubstitutionCount.h"

namespace bpp
{
/**
 * @brief Laplace estimate of the substitution count.
 *
 * This method uses Laplace transforms, as described in
 * Dutheil J, Pupko T, Jean-Marie A, Galtier N.
 * A model-based approach for detecting coevolving positions in a molecule.
 * Mol Biol Evol. 2005 Sep;22(9):1919-28.
 *
 * @see UniformizationSubstitutionCount
 * @see DecompositionSubstitutionCount
 * @author Julien Dutheil
 */
class LaplaceSubstitutionCount :
  public AbstractSubstitutionCount
{
private:
  std::shared_ptr<const SubstitutionModelInterface> model_;
  size_t cutOff_;
  mutable double currentLength_;
  mutable RowMatrix<double> m_;

public:
  LaplaceSubstitutionCount(
      std::shared_ptr<const SubstitutionModelInterface> model,
      size_t cutOff) :
    AbstractSubstitutionCount(std::make_shared<TotalSubstitutionRegister>(model->getStateMap())),
    model_        (model),
    cutOff_       (cutOff),
    currentLength_(0),
    m_            (model->getNumberOfStates(), model->getNumberOfStates())
  {}

  LaplaceSubstitutionCount(std::shared_ptr<const StateMapInterface> stateMap, size_t cutOff) :
    AbstractSubstitutionCount(std::make_shared<TotalSubstitutionRegister>(stateMap)),
    model_        (0),
    cutOff_       (cutOff),
    currentLength_(0),
    m_            (stateMap->getNumberOfModelStates(), stateMap->getNumberOfModelStates())
  {}

  LaplaceSubstitutionCount(const LaplaceSubstitutionCount& asc) :
    AbstractSubstitutionCount(asc),
    model_        (asc.model_),
    cutOff_       (asc.cutOff_),
    currentLength_(asc.currentLength_),
    m_            (asc.m_)
  {}

  LaplaceSubstitutionCount& operator=(const LaplaceSubstitutionCount& asc)
  {
    AbstractSubstitutionCount::operator=(asc);
    model_         = asc.model_;
    cutOff_        = asc.cutOff_;
    currentLength_ = asc.currentLength_;
    m_             = asc.m_;
    return *this;
  }

  virtual ~LaplaceSubstitutionCount() {}

  LaplaceSubstitutionCount* clone() const override { return new LaplaceSubstitutionCount(*this); }

public:
  double getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type = 1) const override;

  std::unique_ptr< Matrix<double>> getAllNumbersOfSubstitutions(double length, size_t type = 1) const override;

  void storeAllNumbersOfSubstitutions(double length, size_t type, Eigen::MatrixXd& mat) const override;

  std::vector<double> getNumberOfSubstitutionsPerType(size_t initialState, size_t finalState, double length) const override
  {
    std::vector<double> v(0);
    v[0] = getNumberOfSubstitutions(initialState, finalState, length, 0);
    return v;
  }

  void setSubstitutionModel(std::shared_ptr<const SubstitutionModelInterface> model) override;

  /*
   *@param reg pointer to a SubstitutionRegister
   *
   */
  void setSubstitutionRegister(std::shared_ptr<const SubstitutionRegisterInterface> reg) override
  {
    throw Exception("LaplaceSubstitutionCount::setSubstitutionRegister. This SubstitutionsCount only works with a TotalSubstitutionRegister.");
  }

protected:
  void computeCounts(double length) const;
  void substitutionRegisterHasChanged() override {}
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_LAPLACESUBSTITUTIONCOUNT_H
