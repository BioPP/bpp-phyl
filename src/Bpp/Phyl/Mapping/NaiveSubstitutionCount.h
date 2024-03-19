// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_NAIVESUBSTITUTIONCOUNT_H
#define BPP_PHYL_MAPPING_NAIVESUBSTITUTIONCOUNT_H

#include <Bpp/Numeric/Matrix/Matrix.h>

#include "WeightedSubstitutionCount.h"

namespace bpp
{
/**
 * @brief Naive substitution count.
 *
 * This substitution count is defined as follow:
 * - 0 if @f$i = j@f$,
 * - 1 if @f$i \neq j @f$.
 *
 * Reference (for instance):
 * TuffÃÂ©ry P, Darlu P.
 * Exploring a phylogenetic approach for the detection of correlated substitutions in proteins.
 * Mol Biol Evol. 2000 Nov;17(11):1753-9
 *
 * @author Julien Dutheil
 */
class NaiveSubstitutionCount :
  public AbstractSubstitutionCount,
  public AbstractWeightedSubstitutionCount
{
private:
  bool allowSelf_;
  std::vector<int> supportedChars_;

public:
  /**
   * @brief Build a new simple substitution count.
   *
   * @param model The substitution model for which this substitution count is parametrized.
   *              The model is not used in the calculation, only for specifying the modeled states.
   * @param reg A pointer toward a substitution register object which discribes the type of substitutions to map.
   * @param allowSelf Tells if "self" mutations, from X to X should be counted together with the ones of type X to Y where X and Y are in the same category, if relevent.
   * The default is "no", to be consistent with other types of substitution counts which account for multiple substitutions, in which case it does not make sense to count "X to X".
   * @param weights the weights of the counts
   */
  NaiveSubstitutionCount(
      std::shared_ptr<const SubstitutionModelInterface> model,
      std::shared_ptr<const SubstitutionRegisterInterface> reg,
      bool allowSelf = false,
      std::shared_ptr<const AlphabetIndex2> weights = nullptr) :
    AbstractSubstitutionCount(reg),
    AbstractWeightedSubstitutionCount(weights),
    allowSelf_(allowSelf),
    supportedChars_(model->getAlphabetStates()) {}

  NaiveSubstitutionCount(
      std::shared_ptr<const StateMapInterface> stateMap,
      std::shared_ptr<const SubstitutionRegisterInterface> reg,
      bool allowSelf = false,
      std::shared_ptr<const AlphabetIndex2> weights = nullptr) :
    AbstractSubstitutionCount(reg),
    AbstractWeightedSubstitutionCount(weights),
    allowSelf_(allowSelf),
    supportedChars_(stateMap->getAlphabetStates()) {}

  virtual ~NaiveSubstitutionCount() {}

  NaiveSubstitutionCount* clone() const override
  {
    return new NaiveSubstitutionCount(*this);
  }

public:
  double getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type = 1) const override
  {
    if (initialState == finalState && !allowSelf_)
      return 0;
    else
    {
      int alphabetState1 = supportedChars_[initialState];
      int alphabetState2 = supportedChars_[finalState];
      return register_->getType(initialState, finalState) == type ? (weights_ ? weights_->getIndex(alphabetState1, alphabetState2) : 1.) : 0.;
    }
  }

  std::unique_ptr< Matrix<double>> getAllNumbersOfSubstitutions(double length, size_t type = 1) const override;

  void storeAllNumbersOfSubstitutions(double length, size_t type, Eigen::MatrixXd& mat) const override;

  std::vector<double> getNumberOfSubstitutionsPerType(size_t initialState, size_t finalState, double length) const override
  {
    std::vector<double> v(getNumberOfSubstitutionTypes());
    for (size_t t = 1; t <= getNumberOfSubstitutionTypes(); ++t)
    {
      v[t - 1] = getNumberOfSubstitutions(initialState, finalState, length, t);
    }
    return v;
  }

  void setSubstitutionModel(std::shared_ptr<const SubstitutionModelInterface> model) override
  {
    if (model)
      supportedChars_ = model->getAlphabetStates();
  }

private:
  void substitutionRegisterHasChanged() override {}
  void weightsHaveChanged() override {}
};

/**
 * @brief Labelling substitution count.
 *
 * This substitution count return a distinct number for each possible mutation.
 * - 0 if @f$i = j@f$,
 * - @f$a(i,j)@f$ if @f$i \neq j @f$, where 'a' is an index giving a unique value for each combination of i and j.
 */
class LabelSubstitutionCount :
  public AbstractSubstitutionCount
{
private:
  LinearMatrix<double> label_;
  std::vector<int> supportedChars_;

public:
  LabelSubstitutionCount(std::shared_ptr<const SubstitutionModelInterface> model);

  LabelSubstitutionCount(std::shared_ptr<const StateMapInterface> statemap);

  virtual ~LabelSubstitutionCount() {}

  LabelSubstitutionCount* clone() const override { return new LabelSubstitutionCount(*this); }

public:
  double getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type = 1) const override
  {
    return label_(initialState, finalState);
  }

  std::unique_ptr< Matrix<double>> getAllNumbersOfSubstitutions(double length, size_t type = 1) const override
  {
    return std::unique_ptr< Matrix<double>>(label_.clone());
  }

  void storeAllNumbersOfSubstitutions(double length, size_t type, Eigen::MatrixXd& mat) const override
  {
    auto nbStates = supportedChars_.size();

    mat.resize(Eigen::Index(nbStates), Eigen::Index(nbStates));

    for (size_t i = 0; i < nbStates; i++)
    {
      for (size_t j = 0; j < nbStates; j++)
      {
        mat(Eigen::Index(i), Eigen::Index(j)) = label_(i, j);
      }
    }
  }


  std::vector<double> getNumberOfSubstitutionsPerType(size_t initialState, size_t finalState, double length) const override
  {
    std::vector<double> v(1);
    v[0] = label_(initialState, finalState);
    return v;
  }

  void setSubstitutionModel(std::shared_ptr<const SubstitutionModelInterface> model) override
  {
    if (model)
      supportedChars_ = model->getAlphabetStates();
  }

  void setSubstitutionRegister(std::shared_ptr<const SubstitutionRegisterInterface> reg) override
  {
    throw Exception("OneJumpSubstitutionCount::setSubstitutionRegister. This SubstitutionsCount only works with a TotalSubstitutionRegister.");
  }

private:
  void substitutionRegisterHasChanged() override {}
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_NAIVESUBSTITUTIONCOUNT_H
