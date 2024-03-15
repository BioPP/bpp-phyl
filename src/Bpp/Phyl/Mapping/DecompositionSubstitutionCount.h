// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_DECOMPOSITIONSUBSTITUTIONCOUNT_H
#define BPP_PHYL_MAPPING_DECOMPOSITIONSUBSTITUTIONCOUNT_H

#include <Bpp/Numeric/Matrix/Matrix.h>

#include "DecompositionMethods.h"
#include "SubstitutionDistance.h"
#include "WeightedSubstitutionCount.h"

namespace bpp
{
/**
 * @brief Analytical substitution count using the eigen decomposition method.
 *
 * The codes is adapted from the original R code by Paula Tataru and Asger Hobolth.
 *
 * @author Julien Dutheil
 */

class DecompositionSubstitutionCount :
  public AbstractSubstitutionCount,
  public AbstractWeightedSubstitutionCount,
  public AbstractSubstitutionDistance,
  public DecompositionMethods
{
private:
  mutable std::vector< RowMatrix<double> > counts_;
  mutable double currentLength_;

public:
  DecompositionSubstitutionCount(
      std::shared_ptr<const SubstitutionModelInterface> model,
      std::shared_ptr<const SubstitutionRegisterInterface> reg,
      std::shared_ptr<const AlphabetIndex2> weights = nullptr,
      std::shared_ptr<const AlphabetIndex2> distances = nullptr);

  DecompositionSubstitutionCount(
      std::shared_ptr<const SubstitutionRegisterInterface> reg,
      std::shared_ptr<const AlphabetIndex2> weights = nullptr,
      std::shared_ptr<const AlphabetIndex2> distances = nullptr);

  DecompositionSubstitutionCount(const DecompositionSubstitutionCount& dsc) :
    AbstractSubstitutionCount(dsc),
    AbstractWeightedSubstitutionCount(dsc),
    AbstractSubstitutionDistance(dsc),
    DecompositionMethods(dsc),
    counts_(dsc.counts_),
    currentLength_(dsc.currentLength_)
  {}

  DecompositionSubstitutionCount& operator=(const DecompositionSubstitutionCount& dsc)
  {
    AbstractSubstitutionCount::operator=(dsc);
    AbstractWeightedSubstitutionCount::operator=(dsc);
    AbstractSubstitutionDistance::operator=(dsc);
    DecompositionMethods::operator=(dsc);
    counts_         = dsc.counts_;
    currentLength_  = dsc.currentLength_;
    return *this;
  }

  virtual ~DecompositionSubstitutionCount() {}

  DecompositionSubstitutionCount* clone() const override { return new DecompositionSubstitutionCount(*this); }

public:
  double getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type = 1) const override;

  std::unique_ptr< Matrix<double> > getAllNumbersOfSubstitutions(double length, size_t type = 1) const override;

  void storeAllNumbersOfSubstitutions(double length, size_t type, Eigen::MatrixXd& mat) const override;

  std::vector<double> getNumberOfSubstitutionsPerType(size_t initialState, size_t finalState, double length) const override;

  /**
   * @brief Set the substitution model.
   *
   * @param model A pointer toward the substitution model to use.
   */
  void setSubstitutionModel(std::shared_ptr<const SubstitutionModelInterface> model) override;

protected:
  void initCounts_();

  void computeCounts_(double length) const;

  void substitutionRegisterHasChanged() override;

  void weightsHaveChanged() override;

  void distancesHaveChanged() override;

private:
  void fillBMatrices_();

  void setDistanceBMatrices_();
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_DECOMPOSITIONSUBSTITUTIONCOUNT_H
