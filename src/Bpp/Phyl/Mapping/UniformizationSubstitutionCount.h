// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_UNIFORMIZATIONSUBSTITUTIONCOUNT_H
#define BPP_PHYL_MAPPING_UNIFORMIZATIONSUBSTITUTIONCOUNT_H

#include <Bpp/Numeric/Matrix/Matrix.h>

#include "SubstitutionDistance.h"
#include "WeightedSubstitutionCount.h"

namespace bpp
{
/**
 * @brief Analytical (weighted) substitution count using the uniformization method.
 *
 * The code is adapted from the original R code by Paula Tataru and Asger Hobolth.
 *
 * @author Julien Dutheil
 */
class UniformizationSubstitutionCount :
  public AbstractSubstitutionCount,
  public AbstractWeightedSubstitutionCount,
  public AbstractSubstitutionDistance
{
private:
  std::shared_ptr<const SubstitutionModelInterface> model_;
  size_t nbStates_;
  std::vector< RowMatrix<double> > bMatrices_;
  mutable std::vector< RowMatrix<double> > power_;
  mutable std::vector< std::vector< RowMatrix<double> > > s_;
  double miu_;
  mutable std::vector< RowMatrix<double> > counts_;
  mutable double currentLength_;

public:
  UniformizationSubstitutionCount(
      std::shared_ptr<const SubstitutionModelInterface> model,
      std::shared_ptr<const SubstitutionRegisterInterface> reg,
      std::shared_ptr<const AlphabetIndex2> weights = nullptr,
      std::shared_ptr<const AlphabetIndex2> distances = nullptr);

  UniformizationSubstitutionCount(
      const StateMapInterface& stateMap,
      std::shared_ptr<const SubstitutionRegisterInterface> reg,
      std::shared_ptr<const AlphabetIndex2> weights = nullptr,
      std::shared_ptr<const AlphabetIndex2> distances = nullptr);

  UniformizationSubstitutionCount(const UniformizationSubstitutionCount& usc) :
    AbstractSubstitutionCount(usc),
    AbstractWeightedSubstitutionCount(usc),
    AbstractSubstitutionDistance(usc),
    model_(usc.model_),
    nbStates_(usc.nbStates_),
    bMatrices_(usc.bMatrices_),
    power_(usc.power_),
    s_(usc.s_),
    miu_(usc.miu_),
    counts_(usc.counts_),
    currentLength_(usc.currentLength_)
  {}

  UniformizationSubstitutionCount& operator=(const UniformizationSubstitutionCount& usc)
  {
    AbstractSubstitutionCount::operator=(usc);
    AbstractWeightedSubstitutionCount::operator=(usc);
    AbstractSubstitutionDistance::operator=(usc);
    model_          = usc.model_;
    nbStates_       = usc.nbStates_;
    bMatrices_      = usc.bMatrices_;
    power_          = usc.power_;
    s_              = usc.s_;
    miu_            = usc.miu_;
    counts_         = usc.counts_;
    currentLength_  = usc.currentLength_;
    return *this;
  }

  virtual ~UniformizationSubstitutionCount() {}

  UniformizationSubstitutionCount* clone() const override { return new UniformizationSubstitutionCount(*this); }

public:
  double getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type = 1) const override;

  std::unique_ptr< Matrix<double> > getAllNumbersOfSubstitutions(double length, size_t type = 1) const override;

  void storeAllNumbersOfSubstitutions(double length, size_t type, Eigen::MatrixXd& mat) const override;

  std::vector<double> getNumberOfSubstitutionsPerType(size_t initialState, size_t finalState, double length) const override;

  void setSubstitutionModel(std::shared_ptr<const SubstitutionModelInterface> model) override;

protected:
  void computeCounts_(double length) const;
  void substitutionRegisterHasChanged() override;
  void weightsHaveChanged() override;
  void distancesHaveChanged() override;

private:
  void resetBMatrices_();
  void initBMatrices_();
  void fillBMatrices_();

  void setDistanceBMatrices_();
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_UNIFORMIZATIONSUBSTITUTIONCOUNT_H
