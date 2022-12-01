//
// File: UniformizationSubstitutionCount.h
// Authors:
//   Julien Dutheil
// Created: 2011-03-19 13:54:00
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

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
      std::shared_ptr<const StateMapInterface> statemap,
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
