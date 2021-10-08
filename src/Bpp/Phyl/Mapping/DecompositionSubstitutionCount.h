//
// File: DecompositionSubstitutionCount.h
// Authors:
//   Julien Dutheil
// Created: 2011-03-17 16:08:00
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)
  
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
  DecompositionSubstitutionCount(const SubstitutionModel* model, SubstitutionRegister* reg, std::shared_ptr<const AlphabetIndex2> weights = 0, std::shared_ptr<const AlphabetIndex2> distances = 0);

  DecompositionSubstitutionCount(SubstitutionRegister* reg, std::shared_ptr<const AlphabetIndex2> weights = 0, std::shared_ptr<const AlphabetIndex2> distances = 0);

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

  DecompositionSubstitutionCount* clone() const { return new DecompositionSubstitutionCount(*this); }

public:
  double getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type = 1) const;

  Matrix<double>* getAllNumbersOfSubstitutions(double length, size_t type = 1) const;

  void storeAllNumbersOfSubstitutions(double length, size_t type, Eigen::MatrixXd& mat) const;

  std::vector<double> getNumberOfSubstitutionsPerType(size_t initialState, size_t finalState, double length) const;

  /**
   * @brief Set the substitution model.
   *
   * @param model A pointer toward the substitution model to use.
   */

  void setSubstitutionModel(const SubstitutionModel* model);

protected:
  void initCounts_();

  void computeCounts_(double length) const;

  void substitutionRegisterHasChanged();

  void weightsHaveChanged();

  void distancesHaveChanged();

private:
  void fillBMatrices_();

  void setDistanceBMatrices_();
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_DECOMPOSITIONSUBSTITUTIONCOUNT_H
