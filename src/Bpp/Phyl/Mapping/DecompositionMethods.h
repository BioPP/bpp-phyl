// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_DECOMPOSITIONMETHODS_H
#define BPP_PHYL_MAPPING_DECOMPOSITIONMETHODS_H

#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "../Model/SubstitutionModel.h"
#include "SubstitutionRegister.h"

namespace bpp
{
/**
 * @brief Methods useful for analytical substitution count and rewards
 * using the eigen decomposition method.
 *
 * The code is adapted from the original R code by Paula Tataru and
 * Asger Hobolth.
 *
 * @author Julien Dutheil, Laurent Guéguen
 */

class DecompositionMethods
{
protected:
  std::shared_ptr<const SubstitutionModelInterface> model_;
  size_t nbStates_;
  size_t nbTypes_;
  mutable RowMatrix<double> jMat_, jIMat_;

  /**
   * @brief Real and imaginary eigenvectors, for non-reversible
   * computation
   */
  ColMatrix<double> rightEigenVectors_, rightIEigenVectors_;
  RowMatrix<double> leftEigenVectors_, leftIEigenVectors_;

  /**
   * @brief computation matrices
   */
  std::vector< RowMatrix<double> > bMatrices_, insideProducts_, insideIProducts_;

public:
  DecompositionMethods(
      std::shared_ptr<const SubstitutionModelInterface> model,
      std::shared_ptr<const SubstitutionRegisterInterface> reg);

  DecompositionMethods(std::shared_ptr<const SubstitutionRegisterInterface> reg);

  DecompositionMethods(std::shared_ptr<const SubstitutionModelInterface> model);

  DecompositionMethods(const StateMapInterface& stateMap);

  DecompositionMethods(const DecompositionMethods& dm) :
    model_(dm.model_),
    nbStates_(dm.nbStates_),
    nbTypes_(dm.nbTypes_),
    jMat_(dm.jMat_),
    jIMat_(dm.jIMat_),
    rightEigenVectors_(dm.rightEigenVectors_),
    rightIEigenVectors_(dm.rightIEigenVectors_),
    leftEigenVectors_(dm.leftEigenVectors_),
    leftIEigenVectors_(dm.leftIEigenVectors_),
    bMatrices_(dm.bMatrices_),
    insideProducts_(dm.insideProducts_),
    insideIProducts_(dm.insideIProducts_)
  {}

  DecompositionMethods& operator=(const DecompositionMethods& dm)
  {
    model_          = dm.model_;
    nbStates_       = dm.nbStates_;
    nbTypes_        = dm.nbTypes_;
    jMat_           = dm.jMat_;
    jIMat_          = dm.jIMat_;

    rightEigenVectors_ = dm.rightEigenVectors_;
    rightIEigenVectors_ = dm.rightIEigenVectors_;
    leftEigenVectors_ = dm.leftEigenVectors_;
    leftIEigenVectors_ = dm.leftIEigenVectors_;
    bMatrices_      = dm.bMatrices_;
    insideProducts_ = dm.insideProducts_;
    insideIProducts_ =  dm.insideIProducts_;

    return *this;
  }

  DecompositionMethods* clone() const { return new DecompositionMethods(*this); }

  virtual ~DecompositionMethods() {}


  /**
   * @brief Set the substitution model.
   *
   * @param model A pointer toward the substitution model to use.
   */
  void setSubstitutionModel(std::shared_ptr<const SubstitutionModelInterface> model);

protected:
  void initStates_();

  void initBMatrices_();

  void computeProducts_();

  /*
   * @brief Perform the computation of the conditional expectations
   *
   */

  void computeExpectations(RowMatrix<double>& mapping, double length) const;

  void computeExpectations(std::vector< RowMatrix<double> >& mappings, double length) const;

  /**
   * @brief Compute the integral part of the computation
   *
   */

  void jFunction_(const std::vector<double>& lambda, double t, RowMatrix<double>& result) const;

  /**
   * @brief Compute the integral part of the computation, in complex numbers
   *
   */

  void jFunction_(const std::vector<double>& lambda, const std::vector<double>& ilambda, double t, RowMatrix<double>& result, RowMatrix<double>& iresult) const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_DECOMPOSITIONMETHODS_H
