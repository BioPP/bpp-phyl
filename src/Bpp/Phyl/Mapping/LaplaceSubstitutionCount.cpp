//
// File: LaplaceSubstitutionCount.cpp
// Created by: Julien Dutheil
// Created on: Wed Apr 5 11:21 2006
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#include "LaplaceSubstitutionCount.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace bpp;

/******************************************************************************/

void LaplaceSubstitutionCount::computeCounts(double length) const
{
  RowMatrix<double> Q = model_->getGenerator();
  // L is the diagonal matrix with all substitution rates.
  size_t s = Q.getNumberOfRows();
  RowMatrix<double> QL(s, s);
  for (size_t i = 0; i < s; i++)
  {
    for (size_t j = 0; j < s; j++)
    {
      QL(i, j) = ((i == j) ? 0. : Q(i, j));
    }
  }

  MatrixTools::fill(m_, 0.);
  RowMatrix<double> M2(s, s);
  RowMatrix<double> M3(s, s);
  RowMatrix<double> M4(s, s);
  RowMatrix<double> M5(s, s);
  for (size_t n = 1; n < cutOff_; ++n)
  {
    MatrixTools::fill(M2, 0.);
    for (size_t p = 0; p < n; ++p)
    {
      MatrixTools::pow(Q, p, M3);         // Q^p -> M5
      MatrixTools::mult(M3, QL, M4);      // Q^p . QL -> M4
      MatrixTools::pow(Q, n - p - 1, M3); // Q^(n-p-1) -> M3
      MatrixTools::mult(M4, M3, M5);      // Q^p . QL . Q^(n-p-1) -> M5
      MatrixTools::add(M2, M5);
    }
    MatrixTools::scale(M2, pow(length, static_cast<double>(n)) / static_cast<double>(NumTools::fact(n)));
    MatrixTools::add(m_, M2);
  }

  // Now we must divide by pijt:
  RowMatrix<double> P = model_->getPij_t(length);
  for (size_t i = 0; i < s; i++)
  {
    for (size_t j = 0; j < s; j++)
    {
      m_(i, j) /= P(i, j);
    }
  }
}

/******************************************************************************/

double LaplaceSubstitutionCount::getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type) const
{
  if (length == currentLength_)
    return m_(initialState, finalState);
  if (length < 0.000001)
    return initialState == finalState ? 0. : 1.;  // Limit case!
  // Else we need to recompute M:
  computeCounts(length);

  currentLength_ = length;
  return m_(initialState, finalState);
}

/******************************************************************************/

Matrix<double>* LaplaceSubstitutionCount::getAllNumbersOfSubstitutions(double length, size_t type) const
{
  if (length == currentLength_)
    return new RowMatrix<double>(m_);
  if (length < 0.000001) // Limit case!
  {
    size_t s = model_->getAlphabet()->getSize();
    for (size_t i = 0; i < s; i++)
    {
      for (size_t j = 0; j < s; j++)
      {
        m_(i, j) = i == j ? 0. : 1.;
      }
    }
  }
  else
  {
    // Else we need to recompute M:
    computeCounts(length);
  }

  currentLength_ = length;

  return new RowMatrix<double>(m_);
}

/******************************************************************************/

void LaplaceSubstitutionCount::setSubstitutionModel(const SubstitutionModel* model)
{
  model_ = model;
  size_t n = model->getAlphabet()->getSize();
  m_.resize(n, n);
  // Recompute counts:
  computeCounts(currentLength_);
}

/******************************************************************************/

