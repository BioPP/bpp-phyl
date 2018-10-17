//
// File: BinarySubstitutionModel.h
// Created by: Keren Halabi
// Created on: 2018
//

/*
   Copyright or � or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _TWOPARAMETERBINARYSUBSTITUTIONMODEL_H_
#define _TWOPARAMETERBINARYSUBSTITUTIONMODEL_H_

#include "AbstractSubstitutionModel.h"
#include <Bpp/Seq/Alphabet/BinaryAlphabet.h>

namespace bpp
{
/**
 * @brief The Model on two states
 *
 * \f[
 * Q = mu.\begin{pmatrix}
 * -\pi_{0} & \pi_{0}  \\
 * \pi_{1} & -\pi_{1}  \\
 * \end{pmatrix}
 * \f]
 */

class TwoParameterBinarySubstitutionModel :
  public AbstractReversibleSubstitutionModel
{
private:
  double mu_;
  double pi0_;

public:
  TwoParameterBinarySubstitutionModel(const BinaryAlphabet* alpha, double mu = 1., double pi0 = 0.5);

  virtual ~TwoParameterBinarySubstitutionModel() {}

  TwoParameterBinarySubstitutionModel* clone() const { return new TwoParameterBinarySubstitutionModel(*this); }

  std::string getName() const { return "TwoParameterBinary"; }

  size_t getNumberOfStates() const { return 2; }

  void setMuBounds(double lb, double ub);

protected:
  void updateMatrices();

};
} // end of namespace bpp.

#endif  // _BINARYSUBSTITUTIONMODEL_H_

