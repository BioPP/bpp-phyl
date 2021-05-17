//
// File: DataFlowNumeric.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2018-06-07
// Last modified: 2018-07-11
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#include <Bpp/Exceptions.h>

#include "DataFlowCWiseComputing.h"

namespace bpp {
  void failureDeltaNotDerivable (const std::type_info & contextNodeType) {
    throw Exception (prettyTypeName (contextNodeType) +
                     ": does not support derivation for the delta dependency");
  }
  void failureNumericalDerivationNotConfigured () {
    throw Exception ("Numerical derivation of expression is not configured: define the node "
                     "providing the delta value, and choose a computation type.");
  }

  // Precompiled instantiations of numeric nodes

  template class CWiseApply<MatrixLik, MatrixLik, TransitionFunction>;
    
  template class CWiseAdd<double, std::tuple<double, double>>;
  template class CWiseAdd<ExtendedFloat, std::tuple<ExtendedFloat, ExtendedFloat>>;
  template class CWiseAdd<VectorLik, std::tuple<VectorLik, VectorLik>>;
  template class CWiseAdd<RowLik, std::tuple<RowLik, RowLik>>;
  template class CWiseAdd<MatrixLik, std::tuple<MatrixLik, MatrixLik>>;
  template class CWiseAdd<TransitionFunction, std::tuple<TransitionFunction, TransitionFunction>>;

  template class CWiseAdd<RowLik, MatrixLik>;
  template class CWiseAdd<VectorLik, MatrixLik>;
  template class CWiseAdd<DataLik, VectorLik>;
  template class CWiseAdd<DataLik, RowLik>;

  template class CWiseAdd<double, ReductionOf<double>>;
  template class CWiseAdd<ExtendedFloat, ReductionOf<ExtendedFloat>>;
  template class CWiseAdd<VectorLik, ReductionOf<VectorLik>>;
  template class CWiseAdd<RowLik, ReductionOf<RowLik>>;
  template class CWiseAdd<MatrixLik, ReductionOf<MatrixLik>>;
  template class CWiseAdd<TransitionFunction, ReductionOf<TransitionFunction>>;

  template class CWiseMean<VectorLik, ReductionOf<VectorLik>, ReductionOf<double>>;
  template class CWiseMean<RowLik, ReductionOf<RowLik>, ReductionOf<double>>;
  template class CWiseMean<MatrixLik, ReductionOf<MatrixLik>, ReductionOf<double>>;
  template class CWiseMean<double, ReductionOf<double>, ReductionOf<double>>;
  template class CWiseMean<ExtendedFloat, ReductionOf<ExtendedFloat>, ReductionOf<double>>;

  template class CWiseMean<VectorLik, ReductionOf<VectorLik>, Eigen::VectorXd>;
  template class CWiseMean<RowLik, ReductionOf<RowLik>,  Eigen::VectorXd>;
  template class CWiseMean<MatrixLik, ReductionOf<MatrixLik>, Eigen::VectorXd>;
  template class CWiseMean<VectorLik, ReductionOf<VectorLik>, Eigen::RowVectorXd>;
  template class CWiseMean<RowLik, ReductionOf<RowLik>, Eigen::RowVectorXd>;
  template class CWiseMean<MatrixLik, ReductionOf<MatrixLik>, Eigen::RowVectorXd>;

  template class CWiseSub<double, std::tuple<double, double>>;
  template class CWiseSub<ExtendedFloat, std::tuple<ExtendedFloat, ExtendedFloat>>;
  template class CWiseSub<VectorLik, std::tuple<VectorLik, VectorLik>>;
  template class CWiseSub<RowLik, std::tuple<RowLik, RowLik>>;
  template class CWiseSub<MatrixLik, std::tuple<MatrixLik, MatrixLik>>;

  template class CWiseSub<VectorLik, std::tuple<VectorLik, DataLik>>;
  template class CWiseSub<RowLik, std::tuple<RowLik, DataLik>>;

  template class CWiseMul<double, std::tuple<double, double>>;
  template class CWiseMul<double, std::tuple<double, uint>>;
  template class CWiseMul<ExtendedFloat, std::tuple<ExtendedFloat, ExtendedFloat>>;
  template class CWiseMul<ExtendedFloat, std::tuple<ExtendedFloat, uint>>;
  template class CWiseMul<VectorLik, std::tuple<VectorLik, VectorLik>>;
  template class CWiseMul<RowLik, std::tuple<RowLik, RowLik>>;
  template class CWiseMul<MatrixLik, std::tuple<MatrixLik, MatrixLik>>;
  
  template class CWiseMul<RowLik, std::tuple<RowLik, Eigen::RowVectorXi>>;
  template class CWiseMul<VectorLik, std::tuple<VectorLik, Eigen::RowVectorXi>>;
  template class CWiseMul<VectorLik, std::tuple<DataLik, VectorLik>>;
  template class CWiseMul<RowLik, std::tuple<DataLik, RowLik>>;
  template class CWiseMul<MatrixLik, std::tuple<DataLik, MatrixLik>>;
  template class CWiseMul<TransitionFunction, std::tuple<TransitionFunction, TransitionFunction>>;
  template class CWiseMul<TransitionFunction, std::tuple<double, TransitionFunction>>;
    
  template class CWiseMul<double, ReductionOf<double>>;
  template class CWiseMul<ExtendedFloat, ReductionOf<ExtendedFloat>>;
  template class CWiseMul<VectorLik, ReductionOf<VectorLik>>;
  template class CWiseMul<RowLik, ReductionOf<RowLik>>;
  template class CWiseMul<MatrixLik, ReductionOf<MatrixLik>>;

  template class CWiseNegate<double>;
  template class CWiseNegate<ExtendedFloat>;
  template class CWiseNegate<VectorLik>;
  template class CWiseNegate<RowLik>;
  template class CWiseNegate<MatrixLik>;

  template class CWiseInverse<double>;
  template class CWiseInverse<ExtendedFloat>;
  template class CWiseInverse<VectorLik>;
  template class CWiseInverse<RowLik>;
  template class CWiseInverse<MatrixLik>;

  template class CWiseLog<double>;
  template class CWiseLog<ExtendedFloat>;
  template class CWiseLog<VectorLik>;
  template class CWiseLog<RowLik>;
  template class CWiseLog<MatrixLik>;
    
  template class CWiseExp<double>;
  template class CWiseExp<ExtendedFloat>;
  template class CWiseExp<VectorLik>;
  template class CWiseExp<RowLik>;
  template class CWiseExp<MatrixLik>;
    
  template class CWiseConstantPow<double>;
  template class CWiseConstantPow<ExtendedFloat>;
  template class CWiseConstantPow<VectorLik>;
  template class CWiseConstantPow<RowLik>;
  template class CWiseConstantPow<MatrixLik>;

  template class ScalarProduct<DataLik, VectorLik, VectorLik>;
  template class ScalarProduct<DataLik, RowLik, RowLik>;

  template class LogSumExp<DataLik, VectorLik, Eigen::VectorXd>;
  template class LogSumExp<DataLik, RowLik, Eigen::RowVectorXd>;

  template class SumOfLogarithms<VectorLik>;
  template class SumOfLogarithms<RowLik>;

  template class MatrixProduct<ExtendedFloatRowVectorXd, Eigen::RowVectorXd, ExtendedFloatMatrixXd>;
  template class MatrixProduct<Eigen::RowVectorXd, Eigen::RowVectorXd, Eigen::MatrixXd>;
  template class MatrixProduct<MatrixLik, Eigen::MatrixXd, MatrixLik>;
  template class MatrixProduct<MatrixLik, Transposed<Eigen::MatrixXd>, MatrixLik>;

  template class ShiftDelta<double>;
  template class ShiftDelta<VectorLik>;
  template class ShiftDelta<RowLik>;
  template class ShiftDelta<MatrixLik>;

  template class CombineDeltaShifted<double>;
  template class CombineDeltaShifted<VectorLik>;
  template class CombineDeltaShifted<RowLik>;
  template class CombineDeltaShifted<MatrixLik>;
  template class CombineDeltaShifted<TransitionFunction>;

} // namespace bpp
