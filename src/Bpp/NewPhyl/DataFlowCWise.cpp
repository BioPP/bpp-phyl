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

#include "DataFlowCWise.h"

namespace bpp {
  namespace dataflow {
    void failureDeltaNotDerivable (const std::type_info & contextNodeType) {
      throw Exception (prettyTypeName (contextNodeType) +
                       ": does not support derivation for the delta dependency");
    }
    void failureNumericalDerivationNotConfigured () {
      throw Exception ("Numerical derivation of expression is not configured: define the node "
                       "providing the delta value, and choose a computation type.");
    }

    // Precompiled instantiations of numeric nodes
    template class CWiseFill<Eigen::RowVectorXd, double>;
    //template class CWiseFill<Eigen::VectorXd, double>;
    //template class CWiseFill<Eigen::MatrixXd, Eigen::VectorXd>;
    //template class CWiseFill<Eigen::MatrixXd, Eigen::RowVectorXd>;

    template class CWiseAdd<double, std::tuple<double, double>>;
    template class CWiseAdd<Eigen::VectorXd, std::tuple<Eigen::VectorXd, Eigen::VectorXd>>;
    template class CWiseAdd<Eigen::RowVectorXd, std::tuple<Eigen::RowVectorXd, Eigen::RowVectorXd>>;
    template class CWiseAdd<Eigen::MatrixXd, std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>>;

    template class CWiseAdd<double, ReductionOf<double>>;
    template class CWiseAdd<Eigen::VectorXd, ReductionOf<Eigen::VectorXd>>;
    template class CWiseAdd<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>>;
    template class CWiseAdd<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

    template class CWiseMean<Eigen::VectorXd, ReductionOf<Eigen::VectorXd>, Eigen::VectorXd>;
    template class CWiseMean<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>, Eigen::VectorXd>;
    template class CWiseMean<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>, Eigen::VectorXd>;
    template class CWiseMean<Eigen::VectorXd, ReductionOf<Eigen::VectorXd>, Eigen::RowVectorXd>;
    template class CWiseMean<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>, Eigen::RowVectorXd>;
    template class CWiseMean<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>, Eigen::RowVectorXd>;


    template class CWiseSub<double, std::tuple<double, double>>;
    template class CWiseSub<Eigen::VectorXd, std::tuple<Eigen::VectorXd, Eigen::VectorXd>>;
    template class CWiseSub<Eigen::VectorXd, std::tuple<double, Eigen::VectorXd>>;
    template class CWiseSub<Eigen::RowVectorXd, std::tuple<Eigen::RowVectorXd, Eigen::RowVectorXd>>;
    template class CWiseSub<Eigen::RowVectorXd, std::tuple<double, Eigen::RowVectorXd>>;
    template class CWiseSub<Eigen::MatrixXd, std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>>;

    template class CWiseMul<double, std::tuple<double, double>>;
    template class CWiseMul<Eigen::VectorXd, std::tuple<Eigen::VectorXd, Eigen::VectorXd>>;
    template class CWiseMul<Eigen::RowVectorXd, std::tuple<Eigen::RowVectorXd, Eigen::RowVectorXd>>;
    template class CWiseMul<Eigen::MatrixXd, std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>>;
    template class CWiseMul<Eigen::VectorXd, std::tuple<double, Eigen::VectorXd>>;
    template class CWiseMul<Eigen::RowVectorXd, std::tuple<double, Eigen::RowVectorXd>>;
    template class CWiseMul<Eigen::MatrixXd, std::tuple<double, Eigen::MatrixXd>>;

    template class CWiseMul<double, ReductionOf<double>>;
    template class CWiseMul<Eigen::VectorXd, ReductionOf<Eigen::VectorXd>>;
    template class CWiseMul<Eigen::RowVectorXd, ReductionOf<Eigen::RowVectorXd>>;
    template class CWiseMul<Eigen::MatrixXd, ReductionOf<Eigen::MatrixXd>>;

    template class CWiseNegate<double>;
    template class CWiseNegate<Eigen::VectorXd>;
    template class CWiseNegate<Eigen::RowVectorXd>;
    template class CWiseNegate<Eigen::MatrixXd>;

    template class CWiseInverse<double>;
    template class CWiseInverse<Eigen::VectorXd>;
    template class CWiseInverse<Eigen::RowVectorXd>;
    template class CWiseInverse<Eigen::MatrixXd>;

    template class CWiseLog<double>;
    template class CWiseLog<Eigen::VectorXd>;
    template class CWiseLog<Eigen::RowVectorXd>;
    template class CWiseLog<Eigen::MatrixXd>;
    
    template class CWiseExp<double>;
    template class CWiseExp<Eigen::VectorXd>;
    template class CWiseExp<Eigen::RowVectorXd>;
    template class CWiseExp<Eigen::MatrixXd>;
    
    template class CWiseConstantPow<double>;
    template class CWiseConstantPow<Eigen::VectorXd>;
    template class CWiseConstantPow<Eigen::RowVectorXd>;
    template class CWiseConstantPow<Eigen::MatrixXd>;

    template class ScalarProduct<Eigen::VectorXd, Eigen::VectorXd>;
    template class ScalarProduct<Eigen::RowVectorXd, Eigen::RowVectorXd>;

    template class LogSumExp<Eigen::VectorXd, Eigen::VectorXd>;
    template class LogSumExp<Eigen::RowVectorXd, Eigen::RowVectorXd>;

    template class SumOfLogarithms<Eigen::VectorXd>;
    template class SumOfLogarithms<Eigen::RowVectorXd>;

    template class MatrixProduct<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>;
    template class MatrixProduct<Eigen::RowVectorXd, Eigen::RowVectorXd, Eigen::MatrixXd>;
    template class MatrixProduct<Eigen::MatrixXd, Transposed<Eigen::MatrixXd>, Eigen::MatrixXd>;

    template class ShiftDelta<double>;
    template class ShiftDelta<Eigen::VectorXd>;
    template class ShiftDelta<Eigen::RowVectorXd>;
    template class ShiftDelta<Eigen::MatrixXd>;

    template class CombineDeltaShifted<double>;
    template class CombineDeltaShifted<Eigen::VectorXd>;
    template class CombineDeltaShifted<Eigen::RowVectorXd>;
    template class CombineDeltaShifted<Eigen::MatrixXd>;
  } // namespace dataflow
} // namespace bpp
