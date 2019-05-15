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

#include "DataFlowNumeric.h"

namespace bpp {
  std::string to_string (const NoDimension &) {
    return "()";
  }

  std::string to_string (const MatrixDimension & dim) {
    return "(" + std::to_string (dim.rows) + "," + std::to_string (dim.cols) + ")";
  }
  std::size_t hash (const MatrixDimension & dim) {
    std::size_t seed = 0;
    combineHash (seed, dim.rows);
    combineHash (seed, dim.cols);
    return seed;
  }

  namespace numeric {
    void checkDimensionIsSquare (const MatrixDimension & dim) {
      if (dim.rows != dim.cols) {
        throw Exception ("MatrixDimension is not square: " + std::to_string (dim.rows) + "x" +
                         std::to_string (dim.cols));
      }
    }
  } // namespace numeric

  namespace dataflow {
    void checkRecreateWithoutDependencies (const std::type_info & contextNodeType, const NodeRefVec & deps) {
      if (!deps.empty ()) {
        throw Exception (prettyTypeName (contextNodeType) +
                         "recreate called with dependencies, but node does not have dependencies");
      }
    }

    // Precompiled instantiations of numeric nodes
    template class ConstantZero<double>;
    template class ConstantZero<Parameter>;
    template class ConstantZero<Eigen::VectorXd>;
    template class ConstantZero<Eigen::RowVectorXd>;
    template class ConstantZero<Eigen::MatrixXd>;

    template class ConstantOne<double>;
    template class ConstantOne<Parameter>;
    template class ConstantOne<Eigen::VectorXd>;
    template class ConstantOne<Eigen::RowVectorXd>;
    template class ConstantOne<Eigen::MatrixXd>;

    template class NumericConstant<double>;
    template class NumericConstant<Eigen::VectorXd>;
    template class NumericConstant<Eigen::RowVectorXd>;
    template class NumericConstant<Eigen::MatrixXd>;

    template class NumericMutable<double>;
    template class NumericMutable<Eigen::VectorXd>;
    template class NumericMutable<Eigen::RowVectorXd>;
    template class NumericMutable<Eigen::MatrixXd>;

    template class Convert<double, double>;
    template class Convert<Eigen::VectorXd, Eigen::VectorXd>;
    template class Convert<Eigen::RowVectorXd, Eigen::RowVectorXd>;
    template class Convert<Eigen::RowVectorXd, Eigen::RowVectorXi>;
    template class Convert<Eigen::MatrixXd, Eigen::MatrixXd>;
    template class Convert<Eigen::VectorXd, double>;
    template class Convert<Eigen::RowVectorXd, double>;
    template class Convert<Eigen::MatrixXd, double>;
    template class Convert<Eigen::MatrixXd, Transposed<Eigen::MatrixXd>>;
    template class Convert<Eigen::RowVectorXd, Transposed<Eigen::VectorXd>>;
    template class Convert<Eigen::VectorXd, Transposed<Eigen::RowVectorXd>>;
  } // namespace dataflow
} // namespace bpp
