// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Exceptions.h>

#include "DataFlowNumeric.h"

using namespace Eigen;

namespace bpp
{
/*****************************************/
/* copyBppToEigen */


template void copyBppToEigen (const std::vector<ExtendedFloat>& bppVector, Eigen::RowVectorXd& eigenVector);
template void copyBppToEigen (const std::vector<ExtendedFloat>& bppVector, Eigen::VectorXd& eigenVector);
template void copyBppToEigen (const std::vector<ExtendedFloat>& bppVector, ExtendedFloatRowVectorXd& eigenVector);
template void copyBppToEigen (const std::vector<ExtendedFloat>& bppVector, ExtendedFloatVectorXd& eigenVector);

template void copyBppToEigen (const bpp::Vdouble& bppVector, Eigen::RowVectorXd& eigenVector);
template void copyBppToEigen (const bpp::Vdouble& bppVector, Eigen::VectorXd& eigenVector);
template void copyBppToEigen (const bpp::Vdouble& bppVector, ExtendedFloatRowVectorXd& eigenVector);
template void copyBppToEigen (const bpp::Vdouble& bppVector, ExtendedFloatVectorXd& eigenVector);


// template void copyBppToEigen(bpp::Matrix<double> const&, Eigen::Matrix<double, -1, -1, 0, -1, -1>&);


void copyBppToEigen (const std::vector<ExtendedFloatVectorXd>& bppVector, ExtendedFloatMatrixXd& eigenVector)
{
  // Look for largest extendedfloat
  ExtendedFloat::ExtType maxE = std::max_element(bppVector.begin(), bppVector.end(), [](const ExtendedFloatVectorXd& lhs, const ExtendedFloatVectorXd& rhs){
      return lhs.exponent_part() < rhs.exponent_part();
    })->exponent_part();

  eigenVector.exponent_part() = maxE;
  eigenVector.float_part() = Eigen::MatrixXd::NullaryExpr(static_cast<Eigen::Index>(bppVector[0].rows()),
        static_cast<Eigen::Index>(bppVector.size()),
        [&maxE, &bppVector](int i, int j){
      return bppVector[(size_t)j].float_part()(i) * bpp::constexpr_power<double> (bpp::ExtendedFloat::radix, bppVector[(size_t)j].exponent_part () - maxE);
    });
}


void copyBppToEigen (const std::vector<Eigen::VectorXd>& bppVector, Eigen::MatrixXd& eigenMatrix)
{
  // Look for largest extendedfloat
  const auto eigenRows = static_cast<Eigen::Index>(bppVector[0].rows());
  const auto eigenCols = static_cast<Eigen::Index>(bppVector.size());

  eigenMatrix.resize (eigenRows, eigenCols);
  eigenMatrix.fill(0);
  for (Eigen::Index j = 0; j < eigenCols; ++j)
  {
    for (Eigen::Index i = 0; i < eigenRows; ++i)
    {
      eigenMatrix (i, j) = bppVector[(std::size_t)j](i);
    }
  }
}


void copyBppToEigen (const bpp::Matrix<double>& bppMatrix, ExtendedFloatMatrixXd& eigenMatrix)
{
  const auto eigenRows = static_cast<Eigen::Index>(bppMatrix.getNumberOfRows ());
  const auto eigenCols = static_cast<Eigen::Index>(bppMatrix.getNumberOfColumns ());
  eigenMatrix.resize (eigenRows, eigenCols);
  eigenMatrix.fill(0);
  for (Eigen::Index i = 0; i < eigenRows; ++i)
  {
    for (Eigen::Index j = 0; j < eigenCols; ++j)
    {
      eigenMatrix.float_part() (i, j) = bppMatrix (static_cast<std::size_t>(i), static_cast<std::size_t>(j));
    }
  }
#ifdef DEBUG
  std::cerr << "=== copyBppToEigen(" << typeid(bppMatrix).name() << ", " << typeid(eigenMatrix).name() << ") ===" << std::endl;
  std::cerr << &bppMatrix << std::endl;
  std::cerr << eigenRows << "," << eigenCols << std::endl;
  std::cerr << eigenMatrix << std::endl;
  std::cerr << "=== end copyBppToEigen === " << std::endl;
#endif
}

void copyBppToEigen (const bpp::Matrix<double>& bppMatrix, Eigen::MatrixXd& eigenMatrix)
{
  const auto eigenRows = static_cast<Eigen::Index>(bppMatrix.getNumberOfRows ());
  const auto eigenCols = static_cast<Eigen::Index>(bppMatrix.getNumberOfColumns ());
  eigenMatrix.resize (eigenRows, eigenCols);
  eigenMatrix.fill(0);
  for (Eigen::Index i = 0; i < eigenRows; ++i)
  {
    for (Eigen::Index j = 0; j < eigenCols; ++j)
    {
      eigenMatrix (i, j) = bppMatrix (static_cast<std::size_t>(i), static_cast<std::size_t>(j));
    }
  }
#ifdef DEBUG
  std::cerr << "copyBppToEigen(" << typeid(bppMatrix).name() << ", " << typeid(eigenMatrix).name() << ")" << std::endl;
  std::cerr << &bppMatrix << std::endl;
  std::cerr << eigenRows << "," << eigenCols << std::endl;
  std::cerr << eigenMatrix << std::endl;
#endif
}


/*****************************************/
/* copyEigenToBpp */

template void copyEigenToBpp(const ExtendedFloatMatrixXd& eigenMatrix, std::vector<std::vector<double>>& bppMatrix);
template void copyEigenToBpp(const ExtendedFloatMatrixXd& eigenMatrix, std::vector<std::vector<bpp::ExtendedFloat>>& bppMatrix);

template void copyEigenToBpp(const Eigen::MatrixXd& eigenMatrix, std::vector<std::vector<double>>& bppMatrix);
template void copyEigenToBpp(const Eigen::MatrixXd& eigenMatrix, std::vector<std::vector<bpp::ExtendedFloat>>& bppMatrix);

void copyEigenToBpp (const MatrixLik& eigenMatrix, bpp::Matrix<double>& bppMatrix)
{
  const auto eigenRows = static_cast<std::size_t>(eigenMatrix.rows());
  const auto eigenCols = static_cast<std::size_t>(eigenMatrix.cols());

  bppMatrix.resize (eigenRows, eigenCols);
  for (size_t i = 0; i < eigenRows; ++i)
  {
    for (size_t j = 0; j < eigenCols; ++j)
    {
      bppMatrix(i, j) = convert (eigenMatrix (static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)));
    }
  }
}

template void copyEigenToBpp (const RowLik& eigenVector, Vdouble& bppVector);
template void copyEigenToBpp (const VectorLik& eigenVector, Vdouble& bppVector);
template void copyEigenToBpp (const RowLik& eigenVector, std::vector<ExtendedFloat>& bppVector);
template void copyEigenToBpp (const VectorLik& eigenVector, std::vector<ExtendedFloat>& bppVector);


/***************************************************/
/*  to string  */


std::string to_string (const NoDimension&)
{
  return "()";
}

std::string to_string (const MatrixDimension& dim)
{
  return "(" + std::to_string (dim.rows) + "," + std::to_string (dim.cols) + ")";
}

std::size_t hash (const MatrixDimension& dim)
{
  std::size_t seed = 0;
  combineHash (seed, dim.rows);
  combineHash (seed, dim.cols);
  return seed;
}

namespace numeric
{
void checkDimensionIsSquare (const MatrixDimension& dim)
{
  if (dim.rows != dim.cols)
  {
    throw Exception ("MatrixDimension is not square: " + std::to_string (dim.rows) + "x" +
          std::to_string (dim.cols));
  }
}
} // namespace numeric

void checkRecreateWithoutDependencies (const std::type_info& contextNodeType, const NodeRefVec& deps)
{
  if (!deps.empty ())
  {
    throw Exception (prettyTypeName (contextNodeType) +
          "recreate called with dependencies, but node does not have dependencies");
  }
}

// Precompiled instantiations of numeric nodes
template class ConstantZero<unsigned int>;
template class ConstantZero<double>;
template class ConstantZero<char>;
template class ConstantZero<std::string>;

template class ConstantZero<Parameter>;
template class ConstantZero<TransitionFunction>;

template class ConstantZero<Eigen::VectorXd>;
template class ConstantZero<Eigen::RowVectorXd>;
template class ConstantZero<Eigen::MatrixXd>;

template class ConstantZero<ExtendedFloatVectorXd>;
template class ConstantZero<ExtendedFloatRowVectorXd>;
template class ConstantZero<ExtendedFloatMatrixXd>;


template class ConstantOne<unsigned int>;
template class ConstantOne<double>;
template class ConstantOne<char>;
template class ConstantOne<std::string>;

template class ConstantOne<Parameter>;
template class ConstantOne<TransitionFunction>;

template class ConstantOne<Eigen::VectorXd>;
template class ConstantOne<Eigen::RowVectorXd>;
template class ConstantOne<Eigen::MatrixXd>;

template class ConstantOne<ExtendedFloatVectorXd>;
template class ConstantOne<ExtendedFloatRowVectorXd>;
template class ConstantOne<ExtendedFloatMatrixXd>;


template class Identity<double>;
template class Identity<ExtendedFloatMatrixXd>;
template class Identity<Eigen::MatrixXd>;


template class NumericConstant<unsigned int>;
template class NumericConstant<double>;
template class NumericConstant<size_t>;
template class NumericConstant<std::string>;

template class NumericConstant<ExtendedFloatVectorXd>;
template class NumericConstant<ExtendedFloatRowVectorXd>;
template class NumericConstant<ExtendedFloatMatrixXd>;

NumericConstant<char> NodeX('X');

template class NumericConstant<Eigen::VectorXd>;
template class NumericConstant<Eigen::RowVectorXd>;
template class NumericConstant<Eigen::MatrixXd>;

template class NumericMutable<unsigned int>;
template class NumericMutable<double>;

template class NumericMutable<ExtendedFloatVectorXd>;
template class NumericMutable<ExtendedFloatRowVectorXd>;
template class NumericMutable<ExtendedFloatMatrixXd>;

template class NumericMutable<Eigen::VectorXd>;
template class NumericMutable<Eigen::RowVectorXd>;
template class NumericMutable<Eigen::MatrixXd>;

template class Convert<double, double>;

template class Convert<ExtendedFloatVectorXd, double>;
template class Convert<ExtendedFloatVectorXd, ExtendedFloatVectorXd>;
template class Convert<ExtendedFloatVectorXd, Eigen::VectorXd>;
template class Convert<ExtendedFloatVectorXd, Eigen::VectorXi>;
template class Convert<ExtendedFloatVectorXd, Transposed<ExtendedFloatRowVectorXd>>;
template class Convert<ExtendedFloatVectorXd, Transposed<Eigen::RowVectorXd>>;

template class Convert<Eigen::VectorXd, ExtendedFloatVectorXd>;

template class Convert<ExtendedFloatRowVectorXd, double>;
template class Convert<ExtendedFloatRowVectorXd, ExtendedFloatRowVectorXd>;
template class Convert<ExtendedFloatRowVectorXd, Eigen::RowVectorXd>;
template class Convert<ExtendedFloatRowVectorXd, Eigen::RowVectorXi>;
template class Convert<ExtendedFloatRowVectorXd, Transposed<ExtendedFloatVectorXd>>;
template class Convert<ExtendedFloatRowVectorXd, Transposed<Eigen::VectorXd>>;

template class Convert<ExtendedFloatMatrixXd, double>;
template class Convert<ExtendedFloatMatrixXd, ExtendedFloatMatrixXd>;
template class Convert<ExtendedFloatMatrixXd, Eigen::MatrixXd>;
template class Convert<ExtendedFloatMatrixXd, Transposed<ExtendedFloatMatrixXd>>;
template class Convert<ExtendedFloatMatrixXd, Transposed<Eigen::MatrixXd>>;

template class Convert<Eigen::VectorXd, double>;
template class Convert<Eigen::VectorXd, Eigen::VectorXd>;
template class Convert<Eigen::VectorXd, Eigen::VectorXi>;
template class Convert<Eigen::VectorXd, Transposed<Eigen::RowVectorXd>>;

template class Convert<Eigen::RowVectorXd, double>;
template class Convert<Eigen::RowVectorXd, Eigen::RowVectorXd>;
template class Convert<Eigen::RowVectorXd, Eigen::RowVectorXi>;
template class Convert<Eigen::RowVectorXd, Transposed<Eigen::VectorXd>>;

template class Convert<Eigen::MatrixXd, double>;
template class Convert<Eigen::MatrixXd, Eigen::MatrixXd>;
template class Convert<Eigen::MatrixXd, Transposed<Eigen::MatrixXd>>;
} // namespace bpp
