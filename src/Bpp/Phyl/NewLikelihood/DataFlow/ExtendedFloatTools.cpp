#include "ExtendedFloatTools.h"

using namespace bpp;
using namespace std;

void ExtendedFloatTools::scale(MatrixXef& A, ExtendedFloat a, ExtendedFloat b){
  size_t nrows = A.rows();
  size_t ncols = A.cols();
  for (size_t i = 0; i < nrows; i++){
    for (size_t j = 0; j < ncols; j++){
      A(i, j) = a * A(i, j) + b;
    }
  }
}
/******************************************************************************************/
void ExtendedFloatTools::scale(MatrixXef& A, double a, double b){
  ExtendedFloat a_Ef = ExtendedFloat{a};
  ExtendedFloat b_Ef = ExtendedFloat{b};
  a_Ef.normalize();
  b_Ef.normalize();
  scale(A, a_Ef, b_Ef);

}
/******************************************************************************************/
void ExtendedFloatTools::getIdentityMatrix(MatrixXef &outMat, size_t size){
    outMat.resize(size, size);
    for (size_t i = 0; i < size; i++){
        for (size_t j = 0; j < size; j++){
            if (i != j){
                outMat(i,j) = 0;
                outMat(i,j).normalize();
            }else{
                outMat(i,j) = 1;
                outMat(i,j).normalize();
            }
        }
    }
    return;
}
/******************************************************************************************/
MatrixXef ExtendedFloatTools::convertBppMatrixToEigenEF(const Matrix <double>& mat){
    size_t nCols = mat.getNumberOfColumns();
    size_t nRows = mat.getNumberOfRows();
    MatrixXef matEf(nRows, nCols);
    for (size_t i = 0; i < nRows; i++){
        for (size_t j = 0; j < nCols; j++){
            ExtendedFloat ef = mat(i, j);
            ef.normalize();
            matEf(i, j) = ef;
        }
    }
    return matEf;
}
/******************************************************************************************/
Eigen::MatrixXd ExtendedFloatTools::convertBppMatrixToEigenDB(const Matrix <double>& mat){
    size_t nCols = mat.getNumberOfColumns();
    size_t nRows = mat.getNumberOfRows();
    Eigen::MatrixXd matEigen(nRows, nCols);
    for (size_t i = 0; i < nRows; i++){
        for (size_t j = 0; j < nCols; j++){
            matEigen(i, j) = mat(i, j);
        }
    }
    return matEigen;
}
/******************************************************************************************/
Eigen::MatrixXd ExtendedFloatTools::convertMatToDouble (const MatrixXef & eigenMatrix){
    size_t nCols = eigenMatrix.cols();
    size_t nRows = eigenMatrix.rows();
    Eigen::MatrixXd dbMat(nRows, nCols);
    for (size_t i = 0; i < nRows; i ++){
        for (size_t j = 0; j < nCols; j++){
            ExtendedFloat val = eigenMatrix(i,j);          
            dbMat(i, j) = convert(val);
        }
    }
    return dbMat;
}
/******************************************************************************************/
void ExtendedFloatTools::convertMatToDouble (const MatrixXef & eigenMatrix, Eigen::MatrixXd & dbMat){
    size_t nCols = eigenMatrix.cols();
    size_t nRows = eigenMatrix.rows();
    dbMat.resize(nRows, nCols);
    for (size_t i = 0; i < nRows; i ++){
        for (size_t j = 0; j < nCols; j++){
            ExtendedFloat val = eigenMatrix(i,j);          
            dbMat(i, j) = convert(val);
        }
    }
}
/******************************************************************************************/
void ExtendedFloatTools::convertDoubleMatToExtendedFloat(const Eigen::MatrixXd & dbMat, MatrixXef & eigenMatrix){
    size_t nCols = dbMat.cols();
    size_t nRows = dbMat.rows();
    eigenMatrix.resize(nRows, nCols);
    for (size_t i = 0; i < nRows; i ++){
        for (size_t j = 0; j < nCols; j++){
            ExtendedFloat val = ExtendedFloat{dbMat(i,j)};
            val.normalize();
            eigenMatrix(i, j) = val;
           
        }
    }
}
/******************************************************************************************/
ExtendedFloat ExtendedFloatTools::exp(ExtendedFloat &ef){
    double exp = convert(ef);
    double e = std::exp(1);
    ExtendedFloat e_ef = ExtendedFloat{e};
    e_ef.normalize();
    ExtendedFloat res = e_ef.pow(exp);
    return res;

}
/******************************************************************************************/
void ExtendedFloatTools::copyToEigenVector(const std::vector <double> & stdVector, VectorXef & eigenEfVector){
    eigenEfVector.resize(stdVector.size());
    for (size_t i = 0; i < stdVector.size(); i++){
        eigenEfVector[i] = stdVector[i];
        eigenEfVector[i].normalize();
    }
}
/******************************************************************************************/
void ExtendedFloatTools::Taylor(const MatrixXef& A, size_t p, std::vector<MatrixXef> & vO){
    size_t n = A.rows();
    if (n != static_cast<size_t>(A.cols()))
        throw DimensionException("test_ExtendedFloat::pow(). nrows != ncols.", static_cast<size_t>(A.cols()), static_cast<size_t>(A.rows()));
    vO.resize(p+1);
    getIdentityMatrix(vO[0],n);
    vO[1] = A;   
    for (size_t i = 1; i < p; i++)
    {
        //mult(vO[i], A, vO[i+1]);
        vO[i+1] = vO[i] * A;
    }
}
/******************************************************************************************/
VectorXef ExtendedFloatTools::exp(VectorXef v_in){
    size_t nRows = v_in.rows();
    VectorXef v_out(nRows);
    for (size_t i = 0; i < nRows; i++){
        v_out[i] = exp(v_in[i]);
    }
    return v_out;
}
/******************************************************************************************/
void ExtendedFloatTools::mult(const MatrixXef& A, const VectorXef& D, const VectorXef& U, const VectorXef& L, const MatrixXef& B, MatrixXef& O)
{
    size_t ncA = A.cols();
    size_t nrA = A.rows();
    size_t nrB = B.rows();
    size_t ncB = B.cols();
    if (ncA != nrB) throw DimensionException("MatrixTools::mult(). nrows B != ncols A.", nrB, ncA);
    if (ncA != static_cast<size_t>(D.rows())) throw DimensionException("MatrixTools::mult(). Vector size is not equal to matrix size.", static_cast<size_t>(D.rows()), ncA);
    if (ncA != static_cast<size_t>(U.rows())+1) throw DimensionException("MatrixTools::mult(). Vector size is not equal to matrix size-1.", static_cast<size_t>(U.rows()), ncA);
    if (ncA != static_cast<size_t>(L.rows())+1) throw DimensionException("MatrixTools::mult(). Vector size is not equal to matrix size-1.", static_cast<size_t>(L.rows()), ncA);
    O.resize(nrA, ncB);
    for (size_t i = 0; i < nrA; i++)
    {
        for (size_t j = 0; j < ncB; j++)
        {
          O(i, j) = A(i, 0) * D[0] * B(0, j);
          if (nrB>1)
            O(i, j) += A(i,0) * U[0] * B(1,j);
          for (size_t k = 1; k < ncA-1; k++)
          {
            O(i, j) += A(i, k) * (L[k-1] * B(k-1, j) + D[k] * B(k, j) + U[k] * B(k+1,j));
          }
          if (ncA>=2)
            O(i, j) += A(i, ncA-1) * L[ncA-2] * B(ncA-2, j);
          O(i,j) += A(i, ncA-1) * D[ncA-1] * B(ncA-1, j);
        }
    }
}
/******************************************************************************************/
void ExtendedFloatTools::add(MatrixXef& A, ExtendedFloat& x, const MatrixXef& B)
{
    size_t ncA = A.cols();
    size_t nrA = A.rows();
    size_t nrB = B.rows();
    size_t ncB = B.cols();
    if (ncA != ncB) throw DimensionException("test_ExtendedFloat::operator+(). A and B must have the same number of columns.", ncB, ncA);
    if (nrA != nrB) throw DimensionException("test_ExtendedFloat::operator+(). A and B must have the same number of rows.", nrB, nrA);
      
    for (size_t i = 0; i < nrA; i++){
        for (size_t j = 0; j < ncA; j++){
          A(i, j) += x*B(i, j);
        }
    }
}
/******************************************************************************************/
void ExtendedFloatTools::add(MatrixXef& A, double& x, const MatrixXef& B)
{
    ExtendedFloat x_ef = convertToExtendedFloat(x);
    add(A, x_ef, B);

}
/******************************************************************************************/
ExtendedFloat ExtendedFloatTools::abs(ExtendedFloat ef){
  ExtendedFloat zero = ExtendedFloat{0};
  ExtendedFloat minusOne = ExtendedFloat{-1};
  zero.normalize();
  minusOne.normalize();
  return (ef < zero) ? (minusOne * ef) : (ef);
}
/******************************************************************************************/

void ExtendedFloatTools::mult(MatrixXef &MatA, MatrixXef &MatB, MatrixXef &outMat){
  size_t nrMatA = MatA.rows();
  size_t ncMatB = MatB.cols();
  size_t nrMatB = MatB.rows();
  size_t ncMatA = MatA.cols();
  if (ncMatA != nrMatB) throw DimensionException("MatrixTools::mult(). nrows B != ncols A.", nrMatB, ncMatA);
  outMat.resize(nrMatA, ncMatB);
  for (size_t i = 0; i < nrMatA; i++){
    for (size_t j = 0; j < ncMatB; j++){
      outMat(i, j) = 0;
      for (size_t k = 0; k < ncMatA; k++){
        outMat(i, j) += MatA(i, k) * MatB(k, j);
      }
    }
  }
}
/******************************************************************************************/
void ExtendedFloatTools::copy(const MatrixXef &inMat, MatrixXef &outMat){
  size_t nrows = inMat.rows();
  size_t ncols = inMat.cols();
  outMat.resize(nrows, ncols);
  for (size_t i = 0; i < nrows; i++){
    for (size_t j = 0; j < ncols; j++){
      outMat(i, j) = inMat(i,j);
    }
  }
}
/******************************************************************************************/
ExtendedFloat ExtendedFloatTools::convertToExtendedFloat(double d){
    ExtendedFloat ef = ExtendedFloat{d};
    ef.normalize();
    return ef;
}