//
// File: ExtendedFloatTools.h
// Created by: Anat Shafir
// Created on: 2021
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

#ifndef _EXTENDEDFLOATTOOLS_H_
#define _EXTENDEDFLOATTOOLS_H_

#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowNumeric.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloat.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <algorithm>
#include <cassert>
#include <string>
#include <tuple>
#include <list>
#include <type_traits>
#include <Bpp/Exceptions.h>
#include <cmath>


namespace bpp
{
/**
 * @brief Tools to use the ExtendedFloat instance to prevent underflow issues
 */
typedef Eigen::Matrix<ExtendedFloat, Eigen::Dynamic, Eigen::Dynamic> MatrixXef;
typedef Eigen::Matrix<ExtendedFloat, Eigen::Dynamic, 1> VectorXef;
typedef Eigen::Matrix<ExtendedFloat, 1, Eigen::Dynamic> RowVectorXef;
class ExtendedFloatTools{



public:


  ExtendedFloatTools(){}
  virtual ~ExtendedFloatTools() {}
  
  
public:
    static ExtendedFloat convertToExtendedFloat(double d);
    static void scale(MatrixXef& A, ExtendedFloat a, ExtendedFloat b = 0);
    static void scale(MatrixXef& A, double a, double b = 0);
    static void getIdentityMatrix(MatrixXef &outMat, size_t size);
    static MatrixXef convertBppMatrixToEigenEF(const Matrix <double>& mat);
    static Eigen::MatrixXd convertBppMatrixToEigenDB(const Matrix <double>& mat);
    static Eigen::MatrixXd convertMatToDouble (const MatrixXef & eigenMatrix);
    static void convertMatToDouble (const MatrixXef & eigenMatrix, Eigen::MatrixXd & dbMat);
    static void convertDoubleMatToExtendedFloat(const Eigen::MatrixXd & dbMat, MatrixXef & eigenMatrix);
    static ExtendedFloat exp(ExtendedFloat &ef);
    static void copyToEigenVector(const std::vector <double> & stdVector, VectorXef & eigenEfVector);
    static void Taylor(const MatrixXef& A, size_t p, std::vector<MatrixXef> & vO);
    static VectorXef exp(VectorXef v_in);
    /*
    * @param A [in] The first matrix.
    * @param D [in] The diagonal matrix (only diagonal elements in a vector)
    * @param U [in] The upper diagonal matrix (only upper diagonal elements in a vector)
    * @param L [in] The lower diagonal matrix (only lower diagonal elements in a vector)
    * @param B [in] The second matrix.
    * @param O [out] The result matrix.
    * @throw DimensionException If matrices have not the appropriate size.
    */
    static void mult(const MatrixXef& A, const VectorXef& D, const VectorXef& U, const VectorXef& L, const MatrixXef& B, MatrixXef& O);
    /**
    * @brief Add matrix x.B to matrix A.
    *
    * @param A [in,out] Matrix A
    * @param x [in] Scalar x
    * @param B [in] Matrix B
    * @throw DimensionException If A and B have note the same size.
    */
    static void add(MatrixXef& A, ExtendedFloat& x, const MatrixXef& B);
    static void add(MatrixXef& A, double& x, const MatrixXef& B);
    static ExtendedFloat abs(ExtendedFloat ef);
    static void mult(MatrixXef &MatA, MatrixXef &MatB, MatrixXef &outMat);
    static void copy(const MatrixXef &inMat, MatrixXef &outMat);

};
} // end of namespace bpp.

#endif  // _EXTENDEDFLOATTOOLS_H_

