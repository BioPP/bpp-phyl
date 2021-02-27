#include <Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloatTools.h>


using namespace std;
using namespace bpp;




// void scale(MatrixXef& A, ExtendedFloat a, ExtendedFloat b = 0);
// void scale(MatrixXef& A, double a, double b = 0);
// void getIdentityMatrix(MatrixXef &outMat, size_t size){
//     outMat.resize(size, size);
//     for (size_t i = 0; i < size; i++){
//         for (size_t j = 0; j < size; j++){
//             if (i != j){
//                 outMat(i,j) = 0;
//                 outMat(i,j).normalize();
//             }else{
//                 outMat(i,j) = 1;
//                 outMat(i,j).normalize();
//             }
//         }
//     }
//     return;
// }
// /***************************************************************************/
// MatrixXef convertBppMatrixToEigenEF(Matrix <double>& mat){
//     size_t nCols = mat.getNumberOfColumns();
//     size_t nRows = mat.getNumberOfRows();
//     MatrixXef matEf(nRows, nCols);
//     for (size_t i = 0; i < nRows; i++){
//         for (size_t j = 0; j < nCols; j++){
//             ExtendedFloat ef = mat(i, j);
//             ef.normalize();
//             matEf(i, j) = ef;
//         }
//     }
//     return matEf;
// }
// /***************************************************************************/
// Eigen::MatrixXd convertBppMatrixToEigenDB(Matrix <double>& mat){
//     size_t nCols = mat.getNumberOfColumns();
//     size_t nRows = mat.getNumberOfRows();
//     Eigen::MatrixXd matEigen(nRows, nCols);
//     for (size_t i = 0; i < nRows; i++){
//         for (size_t j = 0; j < nCols; j++){
//             matEigen(i, j) = mat(i, j);
//         }
//     }
//     return matEigen;
// }
// /***************************************************************************/
// Eigen::MatrixXd convertMatToDouble (MatrixXef & eigenMatrix){
//     size_t nCols = eigenMatrix.cols();
//     size_t nRows = eigenMatrix.rows();
//     Eigen::MatrixXd dbMat(nRows, nCols);
//     for (size_t i = 0; i < nRows; i ++){
//         for (size_t j = 0; j < nCols; j++){
//             ExtendedFloat val = eigenMatrix(i,j);          
//             dbMat(i, j) = convert(val);
//         }
//     }
//     return dbMat;
// }

// void convertMatToDouble (MatrixXef & eigenMatrix, Eigen::MatrixXd & dbMat){
//     size_t nCols = eigenMatrix.cols();
//     size_t nRows = eigenMatrix.rows();
//     dbMat.resize(nRows, nCols);
//     for (size_t i = 0; i < nRows; i ++){
//         for (size_t j = 0; j < nCols; j++){
//             ExtendedFloat val = eigenMatrix(i,j);          
//             dbMat(i, j) = convert(val);
//         }
//     }
// }
// void convertDoubleMatToExtendedFloat(Eigen::MatrixXd & dbMat, MatrixXef & eigenMatrix){
//     size_t nCols = dbMat.cols();
//     size_t nRows = dbMat.rows();
//     eigenMatrix.resize(nRows, nCols);
//     for (size_t i = 0; i < nRows; i ++){
//         for (size_t j = 0; j < nCols; j++){
//             ExtendedFloat val = ExtendedFloat{dbMat(i,j)};
//             val.normalize();
//             eigenMatrix(i, j) = val;
           
//         }
//     }
// }
// ExtendedFloat exp(ExtendedFloat &ef){
//     double exp = convert(ef);
//     double e = std::exp(1);
//     ExtendedFloat e_ef = ExtendedFloat{e};
//     e_ef.normalize();
//     ExtendedFloat res = e_ef.pow(exp);
//     return res;

// }


// void copyToEigenVector(std::vector <double> & stdVector, VectorXef & eigenEfVector){
//     eigenEfVector.resize(stdVector.size());
//     for (size_t i = 0; i < stdVector.size(); i++){
//         eigenEfVector[i] = stdVector[i];
//         eigenEfVector[i].normalize();
//     }
// }
// void Taylor(const MatrixXef& A, size_t p, std::vector<MatrixXef> & vO){
//     size_t n = A.rows();
//     if (n != static_cast<size_t>(A.cols()))
//         throw DimensionException("test_ExtendedFloat::pow(). nrows != ncols.", static_cast<size_t>(A.cols()), static_cast<size_t>(A.rows()));
//     vO.resize(p+1);
//     getIdentityMatrix(vO[0],n);
//     vO[1] = A;   
//     for (size_t i = 1; i < p; i++)
//     {
//         //mult(vO[i], A, vO[i+1]);
//         vO[i+1] = vO[i] * A;
//     }
// }
// VectorXef exp(VectorXef v_in){
//     size_t nRows = v_in.rows();
//     VectorXef v_out(nRows);
//     for (size_t i = 0; i < nRows; i++){
//         v_out[i] = exp(v_in[i]);
//     }
//     return v_out;
// }

// /*
// * @param A [in] The first matrix.
// * @param D [in] The diagonal matrix (only diagonal elements in a vector)
// * @param U [in] The upper diagonal matrix (only upper diagonal elements in a vector)
// * @param L [in] The lower diagonal matrix (only lower diagonal elements in a vector)
// * @param B [in] The second matrix.
// * @param O [out] The result matrix.
// * @throw DimensionException If matrices have not the appropriate size.
// */
// void mult(const MatrixXef& A, const VectorXef& D, const VectorXef& U, const VectorXef& L, const MatrixXef& B, MatrixXef& O)
// {
//     size_t ncA = A.cols();
//     size_t nrA = A.rows();
//     size_t nrB = B.rows();
//     size_t ncB = B.cols();
//     if (ncA != nrB) throw DimensionException("MatrixTools::mult(). nrows B != ncols A.", nrB, ncA);
//     if (ncA != static_cast<size_t>(D.rows())) throw DimensionException("MatrixTools::mult(). Vector size is not equal to matrix size.", static_cast<size_t>(D.rows()), ncA);
//     if (ncA != static_cast<size_t>(U.rows())+1) throw DimensionException("MatrixTools::mult(). Vector size is not equal to matrix size-1.", static_cast<size_t>(U.rows()), ncA);
//     if (ncA != static_cast<size_t>(L.rows())+1) throw DimensionException("MatrixTools::mult(). Vector size is not equal to matrix size-1.", static_cast<size_t>(L.rows()), ncA);
//     O.resize(nrA, ncB);
//     for (size_t i = 0; i < nrA; i++)
//     {
//         for (size_t j = 0; j < ncB; j++)
//         {
//           O(i, j) = A(i, 0) * D[0] * B(0, j);
//           if (nrB>1)
//             O(i, j) += A(i,0) * U[0] * B(1,j);
//           for (size_t k = 1; k < ncA-1; k++)
//           {
//             O(i, j) += A(i, k) * (L[k-1] * B(k-1, j) + D[k] * B(k, j) + U[k] * B(k+1,j));
//           }
//           if (ncA>=2)
//             O(i, j) += A(i, ncA-1) * L[ncA-2] * B(ncA-2, j);
//           O(i,j) += A(i, ncA-1) * D[ncA-1] * B(ncA-1, j);
//         }
//     }
// }
// /**
// * @brief Add matrix x.B to matrix A.
// *
// * @param A [in,out] Matrix A
// * @param x [in] Scalar x
// * @param B [in] Matrix B
// * @throw DimensionException If A and B have note the same size.
// */
// void add(MatrixXef& A, ExtendedFloat& x, const MatrixXef& B)
// {
//     size_t ncA = A.cols();
//     size_t nrA = A.rows();
//     size_t nrB = B.rows();
//     size_t ncB = B.cols();
//     if (ncA != ncB) throw DimensionException("test_ExtendedFloat::operator+(). A and B must have the same number of columns.", ncB, ncA);
//     if (nrA != nrB) throw DimensionException("test_ExtendedFloat::operator+(). A and B must have the same number of rows.", nrB, nrA);
      
//     for (size_t i = 0; i < nrA; i++){
//         for (size_t j = 0; j < ncA; j++){
//           A(i, j) += x*B(i, j);
//         }
//     }
// }
// ExtendedFloat abs(ExtendedFloat ef){
//   ExtendedFloat zero = ExtendedFloat{0};
//   ExtendedFloat minusOne = ExtendedFloat{-1};
//   zero.normalize();
//   minusOne.normalize();
//   return (ef < zero) ? (minusOne * ef) : (ef);
// }


// void mult(MatrixXef &MatA, MatrixXef &MatB, MatrixXef &outMat){
//   size_t nrMatA = MatA.rows();
//   size_t ncMatB = MatB.cols();
//   size_t nrMatB = MatB.rows();
//   size_t ncMatA = MatA.cols();
//   if (ncMatA != nrMatB) throw DimensionException("MatrixTools::mult(). nrows B != ncols A.", nrMatB, ncMatA);
//   outMat.resize(nrMatA, ncMatB);
//   for (size_t i = 0; i < nrMatA; i++){
//     for (size_t j = 0; j < ncMatB; j++){
//       outMat(i, j) = 0;
//       for (size_t k = 0; k < ncMatA; k++){
//         outMat(i, j) += MatA(i, k) * MatB(k, j);
//       }
//     }
//   }
// }

// void copy(MatrixXef &inMat, MatrixXef &outMat){
//   size_t nrows = inMat.rows();
//   size_t ncols = inMat.cols();
//   outMat.resize(nrows, ncols);
//   for (size_t i = 0; i < nrows; i++){
//     for (size_t j = 0; j < ncols; j++){
//       outMat(i, j) = inMat(i,j);
//     }
//   }
// }

// void scale(MatrixXef& A, ExtendedFloat a, ExtendedFloat b){
//   size_t nrows = A.rows();
//   size_t ncols = A.cols();
//   for (size_t i = 0; i < nrows; i++){
//     for (size_t j = 0; j < ncols; j++){
//       A(i, j) = a * A(i, j) + b;
//     }
//   }
// }
// void scale(MatrixXef& A, double a, double b){
//   ExtendedFloat a_Ef = ExtendedFloat{a};
//   ExtendedFloat b_Ef = ExtendedFloat{b};
//   a_Ef.normalize();
//   b_Ef.normalize();
//   scale(A, a_Ef, b_Ef);

// }

// void updateEigenMatrices(
//     RowMatrix<double> &generator,
//     RowMatrix <double> &rightEigenVectors,
//     RowMatrix <double>& leftEigenVectors,
//     vector<double>& eigenValues,
//     vector <double>& iEigenValues,
//     vector <MatrixXef>& vPowExp,
//     bool isDiagonalizable,
//     bool isNonSingular)

// {
//     MatrixXef generator_Ef = convertBppMatrixToEigenEF(generator);
//     //in reality there is a data member eigenDecompose_, and the following line does not exist
//     bool eigenDecompose = true;

//   // Compute eigen values and vectors:
//   //if (enableEigenDecomposition())
//     if (eigenDecompose){
//     // Look for null lines (such as stop lines)
//     // ie null diagonal elements

//     size_t nbStop=0;
//     //size_t salph = getNumberOfStates();
//     size_t salph = generator_Ef.rows();
//     vector<bool> vnull(salph); // vector of the indices of lines with
//                                // only zeros

//     for (size_t i = 0; i < salph; i++)
//     {
//       if (abs(generator(i, i)) < NumConstants::TINY())
//       {
//         nbStop++;
//         vnull[i]=true;
//       }
//       else
//         vnull[i]=false;
//     }
        
//     if (nbStop != 0)
//     {
//       size_t salphok=salph - nbStop;
      
//       RowMatrix <double> gk(salphok, salphok);
//       size_t gi = 0, gj = 0;

//       for (size_t i = 0; i < salph; i++)
//       {
//         if (!vnull[i])
//         {
//           gj = 0;
//           for (size_t j = 0; j < salph; j++)
//           {
//             if (!vnull[j])
//             {
//               gk(i - gi, j - gj) = generator(i, j);
//             }
//             else
//               gj++;
//           }
//         }
//         else
//           gi++;
//       }

//       EigenValue<double> ev(gk);
//       eigenValues = ev.getRealEigenValues();
//       iEigenValues = ev.getImagEigenValues();

//       for (size_t i = 0; i < nbStop; i++)
//       {
//         eigenValues.push_back(0);
//         iEigenValues.push_back(0);
//       }

//       RowMatrix<double> rev = ev.getV();
//       rightEigenVectors.resize(salph, salph);
//       gi = 0;
//       for (size_t i = 0; i < salph; i++)
//       {
//         if (vnull[i])
//         {
//           gi++;
//           for (size_t j = 0; j < salph; j++)
//           {
//             rightEigenVectors(i, j) = 0;
//           }

//           rightEigenVectors(i, salphok + gi - 1) = 1;
//         }
//         else
//         {
//           for (size_t j = 0; j < salphok; j++)
//           {
//             rightEigenVectors(i, j) = rev(i - gi, j);
//           }

//           for (size_t j = salphok; j < salph; j++)
//           {
//             rightEigenVectors(i, j) = 0;
//           }
//         }
//       }
//     }
//     else
//     {
//       EigenValue<double> ev(generator);
//       rightEigenVectors = ev.getV();
//       eigenValues = ev.getRealEigenValues();
//       iEigenValues = ev.getImagEigenValues();
//       nbStop = 0;
//     }

//     /// Now check inversion and diagonalization
//     try
//     {
//       MatrixTools::inv(rightEigenVectors, leftEigenVectors);

//       // is it diagonalizable ?
//       isDiagonalizable = true;

//       //if (!dynamic_cast<ReversibleSubstitutionModel*>(this))
//       //{ we know it is irreversible model
//       for (auto& vi: iEigenValues){
//           if (abs(vi) > NumConstants::TINY()){
//               isDiagonalizable = false;
//               break;
//           }
//       }     
//       // looking for the vector of 0 eigenvalues

//       vector<size_t> vNullEv;
//       for (size_t i = 0; i< salph - nbStop; i++)
//         if ((abs(eigenValues[i]) < NumConstants::SMALL()) && (abs(iEigenValues[i]) < NumConstants::SMALL()))
//           vNullEv.push_back(i);
      

//       // pb to find unique null eigenvalue      
//       isNonSingular =(vNullEv.size()==1);

//       size_t nulleigen;
      
//       double val;
//       if (!isNonSingular)
//       {
//         //look or check which non-stop right eigen vector elements are
//         //equal.
//         for (auto cnull : vNullEv)
//         {
//           size_t i = 0;
//           while (vnull[i])
//             i++;
          
//           val = rightEigenVectors(i, cnull);
//           i++;
          
//           while (i < salph)
//           {
//             if (!vnull[i])
//             {
//               if (abs(rightEigenVectors(i, cnull) - val) > NumConstants::SMALL())
//                 break;
//             }
//             i++;
//           }
          
//           if (i >= salph)
//           {
//             isNonSingular = true;
//             nulleigen=cnull;
//             break;
//           }
//         }
//       }
//       else
//         nulleigen=vNullEv[0];
      
//       if (isNonSingular)
//       {
//         eigenValues[nulleigen] = 0; // to avoid approximation errors on long long branches
//         iEigenValues[nulleigen] = 0; // to avoid approximation errors on long long branches


//       }
//       else
//       {
//         //ApplicationTools::displayMessage("AbstractSubstitutionModel::updateMatrices : Unable to find eigenvector for eigenvalue 0. Taylor series used instead.");
//         isDiagonalizable = false;
//       }
//     }
//     // if rightEigenVectors_ is singular
//     catch (ZeroDivisionException& e)
//     {
//       //ApplicationTools::displayMessage("AbstractSubstitutionModel::updateMatrices : Singularity during diagonalization. Taylor series used instead.");
//       isNonSingular = false;
//       isDiagonalizable = false;
//     }

//     if (vPowExp.size() == 0)
//       vPowExp.resize(30);

      

//     getIdentityMatrix(vPowExp[0], salph);

//     //MatrixTools::Taylor(generator_, 30, vPowExp_);
//     Taylor(generator_Ef, 30, vPowExp);
//   }

// }



int main() {
   cout << "Just a test!!" << endl;
   Eigen::MatrixXd mat1(2,2);
   mat1(0,0) = 1.6;
   mat1(0,1) = 3;
   mat1(1,0) = 5;
   mat1(1,1) = 8;
   Eigen::MatrixXd mat2(2,2);
   mat2(0,0) = 4.9;
   mat2(0,1) = 10.2;
   mat2(1,0) = 0.2;
   mat2(1,1) = 0.08;

   Eigen::MatrixXd dbMatProd = mat1 * mat2;
   MatrixXef efMat1;
   MatrixXef efMat2;
   ExtendedFloatTools::convertDoubleMatToExtendedFloat (mat1,  efMat1);
   ExtendedFloatTools::convertDoubleMatToExtendedFloat (mat2,  efMat2);
   MatrixXef efMatProd = efMat1 * efMat2;
   Eigen::MatrixXd efMatProd_converted = ExtendedFloatTools::convertMatToDouble(efMatProd);
   cout << "Matrix double product: " << dbMatProd << endl;
   cout << "Matrix Ef product: " << efMatProd_converted << endl;
   cout << "Test Mult...." << endl;
   MatrixXef matTestMult;
   ExtendedFloatTools::mult(efMat1, efMat2, matTestMult);
   cout << "product using mult: " << ExtendedFloatTools::convertMatToDouble(matTestMult) << endl;
   MatrixXef matCopy;
   ExtendedFloatTools::copy(matTestMult, matCopy);
   cout << "testing copy: " << ExtendedFloatTools::convertMatToDouble(matCopy) << endl;
   
   // test exponent
   double expOfDB = std::exp(mat1(0,0));
   ExtendedFloat expEf = ExtendedFloatTools::exp(efMat1(0,0));
   double expEfConverted = convert(expEf);
   cout << "Testing exponent..." << endl;
   cout << "Exponent of double is: " << expOfDB << endl;
   cout << "Exponent of extended float is: " << expEfConverted << endl;

   cout << "Checking Taylor function ..." << endl;
   RowMatrix <double> generator;
   size_t size = 5;
   generator.resize(size, size);
   for (size_t i = 0; i < size; i++){
       for (size_t j = 0; j < size; j++){
           if (j == i-1){
               generator(i, j) = 1;

           }else if(j == i+1){
               generator(i, j) = 2;

           }else if (i == j){
               if (i == 0){
                   generator(i, j) = -2;
               }else if (i == size-1){
                   generator(i, j) = -1;
               }else{
                   generator(i, j) = -3;
               }
           }else{
               generator(i, j) = 0;
           }
           
       }
   }
   cout << "Generator double matrix:\n" << ExtendedFloatTools::convertBppMatrixToEigenDB(generator) << endl;
   cout << "***************" << endl;
   MatrixXef generator_ef = ExtendedFloatTools::convertBppMatrixToEigenEF(generator);
   vector <MatrixXef> vPowExp;
   vector <RowMatrix <double>> vPowExpDB;
   vPowExpDB.resize(30);
   vPowExp.resize(30);
   ExtendedFloatTools::Taylor(generator_ef, 30, vPowExp);
   cout << "Generator Ef matrix: \n" << ExtendedFloatTools::convertMatToDouble(generator_ef) << endl;
   cout << "******************" << endl;
   MatrixTools::Taylor(generator, 30, vPowExpDB);
   cout << "Taylor output double: \n" << ExtendedFloatTools::convertBppMatrixToEigenDB(vPowExpDB[29]) << endl;
   cout << "******************" << endl;
   cout << "Taylor output double: \n" << ExtendedFloatTools::convertMatToDouble(vPowExp[29]) << endl;

   cout << "Testing abs...." << endl;
   ExtendedFloat op1 = efMat1(0,0);
   ExtendedFloat op2 = efMat1(0,1);
   ExtendedFloat abs1 = ExtendedFloatTools::abs(op1-op2);
   ExtendedFloat abs2 = ExtendedFloatTools::abs(op2-op1);
   double absDb = fabs(mat1(0,0)-mat1(0,1));
   double abs1Db = abs1.convert(abs1);
   double abs2Db = abs2.convert(abs2);
   cout << "True abs: " << absDb << endl;
   cout << "abs1: " << abs1Db << endl;
   cout << "abs2: " << abs2Db << endl;
   cout << "************" << endl;
   cout << "Test scaling...." << endl;
   ExtendedFloatTools::scale(efMat1, 2);
   cout << "Mat after scaling: " << ExtendedFloatTools::convertMatToDouble(efMat1) << endl;
   cout << "Mat before scaling: " << mat1 << endl;
   

   return 0;
}
