#include <Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloatTools.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloatMatrix.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowCWiseComputing.h>


using namespace std;
using namespace bpp;



int main() {
   Eigen::MatrixXd mat1(2,2);
   mat1(0,0) = 1.6;
   mat1(0,1) = 3;
   mat1(1,0) = 5;
   mat1(1,1) = 8;

   ExtendedFloatMatrix<Eigen::Dynamic, Eigen::Dynamic> EfMatrix;
   cout << "****************" << endl;
   EfMatrix = mat1;
   cout << "Testing assignment operator for matrix ..." << endl;
   cout << EfMatrix << endl;
   ExtendedFloatRowVector v(2);
   v(0) = 8;
   v(1) = 9;
   cout << "****************" << endl;
   cout << "Testing assignment operator for RowVector..." << endl;
   cout << v << endl;
   ExtendedFloatRowVector EfV = v;
   cout << EfV << endl;
   cout << "Testing assignment constructor for matrix..." << endl;
   ExtendedFloatMatrix<Eigen::Dynamic, Eigen::Dynamic> newMat(2,2);
   newMat(0,0) = 1;
   newMat(0,1) = 2;
   newMat(1,0) = 3;
   newMat(1,1) = 4;
   cout << newMat << endl;
   cout << "****************" << endl;
   cout << "Testing matrix operators and eigen specific functions..." << endl;
   size_t nrows = newMat.rows();
   size_t ncols = newMat.cols();
   cout << "Cols: " << ncols << " Rows: " << nrows << endl;
   //ExtendedFloatMatrix<Eigen::Dynamic> newMatTry = newMat;
   ExtendedFloatMatrix<Eigen::Dynamic, Eigen::Dynamic> addSum = newMat + EfMatrix;
   cout << "Sum :" << endl;
   cout << addSum << endl;
   cout << "Identity function:" << endl;
   ExtendedFloatMatrix<Eigen::Dynamic, Eigen::Dynamic> identity = ExtendedFloatMatrix<Eigen::Dynamic, Eigen::Dynamic>::Identity(2,2);
   cout << identity << endl;
   cout << "****************" << endl;
   cout <<"Testing constructor and eigen functions for Vector:" << endl;
   ExtendedFloatVector u1(2);
   u1[0] = 1.5;
   u1[1] = 3;
   u1 = u1 * 2;
   cout << "Constructor and * operator: " << endl;
   cout << u1 << endl;
   ExtendedFloatVector u2 = ExtendedFloatVector::Ones(4);
   cout <<"Ones function: "<< endl;
   cout << u2 << endl;
   Eigen::MatrixXd m_xd = newMat;
   cout << "Coverting to MatrixXd" << endl;
   cout << m_xd << endl;




//    ExtendedFloatMatrix EfMatrix = ExtendedFloatMatrix::Zero(2,2);
//    for (size_t i = 0; i < 2; i++){
//        for (size_t j = 0; j < 2; j++){
//            EfMatrix(i, j) = (double)(i+j);
//            //EfMatrix(i, j) = ExtendedFloat{(double)(i+j)};
//            //EfMatrix(i, j).normalize();
//        }
//    }
//    cout << EfMatrix << endl;
//    cout << "***************************************************" << endl;
//    cout <<"***   ***   ***   ***   ***   ***   ***   ***   ***" << endl;
//    Eigen::MatrixXd m1(3,3);
//    Eigen::MatrixXd m2(3,2);
//    Eigen::MatrixXd maxRes(3,2);
//    double val = 0;
//    for (size_t i = 0; i < 3; i++){
//        for (size_t j = 0; j < 3; j++){
//            val ++;
//            m1(i, j) = val;
//        }
//    }
//    val = 0;
//    for (size_t i = 0; i < 2; i++){
//        for (size_t j = 0; j < 3; j++){
//            val ++;
//            m2(j, i) = val;
//        }
//    }
//    cout << m1 << endl;
//    cout << m2 << endl;

//    for (size_t i = 0; i < 3; i++){
//        for (size_t j = 0; j < 2; j++){
//            Eigen::ArrayXd m1Arr = m1.row(i).array();
           
//            Eigen::ArrayXd m2Arr = m2.col(j).array();
           
//            cout << m1Arr << endl;
//            cout << "***" << endl;
//            cout << m2Arr << endl;
//            cout << "***" << endl;
//            Eigen::ArrayXd prod = m1Arr * m2Arr;
//            cout << prod << endl;
//            cout << "***" << endl;
//            maxRes(i, j) = prod.maxCoeff();
//        }
//        //cwise(maxRes) = cwise (m2) * cwise (m2);
//    }
//    cout << maxRes << endl;

//    cout << "Just a test!!" << endl;
//    Eigen::MatrixXd mat1(2,2);
//    mat1(0,0) = 1.6;
//    mat1(0,1) = 3;
//    mat1(1,0) = 5;
//    mat1(1,1) = 8;
//    Eigen::MatrixXd mat2(2,2);
//    mat2(0,0) = 4.9;
//    mat2(0,1) = 10.2;
//    mat2(1,0) = 0.2;
//    mat2(1,1) = 0.08;

//    Eigen::MatrixXd dbMatProd = mat1 * mat2;
//    MatrixXef efMat1;
//    MatrixXef efMat2;
//    ExtendedFloatTools::convertDoubleMatToExtendedFloat (mat1,  efMat1);
//    ExtendedFloatTools::convertDoubleMatToExtendedFloat (mat2,  efMat2);
//    MatrixXef efMatProd = efMat1 * efMat2;
//    Eigen::MatrixXd efMatProd_converted = ExtendedFloatTools::convertMatToDouble(efMatProd);
//    cout << "Matrix double product: " << dbMatProd << endl;
//    cout << "Matrix Ef product: " << efMatProd_converted << endl;
//    cout << "Test Mult...." << endl;
//    MatrixXef matTestMult;
//    ExtendedFloatTools::mult(efMat1, efMat2, matTestMult);
//    cout << "product using mult: " << ExtendedFloatTools::convertMatToDouble(matTestMult) << endl;
//    MatrixXef matCopy;
//    ExtendedFloatTools::copy(matTestMult, matCopy);
//    cout << "testing copy: " << ExtendedFloatTools::convertMatToDouble(matCopy) << endl;
   
//    // test exponent
//    double expOfDB = std::exp(mat1(0,0));
//    ExtendedFloat expEf = ExtendedFloatTools::exp(efMat1(0,0));
//    double expEfConverted = convert(expEf);
//    cout << "Testing exponent..." << endl;
//    cout << "Exponent of double is: " << expOfDB << endl;
//    cout << "Exponent of extended float is: " << expEfConverted << endl;

//    cout << "Checking Taylor function ..." << endl;
//    RowMatrix <double> generator;
//    size_t size = 5;
//    generator.resize(size, size);
//    for (size_t i = 0; i < size; i++){
//        for (size_t j = 0; j < size; j++){
//            if (j == i-1){
//                generator(i, j) = 1;

//            }else if(j == i+1){
//                generator(i, j) = 2;

//            }else if (i == j){
//                if (i == 0){
//                    generator(i, j) = -2;
//                }else if (i == size-1){
//                    generator(i, j) = -1;
//                }else{
//                    generator(i, j) = -3;
//                }
//            }else{
//                generator(i, j) = 0;
//            }
           
//        }
//    }
//    cout << "Generator double matrix:\n" << ExtendedFloatTools::convertBppMatrixToEigenDB(generator) << endl;
//    cout << "***************" << endl;
//    MatrixXef generator_ef = ExtendedFloatTools::convertBppMatrixToEigenEF(generator);
//    vector <MatrixXef> vPowExp;
//    vector <RowMatrix <double>> vPowExpDB;
//    vPowExpDB.resize(30);
//    vPowExp.resize(30);
//    ExtendedFloatTools::Taylor(generator_ef, 30, vPowExp);
//    cout << "Generator Ef matrix: \n" << ExtendedFloatTools::convertMatToDouble(generator_ef) << endl;
//    cout << "******************" << endl;
//    MatrixTools::Taylor(generator, 30, vPowExpDB);
//    cout << "Taylor output double: \n" << ExtendedFloatTools::convertBppMatrixToEigenDB(vPowExpDB[29]) << endl;
//    cout << "******************" << endl;
//    cout << "Taylor output double: \n" << ExtendedFloatTools::convertMatToDouble(vPowExp[29]) << endl;

//    cout << "Testing abs...." << endl;
//    ExtendedFloat op1 = efMat1(0,0);
//    ExtendedFloat op2 = efMat1(0,1);
//    ExtendedFloat abs1 = ExtendedFloatTools::abs(op1-op2);
//    ExtendedFloat abs2 = ExtendedFloatTools::abs(op2-op1);
//    double absDb = fabs(mat1(0,0)-mat1(0,1));
//    double abs1Db = abs1.convert(abs1);
//    double abs2Db = abs2.convert(abs2);
//    cout << "True abs: " << absDb << endl;
//    cout << "abs1: " << abs1Db << endl;
//    cout << "abs2: " << abs2Db << endl;
//    cout << "************" << endl;
//    cout << "Test scaling...." << endl;
//    ExtendedFloatTools::scale(efMat1, 2);
//    cout << "Mat after scaling: " << ExtendedFloatTools::convertMatToDouble(efMat1) << endl;
//    cout << "Mat before scaling: " << mat1 << endl;

//    cout << "***************************************************" << endl;
//    cout <<"***   ***   ***   ***   ***   ***   ***   ***   ***" << endl;

   

   return 0;
}
