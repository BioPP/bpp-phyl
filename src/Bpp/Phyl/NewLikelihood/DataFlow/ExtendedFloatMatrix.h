#ifndef EXTENDEDFLOATMATRIX_H
#define EXTENDEDFLOATMATRIX_H


#include "ExtendedFloat.h"

namespace bpp {
    template<int R, int C>
    class ExtendedFloatMatrix : public Eigen::Matrix<double, R, C> {

        public:
            ExtendedFloatMatrix(void):Eigen::Matrix<double, R, C>() {
                
            }
            ExtendedFloatMatrix(int rows, int cols):Eigen::Matrix<double, R, C>(rows, cols){

            }
            ExtendedFloatMatrix(int cols):Eigen::Matrix<double, R, C>(cols){
                
            }

            // This constructor allows you to construct ExtendedFloatMatrix from Eigen expressions
            template<typename D>
            ExtendedFloatMatrix(const Eigen::MatrixBase<D>& other): Eigen::Matrix<double, R, C>(other){     
            }

            // This method allows you to assign Eigen expressions to ExtendedFloatMatrix
            template<typename D>
            ExtendedFloatMatrix& operator=(const Eigen::MatrixBase<D>& other){
                this->Eigen::Matrix<double, R, C>::operator=(other);
                return *this;
            }

    };
    // An extension for VectorXd
    class ExtendedFloatVector : public ExtendedFloatMatrix<Eigen::Dynamic, 1>{
        public:
            ExtendedFloatVector(void):ExtendedFloatMatrix<Eigen::Dynamic, 1>() {
                
            }
            ExtendedFloatVector(int rows):ExtendedFloatMatrix<Eigen::Dynamic, 1>(rows){
                
            }

            // This constructor allows you to construct ExtendedFloatMatrix from Eigen expressions
            template<typename D>
            ExtendedFloatVector(const Eigen::MatrixBase<D>& other): ExtendedFloatMatrix<Eigen::Dynamic, 1>(other){     
            }

            // This method allows you to assign Eigen expressions to ExtendedFloatMatrix
            template<typename D>
            ExtendedFloatVector& operator=(const Eigen::MatrixBase<D>& other){
                this->ExtendedFloatMatrix<Eigen::Dynamic, 1>::operator=(other);
                return *this;
            }
    };

    // An extension to RowVectorXd

    class ExtendedFloatRowVector: public ExtendedFloatMatrix<1, Eigen::Dynamic>{
        public:
        ExtendedFloatRowVector(void):ExtendedFloatMatrix<1, Eigen::Dynamic>() {
                
        }
        ExtendedFloatRowVector(int cols):ExtendedFloatMatrix<1, Eigen::Dynamic>(cols){
                
        }

        // This constructor allows you to construct ExtendedFloatMatrix from Eigen expressions
        template<typename D>
        ExtendedFloatRowVector(const Eigen::MatrixBase<D>& other): ExtendedFloatMatrix<1, Eigen::Dynamic>(other){     
        }

        // This method allows you to assign Eigen expressions to ExtendedFloatMatrix
        template<typename D>
        ExtendedFloatRowVector& operator=(const Eigen::MatrixBase<D>& other){
            this->ExtendedFloatMatrix<1, Eigen::Dynamic>::operator=(other);
            return *this;
        }

    };


}

#endif // EXTENDEDFLOATMATRIX_H
