#include "ExtendedFloat.h"

namespace bpp {
    template<int R>
    class ExtendedFloatMatrix : public Eigen::Matrix<double, R, Eigen::Dynamic> {

        public:
            ExtendedFloatMatrix(void):Eigen::Matrix<double, R, Eigen::Dynamic>() {
                
            }
            ExtendedFloatMatrix(int cols):Eigen::Matrix<double, R, Eigen::Dynamic>(cols){

            }
            ExtendedFloatMatrix(int rows, int cols):Eigen::Matrix<double, R, Eigen::Dynamic>(rows, cols){
                
            }
            // This constructor allows you to construct ExtendedFloatMatrix from Eigen expressions
            template<typename D>
            ExtendedFloatMatrix(const Eigen::MatrixBase<D>& other): Eigen::Matrix<double, R, Eigen::Dynamic>(other){     
            }

            // This method allows you to assign Eigen expressions to ExtendedFloatMatrix
            template<typename D>
            ExtendedFloatMatrix& operator=(const Eigen::MatrixBase<D>& other){
                this->Eigen::Matrix<double, R, Eigen::Dynamic>::operator=(other);
                return *this;
            }



    };


}


