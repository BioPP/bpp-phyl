#ifndef DEFINITIONS_H
#define DEFINITIONS_H
#include <Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloatMatrix.h>

namespace bpp {
    // typedef ExtendedFloatMatrix<Eigen::Dynamic, Eigen::Dynamic> MatrixLik;
    // typedef ExtendedFloatRowVector RowLik;
    // typedef ExtendedFloatVector VectorLik;

    typedef Eigen::MatrixXd MatrixLik;
    typedef Eigen::RowVectorXd RowLik;
    typedef Eigen::VectorXd VectorLik;


    
}


#endif // DEFINITIONS_H