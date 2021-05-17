#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloatEigen.h>

namespace bpp {

  typedef ExtendedFloatMatrixXd MatrixLik;
  typedef ExtendedFloatRowVectorXd RowLik;
  typedef ExtendedFloatVectorXd VectorLik;
  typedef ExtendedFloat DataLik;

   // typedef Eigen::MatrixXd MatrixLik;
   // typedef Eigen::RowVectorXd RowLik;
   // typedef Eigen::VectorXd VectorLik;
   // typedef double DataLik;

  typedef std::vector<DataLik> VDataLik;
  typedef std::vector<VDataLik> VVDataLik;
  
}


#endif // DEFINITIONS_H
