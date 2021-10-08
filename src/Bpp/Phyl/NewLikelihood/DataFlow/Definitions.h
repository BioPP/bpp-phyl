//
// File: Definitions.h
// Authors:
//

#ifndef BPP_PHYL_NEWLIKELIHOOD_DATAFLOW_DEFINITIONS_H
#define BPP_PHYL_NEWLIKELIHOOD_DATAFLOW_DEFINITIONS_H

#include <Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloatEigen.h>


namespace bpp
{
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
#endif // BPP_PHYL_NEWLIKELIHOOD_DATAFLOW_DEFINITIONS_H
