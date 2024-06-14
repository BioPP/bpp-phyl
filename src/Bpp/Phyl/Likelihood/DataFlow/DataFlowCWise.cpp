// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Exceptions.h>

#include "DataFlowCWise.h"

namespace bpp
{
// Precompiled instantiations of numeric nodes
template class CWiseFill<RowLik, double>;
template class CWiseFill<VectorLik, double>;
template class CWiseFill<MatrixLik, VectorLik>;
template class CWiseFill<MatrixLik, RowLik>;

template class CWisePattern<RowLik>;
template class CWisePattern<MatrixLik>;

template class CWiseMatching<RowLik, ReductionOf<RowLik>>;
template class CWiseMatching<MatrixLik, ReductionOf<MatrixLik>>;
template class CWiseMatching<Eigen::RowVectorXd, ReductionOf<double>>;
template class CWiseMatching<ExtendedFloatRowVectorXd, ReductionOf<ExtendedFloat>>;
template class CWiseMatching<MatrixLik, ReductionOf<RowLik>>;

template class CWiseCompound<MatrixLik, ReductionOf<RowLik>>;
template class CWiseCompound<MatrixLik, ReductionOf<VectorLik>>;
} // namespace bpp
