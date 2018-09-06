//
// File: Model.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-29
// Last modified: 2017-05-29
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#pragma once
#ifndef BPP_NEWPHYL_MODEL_H
#define BPP_NEWPHYL_MODEL_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowTemplates.h>
#include <Bpp/NewPhyl/LinearAlgebraFwd.h>
#include <Bpp/NewPhyl/PhylogenyTypes.h>
#include <map>
#include <memory>
#include <string>

namespace bpp {
class TransitionModel;

namespace DF {
	// Compute nodes

	// (model, branch length) -> d(transition matrix)
	class TransitionMatrixFromModelBrlenDerivative;
	template <> struct Builder<TransitionMatrixFromModelBrlenDerivative> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, const TransitionMatrixDimension & dim);
	};

	// (model, branch length) -> d2(transition matrix)
	class TransitionMatrixFromModelBrlenSecondDerivative;
	template <> struct Builder<TransitionMatrixFromModelBrlenSecondDerivative> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, const TransitionMatrixDimension & dim);
	};
} // namespace DF
} // namespace bpp

#endif // BPP_NEWPHYL_MODEL_H
