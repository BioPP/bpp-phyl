//
// File: PhylogenyTypes.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-11-27
// Last modified: 2017-11-27
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

#ifndef BPP_NEWPHYL_PHYLOGENYTYPES_H
#define BPP_NEWPHYL_PHYLOGENYTYPES_H

#include <Bpp/NewPhyl/LinearAlgebraFwd.h>
#include <Bpp/NewPhyl/Signed.h>
#include <cassert>
#include <string>

namespace bpp {
/* Likelihood probabilities (final or intermediate) are stored in a matrix.
 * Frequencies for site k are stored in column k.
 * TODO accessors ?
 */
using LikelihoodData = MatrixDouble;

// defines a Dimension<MatrixDouble> compatible struct.
struct LikelihoodDataDimension : public Dimension<MatrixDouble> {
	LikelihoodDataDimension (SizeType nbSitesArg, SizeType nbStatesArg) noexcept
	    : Dimension<MatrixDouble> (nbStatesArg, nbSitesArg) {}
	LikelihoodDataDimension (const Dimension<MatrixDouble> & matDim) noexcept
	    : Dimension<MatrixDouble> (matDim) {}

	SizeType nbStates () const noexcept { return rows; }
	SizeType nbSites () const noexcept { return cols; }
};
std::string to_string (const LikelihoodDataDimension & dim);

// defines a Dimension<MatrixDouble> compatible struct.
using TransitionMatrix = MatrixDouble;
struct TransitionMatrixDimension : public Dimension<MatrixDouble> {
	TransitionMatrixDimension (SizeType nbStatesArg) noexcept
	    : Dimension<MatrixDouble> (nbStatesArg, nbStatesArg) {}
	TransitionMatrixDimension (const Dimension<MatrixDouble> & matDim) noexcept
	    : Dimension<MatrixDouble> (matDim) {
		assert (rows == cols);
	}

	SizeType nbStates () const noexcept { return rows; }
};
std::string to_string (const TransitionMatrixDimension & dim);

} // namespace bpp

#endif // BPP_NEWPHYL_PHYLOGENYTYPES_H
