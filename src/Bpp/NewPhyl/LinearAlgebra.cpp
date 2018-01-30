//
// File: LinearAlgebra.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-10-27
// Last modified: 2017-10-27
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

#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/NewPhyl/LinearAlgebra.h>
#include <Bpp/NewPhyl/LinearAlgebraUtils.h>

namespace bpp {
// Dimensions
std::string to_string (const Dimension<VectorDouble> & dim) {
	return std::to_string (dim.size);
}
std::string to_string (const Dimension<MatrixDouble> & dim) {
	return "(" + std::to_string (dim.rows) + "," + std::to_string (dim.cols) + ")";
}

Dimension<VectorDouble> dimensions (const VectorDouble & v) noexcept {
	return {v.rows ()};
}
Dimension<MatrixDouble> dimensions (const MatrixDouble & m) noexcept {
	return {m.rows (), m.cols ()};
}

// Predicates
bool isExactZero (const double & d) {
	return d == 0.;
}
bool isExactZero (const VectorDouble & v) {
	return v.isZero (0.);
}
bool isExactZero (const MatrixDouble & m) {
	return m.isZero (0.);
}
bool isExactOne (const double & d) {
	return d == 1.;
}
bool isExactOne (const VectorDouble & v) {
	return v.isOnes (0.);
}
bool isExactOne (const MatrixDouble & m) {
	return m.isOnes (0.);
}
bool isExactIdentity (const MatrixDouble & m) {
	return m.rows () == m.cols () && m.isIdentity (0.);
}
} // namespace bpp
