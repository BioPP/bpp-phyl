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

#include <Bpp/NewPhyl/LinearAlgebra.h>

namespace bpp {
// Utils
namespace {
	auto zeroValue (const Dimension<VectorDouble> & dim) -> decltype (VectorDouble::Zero (dim.size)) {
		return VectorDouble::Zero (dim.size);
	}
	auto zeroValue (const Dimension<MatrixDouble> & dim)
	    -> decltype (MatrixDouble::Zero (dim.rows, dim.cols)) {
		return MatrixDouble::Zero (dim.rows, dim.cols);
	}

	// Needed for numericProps, always false.
	constexpr bool isExactIdentity (const VectorDouble &) { return false; }

	// Generate a string describing useful numeric props of Vector/Matrix
	template <typename T> std::string numericProps (const T & t) {
		std::string s{"props{"};
		auto zero = zeroValue (dimensions (t));
		// Property for all elements
		if (isExactZero (t))
			s += "[0]";
		if (isExactOne (t))
			s += "[1]";
		if (isExactIdentity (t))
			s += "[I]";
		// Property on any element
		if (t.array ().isNaN ().any ())
			s += "N";
		if (t.array ().isInf ().any ())
			s += "i";
		if ((t.array () == zero.array ()).any ())
			s += "0";
		if ((t.array () > zero.array ()).any ())
			s += "+";
		if ((t.array () < zero.array ()).any ())
			s += "-";
		s += "}";
		return s;
	}
} // namespace

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

Dimension<double> dimensions (const DF::Value<double> &) noexcept {
	return {};
}
Dimension<VectorDouble> dimensions (const DF::Value<VectorDouble> & node) noexcept {
	return dimensions (node.accessValueConst ());
}
Dimension<MatrixDouble> dimensions (const DF::Value<MatrixDouble> & node) noexcept {
	return dimensions (node.accessValueConst ());
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

namespace DF {
	// Debug info specialisations
	template <> std::string Value<VectorDouble>::debugInfo () const {
		using std::to_string;
		auto & v = this->value_;
		return "size=" + to_string (dimensions (v)) + " " + numericProps (v);
	}
	template <> std::string Value<MatrixDouble>::debugInfo () const {
		using std::to_string;
		auto & m = this->value_;
		return "dim" + to_string (dimensions (m)) + " " + numericProps (m);
	}

	// Constant<VectorDouble> specialisation
	template <> NodeRef Constant<VectorDouble>::derive (const Node &) {
		return Builder<Constant<VectorDouble>>::makeZero (dimensions (*this));
	}
	template <> bool Constant<VectorDouble>::isDerivable (const Node &) const { return true; }
	std::shared_ptr<Constant<VectorDouble>>
	Builder<Constant<VectorDouble>>::makeZero (const Dimension<VectorDouble> & dim) {
		return make (zeroValue (dim));
	}
	std::shared_ptr<Constant<VectorDouble>>
	Builder<Constant<VectorDouble>>::makeOne (const Dimension<VectorDouble> & dim) {
		return make (VectorDouble::Ones (dim.size));
	}

	// Constant<MatrixDouble> specialisation
	template <> NodeRef Constant<MatrixDouble>::derive (const Node &) {
		return Builder<Constant<MatrixDouble>>::makeZero (dimensions (*this));
	}
	template <> bool Constant<MatrixDouble>::isDerivable (const Node &) const { return true; }
	std::shared_ptr<Constant<MatrixDouble>>
	Builder<Constant<MatrixDouble>>::makeZero (const Dimension<MatrixDouble> & dim) {
		return make (zeroValue (dim));
	}
	std::shared_ptr<Constant<MatrixDouble>>
	Builder<Constant<MatrixDouble>>::makeOne (const Dimension<MatrixDouble> & dim) {
		return make (MatrixDouble::Ones (dim.rows, dim.cols));
	}
} // namespace DF
} // namespace bpp
