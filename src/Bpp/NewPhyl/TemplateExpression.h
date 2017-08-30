//
// File: TemplateExpression.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-08-30
// Last modified: 2017-08-30
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

#ifndef BPP_NEWPHYL_TEMPLATEEXPRESSION_H
#define BPP_NEWPHYL_TEMPLATEEXPRESSION_H

#include <type_traits>

struct Ops {
	using ArgumentTypes = std::tuple<double, double>;
	using ResultType = double;
	static void compute (ResultType & r, const double & a, const double & b) { r = a + b; }
};

namespace Expr {
template <typename T> struct Variable {
	using ArgumentTypes = std::tuple<>;
	using ResultType = T;

	Variable (const T & var) : var_ (var) {}

	void compute (ResultType & r) const noexcept { r = var_; }

	const T & var_;
};
template <typename T> Variable<T> var (const T & t) {
	return {t};
}

template <typename Lhs, typename Rhs> struct Addition {
	using LhsResultType = typename Lhs::ResultType;
	using RhsResultType = typename Rhs::ResultType;

	using ArgumentTypes = std::tuple<LhsResultType, RhsResultType>;
	using ResultType = decltype (std::declval<LhsResultType> () + std::declval<RhsResultType> ());

	Addition (const Lhs & lhs, const Rhs & rhs) : lhs_ (lhs), rhs_ (rhs) {}

	void compute (ResultType & r) const {
		LhsResultType lhs_r; // Pas pratique avec vector<...>.
		RhsResultType rhs_r; // Eigen: utilise des template expr aussi, donc devrait être ok.
		lhs_.compute (lhs_r);
		rhs_.compute (rhs_r);
		r = lhs_r + rhs_r;
	}

	const Lhs lhs_;
	const Rhs rhs_;
};
template <typename Lhs, typename Rhs>
Addition<Lhs, Rhs> operator+ (const Lhs & lhs, const Rhs & rhs) {
	return {lhs, rhs};
}
}

#endif // BPP_NEWPHYL_TEMPLATEEXPRESSION_H
