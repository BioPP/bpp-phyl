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

#include <memory>
#include <tuple>
#include <type_traits>

#include <Bpp/NewPhyl/Debug.h>
#include <iostream>

namespace bpp {
namespace Expr {
	// Constant value
	template <typename T> struct Constant {
		using ArgumentTypes = std::tuple<>;
		using ResultType = const T &;

		constexpr Constant (const T & value) : value_ (value) {}

		constexpr ResultType compute () const noexcept { return value_; }
		constexpr bool isConstant () const noexcept { return true; }

		const T value_;
	};
	template <typename T> Constant<T> constant (const T & value) { return {value}; }

	// Reference to an externally provided value
	template <typename T> struct Ref {
		using ArgumentTypes = std::tuple<>;
		using ResultType = const T &;

		Ref (const T & value, bool isConstant = false) : value_ (value), isConstant_ (isConstant) {}

		ResultType compute () const noexcept { return value_; }
		bool isConstant () const noexcept { return isConstant_; }

		const T & value_;
		bool isConstant_;
	};
	template <typename T> Ref<T> ref (const T & value) { return {value}; }

	// Addition of two sub expressions
	template <typename Lhs, typename Rhs> struct Addition {
		using LhsResultType = typename Lhs::ResultType;
		using RhsResultType = typename Rhs::ResultType;

		using ArgumentTypes = std::tuple<LhsResultType, RhsResultType>;
		using ResultType = decltype (std::declval<LhsResultType> () + std::declval<RhsResultType> ());

		Addition (const Lhs & lhs, const Rhs & rhs) : lhs_ (lhs), rhs_ (rhs) {}

		ResultType compute () const { return lhs_.compute () + rhs_.compute (); }
		bool isConstant () const { return lhs_.isConstant () && rhs_.isConstant (); }

		const Lhs lhs_;
		const Rhs rhs_;
	};
	template <typename Lhs, typename Rhs>
	Addition<Lhs, Rhs> operator+ (const Lhs & lhs, const Rhs & rhs) {
		return {lhs, rhs};
	}

	// Assignment
	template <typename T> struct AbstractUpdateValue {
		virtual ~AbstractUpdateValue () = default;
		virtual void update (T & result) = 0;
		virtual bool isConstant () const = 0;
	};
	template <typename T, typename Expr> struct Assignment : public AbstractUpdateValue<T> {
		Assignment (const Expr & expr) : expr_ (expr) {}
		void update (T & result) override final { result = expr_.compute (); }
		bool isConstant () const override final { return expr_.isConstant (); }
		const Expr expr_;
	};

	/* Simplification.
   *
   * Encoded as a list of possible simplification with guards.
   * When simplifying, each sub expression
   */
	template <typename Expression, int opt_num> struct Simplified {
		using Type = Expression;
		static Type simplify (const Expression & expr) { return expr; }
		static bool guard (const Expression &) { return true; }
		using IsDefault = std::true_type;
	};
	template <typename Lhs, typename Rhs> struct Simplified<Addition<Lhs, Rhs>, 0> {
		// If Lhs == const (0) : Lhs + Rhs => Rhs
		using Type = Rhs;
		static Type simplify (const Addition<Lhs, Rhs> & add) { return add.rhs_; }
		static bool guard (const Addition<Lhs, Rhs> & add) {
			return add.lhs_.isConstant () && add.lhs_.compute () == 0;
		}
		using IsDefault = std::false_type;
	};
	template <typename Lhs, typename Rhs> struct Simplified<Addition<Lhs, Rhs>, 1> {
		// If Rhs == const (0) : Lhs + Rhs => Lhs
		using Type = Lhs;
		static Type simplify (const Addition<Lhs, Rhs> & add) { return add.lhs_; }
		static bool guard (const Addition<Lhs, Rhs> & add) {
			return add.rhs_.isConstant () && add.rhs_.compute () == 0;
		}
		using IsDefault = std::false_type;
	};

  /* TODO keep or trash ?
   *
   * Automatic simplification is complex and might be unmaintanable for normal people.
   */
	template <typename T> using UPAssign = std::unique_ptr<AbstractUpdateValue<T>>;
	template <typename T, typename Expr> UPAssign<T> make_assignment (const Expr & expr) {
		return UPAssign<T> (new Assignment<T, Expr> (expr));
	}

	namespace detail {
		template <typename T, typename Expr, int opt_num>
		UPAssign<T> make_assignment_selector (const Expr & expr);

		template <typename T, typename Expr, int opt_num>
		UPAssign<T> make_assignment_select_default (const Expr & expr, std::true_type) {
			return make_assignment<T> (expr); // Default case
		}
		template <typename T, typename Expr, int opt_num>
		UPAssign<T> make_assignment_select_default (const Expr & expr, std::false_type) {
			using Simpl = Simplified<Expr, opt_num>;
			if (Simpl::guard (expr))
				return make_assignment<T> (Simpl::simplify (expr));
			else
				return make_assignment_selector<T, Expr, opt_num + 1> (expr);
		}

		template <typename T, typename Expr, int opt_num>
		UPAssign<T> make_assignment_selector (const Expr & expr) {
			std::cout << "[selector] " << opt_num << " " << demangle (typeid (expr).name ()) << "\n";
			return make_assignment_select_default<T, Expr, opt_num> (
			    expr, typename Simplified<Expr, opt_num>::IsDefault{});
		}
	}
	template <typename T, typename Expr> UPAssign<T> make_simplified_assignment (const Expr & expr) {
		return detail::make_assignment_selector<T, Expr, 0> (expr);
	}
}
}

#endif // BPP_NEWPHYL_TEMPLATEEXPRESSION_H
