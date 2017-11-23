//
// File: DataFlowTemplates.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-10-17
// Last modified: 2017-10-17
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

#ifndef BPP_NEWPHYL_DATAFLOWTEMPLATES_H
#define BPP_NEWPHYL_DATAFLOWTEMPLATES_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <memory>
#include <string>
#include <typeinfo>
#include <utility>

namespace bpp {
namespace DF {
	// TODO move Value (and maybe later value alternatives like SiteValue) here ?

	// Declarations
	template <typename T> class Parameter;
	template <typename T> class Constant;
	template <typename T> using ParameterRef = std::shared_ptr<Parameter<T>>;

	// Error function
	[[noreturn]] void failureComputeWasCalled (const std::type_info & nodeType);

	/* Parameter node.
	 * Leaf of the DataFlow graph, has no dependency.
	 * Represents a mutable value.
	 */
	template <typename T> class Parameter : public Value<T> {
	public:
		template <typename... Args>
		Parameter (Args &&... args) : Value<T> (noDependency, std::forward<Args> (args)...) {
			this->makeValid ();
		}

		template <typename Callable> void modify (Callable && modifier) {
			this->invalidateRecursively ();
			std::forward<Callable> (modifier) (this->value_);
			this->makeValid ();
		}
		void setValue (const T & t) {
			modify ([&t](T & v) { v = t; });
		}
		void setValue (T && t) {
			modify ([&t](T & v) { v = std::move (t); });
		}

		// Defined as default to enable specialisation
		NodeRef derive (const Node & node) override final { return Value<T>::derive (node); }
		bool isDerivable (const Node & node) const override final {
			return Value<T>::isDerivable (node);
		}

	private:
		void compute () override final { failureComputeWasCalled (typeid (Parameter<T>)); }
	};

	/* Constant node.
	 * DataFlow graph leaf, but for an immutable value.
	 */
	template <typename T> class Constant : public Value<T> {
	public:
		template <typename... Args>
		Constant (Args &&... args) : Value<T> (noDependency, std::forward<Args> (args)...) {
			this->makeValid ();
		}

		bool isConstant () const override final { return true; }

		// Defined as default to enable specialisation
		NodeRef derive (const Node & node) override final { return Value<T>::derive (node); }
		bool isDerivable (const Node & node) const override final {
			return Value<T>::isDerivable (node);
		}

	private:
		void compute () override final { failureComputeWasCalled (typeid (Constant<T>)); }
	};

	// Specialisations in DataFlow.cpp
	template <> NodeRef Parameter<double>::derive (const Node & node);
	template <> bool Parameter<double>::isDerivable (const Node & node) const;

	template <> NodeRef Constant<double>::derive (const Node & node);
	template <> bool Constant<double>::isDerivable (const Node & node) const;
	template <> struct Builder<Constant<double>> {
		static std::shared_ptr<Constant<double>> make (double d);
		static std::shared_ptr<Constant<double>> makeZero ();
		static std::shared_ptr<Constant<double>> makeOne ();
	};
} // namespace DF
} // namespace bpp

#endif // BPP_NEWPHYL_DATAFLOWTEMPLATES_H
