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
	// Declarations
	template <typename T> class Mutable;
	template <typename T> class Constant;
	template <typename T> using MutableRef = std::shared_ptr<Mutable<T>>;

	// Error function
	[[noreturn]] void failureComputeWasCalled (const std::type_info & nodeType);

	/** Mutable value node of type T.
	 * Leaf of the DataFlow graph, has no dependency.
	 * Represents a mutable value.
	 *
	 * A Mutable object always has a valid value (do not call invalidate !).
	 * It starts as initially valid.
	 * Modifications will send invalidations to ensure recomputations, then set it as valid again.
	 *
	 * Has no rebuild support: no dependency to change.
	 * Derivation is supported only for double.
	 */
	template <typename T> class Mutable : public Value<T> {
	public:
		/** Constructor.
		 * Forwards arguments to the T constructor.
		 * The T value is set as initially valid.
		 */
		template <typename... Args>
		Mutable (Args &&... args) : Value<T> (noDependency, std::forward<Args> (args)...) {
			this->makeValid ();
		}

		/** General case for modification of the T object.
		 * Takes a callable object (lamda, function pointer) that performs the modification.
		 * It must take a single T& as argument, which will refer to the T object to modify.
		 * The callable is called exactly once.
		 */
		template <typename Callable> void modify (Callable && modifier) {
			this->invalidateRecursively ();
			std::forward<Callable> (modifier) (this->accessValueMutable ());
			this->makeValid ();
		}

		/// Setter with invalidation.
		void setValue (const T & t) {
			modify ([&t](T & v) { v = t; });
		}
		/// Setter with invalidation (movable value version).
		void setValue (T && t) {
			modify ([&t](T & v) { v = std::move (t); });
		}

		// Defined as default to enable specialisation for some types
		NodeRef derive (const Node & node) final { return Value<T>::derive (node); }
		bool isDerivable (const Node & node) const final { return Value<T>::isDerivable (node); }

	private:
		void compute () final { failureComputeWasCalled (typeid (Mutable<T>)); }
	};

	/** Constant value of type T.
	 * DataFlow graph leaf, but for an immutable value, has no dependency.
	 * Is always valid.
	 *
	 * Has no rebuild support: no dependency to change.
	 * Derivation is supported for linear algebra types, always returns a zero T value.
	 */
	template <typename T> class Constant : public Value<T> {
	public:
		/** Constructor.
		 * Forwards arguments to the T constructor, which is called immediately.
		 */
		template <typename... Args>
		Constant (Args &&... args) : Value<T> (noDependency, std::forward<Args> (args)...) {
			this->makeValid ();
		}

		// Override info method
		bool isConstant () const final { return true; }

		// Defined as default to enable specialisation for some types
		NodeRef derive (const Node & node) final { return Value<T>::derive (node); }
		bool isDerivable (const Node & node) const final { return Value<T>::isDerivable (node); }

	private:
		void compute () final { failureComputeWasCalled (typeid (Constant<T>)); }
	};

	/* Specialisations in DataFlowNumeric.cpp
   * Mutable<double>: enable derivation.
   * Constant<double>: enable derivation.
   */
	template <> NodeRef Mutable<double>::derive (const Node & node);
	template <> bool Mutable<double>::isDerivable (const Node & node) const;

	template <> NodeRef Constant<double>::derive (const Node & node);
	template <> bool Constant<double>::isDerivable (const Node & node) const;
} // namespace DF
} // namespace bpp

#endif // BPP_NEWPHYL_DATAFLOWTEMPLATES_H
