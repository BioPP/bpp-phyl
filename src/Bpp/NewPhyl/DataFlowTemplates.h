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
	template <typename T> struct Parameter;
	template <typename T> using ParameterRef = std::shared_ptr<Parameter<T>>;

	// Error function
	[[noreturn]] void failureComputeWasCalled (const std::type_info & nodeType);

	/* Parameter node.
	 */
	template <typename T> struct Parameter : public Value<T> {
		template <typename... Args>
		Parameter (Args &&... args) : Value<T> (noDependency, std::forward<Args> (args)...) {
			this->makeValid ();
		}

		void setValue (T t) noexcept {
			this->invalidate ();
			this->value_ = std::move (t);
			this->makeValid ();
		}

		void compute () override final { failureComputeWasCalled (typeid (Parameter<T>)); }

		template <typename... Args> static std::shared_ptr<Parameter<T>> create (Args &&... args) {
			return std::make_shared<Parameter<T>> (std::forward<Args> (args)...);
		}

		/* Specialisable functions.
		 * They can be specialised for specific T types (require existing declaration).
		 * We must provide a default version for the declaration: same as parent.
		 */
		NodeRef derive (const Node & node) override final { return Value<T>::derive (node); }
	};

	// Specialisations in DataFlowNumeric.cpp
	template <> NodeRef Parameter<double>::derive (const Node & node);

	/* Constant node declaration.
	 * Definitions are in their respective headers.
	 * TODO use member specialisation with forward declared wrapper types for vector/matrix
	 */
	// template <typename T> struct Constant;
} // namespace DF
} // namespace bpp

#endif // BPP_NEWPHYL_DATAFLOWTEMPLATES_H
