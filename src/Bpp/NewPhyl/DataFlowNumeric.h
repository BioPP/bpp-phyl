//
// File: DataFlowNumeric.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-09-15
// Last modified: 2017-09-15
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

#ifndef BPP_NEWPHYL_DATAFLOWNUMERIC_H
#define BPP_NEWPHYL_DATAFLOWNUMERIC_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Debug.h> // description
#include <string>              // description
#include <typeinfo>
#include <utility>

namespace bpp {
// FIXME improve (test with eigen::vec to decide impl)
template <typename T> struct NumericInfo { using Derivable = std::false_type; };
template <> struct NumericInfo<double> {
	using Derivable = std::true_type;
	static double zero () noexcept { return 0.0; }
	static double one () noexcept { return 1.0; }
};

namespace DF {
	// Fwd declaration
	template <typename T> class Parameter;

	// Typedefs
	template <typename T> using ParameterRef = std::shared_ptr<Parameter<T>>;

	struct NumericProperties {
		bool isConstant{false};
		bool isConstantZero{false};
	};

	/* Constant value.
	 */
	template <typename T> class Constant : public Value<T> {
	public:
		template <typename... Args>
		Constant (Args &&... args) : Value<T> (noDependency, std::forward<Args> (args)...) {
			this->makeValid ();
		}

		void compute () override final { failureComputeWasCalled (typeid (Constant<T>)); }

		// Deriving a constant returns 0
		NodeRef derive (const Node &) override final {
			return createNode<Constant<T>> (NumericInfo<T>::zero ());
		}
    void numericProperties (NumericProperties & props) const override {
      Value<T>::numericProperties (props);
      props.isConstant = true;
    }

		std::string description () const override final {
			return "Constant<" + prettyTypeName<T> () + ">(" + debug_to_string (this->value ()) + ")";
		}
	};

	/* Parameter node.
	 */
	template <typename T> class Parameter : public Value<T> {
	public:
		template <typename... Args>
		Parameter (Args &&... args) : Value<T> (noDependency, std::forward<Args> (args)...) {
			this->makeValid ();
		}

		void setValue (T t) noexcept {
			this->invalidate ();
			this->value_ = std::move (t);
			this->makeValid ();
		}

		NodeRef derive (const Node & variable) override final;

		std::string description () const final { return "Parameter<" + prettyTypeName<T> () + ">"; }

	private:
		void compute () override final { failureComputeWasCalled (typeid (Parameter<T>)); }
	};

	// Derivation only make sense for some types (like not for Sequence*)
	// Use template trick to only generate an error for unsupported types
	// FIXME improve too
	template <typename T>
	NodeRef deriveParameterHelper (const Parameter<T> & param, const Node & variable,
	                               std::true_type) {
		auto && v = (&variable == &param) ? NumericInfo<T>::one () : NumericInfo<T>::zero ();
		// FIXME one / zero depends on dynamic size for vectors !
		return createNode<Constant<T>> (v);
	}
	template <typename T>
	NodeRef deriveParameterHelper (const Parameter<T> &, const Node &, std::false_type) {
		failureDerivationNotSupportedForType (typeid (T));
	}
	template <typename T> NodeRef Parameter<T>::derive (const Node & variable) {
		return deriveParameterHelper<T> (*this, variable, typename NumericInfo<T>::Derivable{});
	}
}
}

#endif // BPP_NEWPHYL_DATAFLOWNUMERIC_H
