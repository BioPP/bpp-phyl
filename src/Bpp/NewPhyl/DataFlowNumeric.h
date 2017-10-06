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
#include <Bpp/NewPhyl/DataFlowBuilder.h>
#include <Bpp/NewPhyl/DataFlowTemplateUtils.h>
#include <Bpp/NewPhyl/Debug.h> // description
#include <string>              // description
#include <type_traits>
#include <typeinfo>
#include <utility>

namespace bpp {
namespace DF {
	// Fwd declaration
	template <typename T> class Parameter;

	// Error functions
	[[noreturn]] void failureComputeWasCalled (const std::type_info & paramType);
	[[noreturn]] void failureDerivationNotSupportedForParameterType (const std::type_info & type);

	// Typedefs
	template <typename T> using ParameterRef = std::shared_ptr<Parameter<T>>;

	/* Numeric info.
	 *
	 * struct NumericProperties is used by node to give numeric properties about them.
	 * Compute class can override their numericProperties() method.
	 *
	 * T createZeroValue (const T&) is used to create a "0" for type T.
	 * This is useful for vector<T>, to create a vector of 0s with the same size as the source vector.
	 * TODO work in progress, may need refactoring
	 */
	struct NumericProperties {
		bool isConstant{false};
		bool isConstantZero{false};
		bool isConstantOne{false};
	};

	template <typename T> using IsParameterTypeDerivable = std::is_floating_point<T>;

	template <typename T, typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
	T createZeroValue (const T &) {
		return 0.;
	}
	template <typename T, typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
	T createOneValue (const T &) {
		return 1.;
	}
	template <typename T, typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
	bool isZeroValue (const T & v) {
		return v == T (0.);
	}
	template <typename T, typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
	bool isOneValue (const T & v) {
		return v == T (1.);
	}

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
			return createNode<Constant<T>> (createZeroValue (this->accessValue ()));
		}
		NumericProperties numericProperties () const override {
			auto props = Value<T>::numericProperties ();
			props.isConstant = true;
			props.isConstantZero = isZeroValue (this->accessValue ());
			props.isConstantOne = isOneValue (this->accessValue ());
			return props;
		}

		std::string description () const override final {
			return "Constant<" + prettyTypeName<T> () + ">(" + debug_to_string (this->accessValue ()) +
			       ")";
		}

		struct Builder : public AbstractBuilder {
			NodeRef build () const override { return createNode<Constant<T>> (constant); }

			T constant;
		};
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

	private:
		void compute () override final { failureComputeWasCalled (typeid (Parameter<T>)); }
	};

	/* Derive for parameters is special.
	 * It is only defined for real types (double, floats).
	 * Deriving a parameter should return a 1 if we derive from self, or 0 (derivation of constant).
	 * For other types, the function should throw an exception.
	 * Thus tag dispatching is done to route the call to derive() to the right implementation.
	 */
	template <typename T>
	NodeRef deriveParameterImpl (const Parameter<T> & parameter, const Node & variable,
	                             std::true_type) {
		if (&parameter == &variable)
			return createNode<Constant<T>> (createOneValue (parameter.accessValue ()));
		else
			return createNode<Constant<T>> (createZeroValue (parameter.accessValue ()));
	}
	template <typename T>
	NodeRef deriveParameterImpl (const Parameter<T> &, const Node &, std::false_type) {
		failureDerivationNotSupportedForParameterType (typeid (T));
	}
	template <typename T> NodeRef Parameter<T>::derive (const Node & variable) {
		return deriveParameterImpl (*this, variable, IsParameterTypeDerivable<T>{});
	}

	////////////////////////////////////// FIXME structures for double
	// Move away later !

	// Add double
	struct AddDouble : public Value<double> {
		using Dependencies = ReductionOfValue<double>;
		struct Builder;

		AddDouble (NodeRefVec && deps) : Value<double> (std::move (deps)) { checkDependencies (*this); }
		std::string description () const override final { return "+"; }
		void compute () override final {
			double a = 0.;
			for (const auto & dep : this->dependencies ())
				a += accessValueUnsafe<double> (*dep);
			this->value_ = a;
		}
	};
	struct AddDouble::Builder : public AbstractBuilder {
		NodeRef build () const override final {
			NodeRefVec deps;
			for (auto & builder : this->dependencies)
				deps.emplace_back (builder->build ());
			return createNode<AddDouble> (std::move (deps));
		}
	};

	// Multiply double
	struct MulDouble : public Value<double> {
		using Dependencies = ReductionOfValue<double>;
		struct Builder;

		MulDouble (NodeRefVec && deps) : Value<double> (std::move (deps)) { checkDependencies (*this); }
		std::string description () const final { return "*"; }
		void compute () override final {
			double a = 1.;
			for (const auto & dep : this->dependencies ())
				a *= accessValueUnsafe<double> (*dep);
			this->value_ = a;
		}
	};
	struct MulDouble::Builder : public AbstractBuilder {
		NodeRef build () const override final {
			NodeRefVec deps;
			for (auto & builder : this->dependencies)
				deps.emplace_back (builder->build ());
			return createNode<MulDouble> (std::move (deps));
		}
	};
}
}

#endif // BPP_NEWPHYL_DATAFLOWNUMERIC_H
