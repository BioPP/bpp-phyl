//
// File: Optimizer.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-07-07
// Last modified: 2017-07-07
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

#ifndef BPP_NEWPHYL_OPTIMIZER_H
#define BPP_NEWPHYL_OPTIMIZER_H

#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/Optional.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>
#include <unordered_map>
#include <utility>

namespace bpp {

// Wraps a DF::Parameter<double> data flow node as a bpp::Parameter
class DataFlowParameter : public Parameter {
public:
	DataFlowParameter (const std::string & name, std::shared_ptr<DF::ParameterDouble> existingParam)
	    : Parameter (name, {}), dfParam_ (std::move (existingParam)) {}
	DataFlowParameter (const std::string & name, double initialValue)
	    : DataFlowParameter (name, DF::ParameterDouble::create (initialValue)) {}

	// Parameter boilerplate
	DataFlowParameter * clone () const override { return new DataFlowParameter (*this); }

	// Override value access
	double getValue () const override { return dfParam_->accessValue (); }
	void setValue (double v) override { dfParam_->setValue (v); }
	// TODO care about the listeners

	const std::shared_ptr<DF::ParameterDouble> & getDataFlowParameter () const noexcept {
		return dfParam_;
	}

private:
	std::shared_ptr<DF::ParameterDouble> dfParam_;
};

/*
 * TODO use AbstractParametrizable (provides no way of adding parameters)
 */
class DataFlowFunction : public DerivableSecondOrder {
private:
	DF::ValueRef<double> dfFunction_;
	ParameterList variables_;

	struct StringPairHash {
		std::size_t operator() (const std::pair<std::string, std::string> & p) const {
			std::hash<std::string> strHash{};
			return strHash (p.first) ^ (strHash (p.second) << 1);
		}
	};
	mutable std::unordered_map<std::string, DF::ValueRef<double>> firstOrderDerivativeNodes_;
	mutable std::unordered_map<std::pair<std::string, std::string>, DF::ValueRef<double>,
	                           StringPairHash>
	    secondOrderDerivativeNodes_;

public:
	DataFlowFunction (DF::ValueRef<double> dfNode, const ParameterList & variables)
	    : dfFunction_ (std::move (dfNode)), variables_ (variables) {}

	// Boilerplate
	DataFlowFunction * clone () const override { return new DataFlowFunction (*this); }

	// bpp::Parametrizable (prefix unused FIXME?)
	bool hasParameter (const std::string & name) const override {
		return variables_.hasParameter (name);
	}
	const ParameterList & getParameters () const override { return variables_; }
	const Parameter & getParameter (const std::string & name) const override {
		return variables_.getParameter (name);
	}
	double getParameterValue (const std::string & name) const override {
		return variables_.getParameterValue (name);
	}
	void setAllParametersValues (const ParameterList & params) override {
		return variables_.setAllParametersValues (params);
	}
	void setParameterValue (const std::string & name, double value) override {
		variables_.setParameterValue (name, value);
	}
	void setParametersValues (const ParameterList & params) override {
		variables_.setParametersValues (params);
	}
	bool matchParametersValues (const ParameterList & params) override {
		return variables_.matchParametersValues (params);
	}
	std::size_t getNumberOfParameters () const override { return variables_.size (); }
	void setNamespace (const std::string &) override {}
	std::string getNamespace () const override { return {}; }
	std::string getParameterNameWithoutNamespace (const std::string & name) const override {
		return name;
	}

	// bpp::Function
	void setParameters (const ParameterList & params) override {
		variables_.setParametersValues (params);
	}
	double getValue () const override { return dfFunction_->getValue (); }

	// bpp::DerivableFirstOrder
	void enableFirstOrderDerivatives (bool) override {}
	bool enableFirstOrderDerivatives () const override { return true; }
	double getFirstOrderDerivative (const std::string & variable) const override {
		auto it = firstOrderDerivativeNodes_.find (variable);
		DF::ValueRef<double> node;
		if (it != firstOrderDerivativeNodes_.end ()) {
			node = it->second;
		} else {
			node = DF::convertRef<DF::Value<double>> (
			    dfFunction_->derive (*getDataFlowParameter (variable)));
			firstOrderDerivativeNodes_.insert (std::make_pair (variable, node));
		}
		return node->getValue ();
	}

	// bpp::DerivableSecondOrder
	void enableSecondOrderDerivatives (bool) override {}
	bool enableSecondOrderDerivatives () const override { return true; }
	double getSecondOrderDerivative (const std::string & variable) const override {
		return getSecondOrderDerivative (variable, variable);
	}
	double getSecondOrderDerivative (const std::string & variable1,
	                                 const std::string & variable2) const override {
		auto mapKey = std::make_pair (variable1, variable2);
		auto it = secondOrderDerivativeNodes_.find (mapKey);
		DF::ValueRef<double> node;
		if (it != secondOrderDerivativeNodes_.end ()) {
			node = it->second;
		} else {
			node =
			    DF::convertRef<DF::Value<double>> (dfFunction_->derive (*getDataFlowParameter (variable1))
			                                           ->derive (*getDataFlowParameter (variable2)));
			secondOrderDerivativeNodes_.insert (std::make_pair (mapKey, node));
		}
		return node->getValue ();
	}

	// Debug introspection FIXME
	Vector<DF::NamedNodeRef> getAllNamedNodes (const std::string & funcName) const {
		Vector<DF::NamedNodeRef> namedNodes;
		namedNodes.emplace_back (DF::NamedNodeRef{dfFunction_, funcName});
		for (auto & derivative : firstOrderDerivativeNodes_)
			namedNodes.emplace_back (
			    DF::NamedNodeRef{derivative.second, "d(" + funcName + ")/d(" + derivative.first + ")"});
		for (auto & derivative : secondOrderDerivativeNodes_)
			namedNodes.emplace_back (
			    DF::NamedNodeRef{derivative.second, "d2(" + funcName + ")/d(" + derivative.first.first +
			                                            ")d(" + derivative.first.second + ")"});
		return namedNodes;
	}

private:
	const std::shared_ptr<DF::ParameterDouble> &
	getDataFlowParameter (const std::string & name) const {
		return dynamic_cast<const DataFlowParameter &> (getParameter (name)).getDataFlowParameter ();
	}
};
} // namespace bpp

#endif // BPP_NEWPHYL_OPTIMIZER_H
