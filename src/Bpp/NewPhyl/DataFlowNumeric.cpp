//
// File: DataFlowNumeric.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-10-09 00:00:00
// Last modified: 2017-10-10
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
#include <Bpp/NewPhyl/DataFlowTemplateUtils.h>
#include <Bpp/NewPhyl/Range.h>
#include <memory>
#include <string>
#include <typeinfo>
#include <utility>

namespace bpp {
namespace DF {
	// ConstantDouble
	ConstantDouble::ConstantDouble (double d) : Value<double> (noDependency, d) {
		this->makeValid ();
	}
	void ConstantDouble::compute () { failureComputeWasCalled (typeid (ConstantDouble)); }
	std::string ConstantDouble::description () const {
		return "double(" + std::to_string (this->accessValue ()) + ")";
	}
	bool ConstantDouble::isConstant () const { return true; }
	NodeRef ConstantDouble::derive (const Node &) { return zero; }
	std::shared_ptr<ConstantDouble> ConstantDouble::zero = std::make_shared<ConstantDouble> (0.);
	std::shared_ptr<ConstantDouble> ConstantDouble::one = std::make_shared<ConstantDouble> (1.);
	std::shared_ptr<ConstantDouble> ConstantDouble::create (double d) {
		if (d == 0.) {
			return zero;
		} else if (d == 1.) {
			return one;
		} else {
			return std::make_shared<ConstantDouble> (d);
		}
	}

	// ParameterDouble
	ParameterDouble::ParameterDouble (double d) : Value<double> (noDependency, d) {
		this->makeValid ();
	}
	void ParameterDouble::compute () { failureComputeWasCalled (typeid (ParameterDouble)); }
	std::string ParameterDouble::description () const { return "Parameter<double>"; }
	NodeRef ParameterDouble::derive (const Node & node) {
		if (&node == static_cast<const Node *> (this)) {
			return ConstantDouble::one;
		} else {
			return ConstantDouble::zero;
		}
	}
	void ParameterDouble::setValue (double d) {
		this->invalidate ();
		this->value_ = d;
		this->makeValid ();
	}
	std::shared_ptr<ParameterDouble> ParameterDouble::create (double d) {
		return std::make_shared<ParameterDouble> (d);
	}

	// AddDouble
	AddDouble::AddDouble (NodeRefVec && deps) : Value<double> (std::move (deps)) {
		checkDependencies (*this);
	}
	void AddDouble::compute () {
		callWithValues (*this, [](double & r) { r = 0.; }, [](double & r, double d) { r += d; });
	}
	std::string AddDouble::description () const { return "double + double"; }
	NodeRef AddDouble::derive (const Node & node) {
		return AddDouble::create (this->dependencies ().map (
		    [&node](const NodeRef & dep) { return dep->derive (node); }));
	}
	ValueRef<double> AddDouble::create (NodeRefVec && deps) {
    checkDependencies<Dependencies> (deps, typeid (AddDouble));
		// Remove '0s' from deps
		removeDependenciesIf (deps, predicateIsConstantValueMatching<double> (0.));
		// Node choice
		if (deps.size () == 1) {
			return convertRef<Value<double>> (deps[0]);
		} else if (deps.size () == 0) {
			return ConstantDouble::zero;
		} else {
			return std::make_shared<AddDouble> (std::move (deps));
		}
	}

	// MulDouble
	MulDouble::MulDouble (NodeRefVec && deps) : Value<double> (std::move (deps)) {
		checkDependencies (*this);
	}
	void MulDouble::compute () {
		callWithValues (*this, [](double & r) { r = 1.; }, [](double & r, double d) { r *= d; });
	}
	std::string MulDouble::description () const { return "double * double"; }
	NodeRef MulDouble::derive (const Node & node) {
		NodeRefVec addDeps;
		for (auto i : bpp::index_range (this->dependencies ())) {
			NodeRefVec mulDeps = this->dependencies ();
			mulDeps[i] = this->dependencies ()[i]->derive (node);
			addDeps.emplace_back (MulDouble::create (std::move (mulDeps)));
		}
		return AddDouble::create (std::move (addDeps));
	}
	ValueRef<double> MulDouble::create (NodeRefVec && deps) {
    checkDependencies<Dependencies> (deps, typeid (MulDouble));
		// Return 0 if any dep is 0
		if (std::any_of (deps.begin (), deps.end (), predicateIsConstantValueMatching<double> (0.))) {
			return ConstantDouble::zero;
		}
		// Remove any 1s
		removeDependenciesIf (deps, predicateIsConstantValueMatching<double> (1.));
		// Node choice
		if (deps.size () == 1) {
			return convertRef<Value<double>> (deps[0]);
		} else if (deps.size () == 0) {
			return ConstantDouble::one;
		} else {
			return std::make_shared<MulDouble> (std::move (deps));
		}
	}
} // namespace DF
} // namespace bpp
