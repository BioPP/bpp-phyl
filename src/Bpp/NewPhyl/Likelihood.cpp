//
// File: Likelihood.cpp
// Authors:
// Created: 2017-06-06
// Last modified: 2017-06-06
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

#include <Bpp/Exceptions.h>
#include <Bpp/NewPhyl/Likelihood.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>

namespace bpp {
namespace dataflow {
	/* TODO add function to generate initial matrix
	 *
	 * use correct way to get initial state.
	 * Newlik code suggests: use an AlignedValuesContainer.
	 * Call getStateValueAt (site index, sequence index, model state tested).
	 *
	 * Models can have different state spaces for the same sequence.
	 * DFG creation should take a StateMap as argument.
	 * Check that every model used is valid with the StateMap (same size).
	 *
	 * Remove this node. Directly build Constant<MatrixDouble> in Phylogeny.
	 * Build these nodes for an AlignedValuesContainer, sequenceName, statemap.
	 * Rename SequenceNodeAccess to InitialLikelihoodNodeProvider.
	 * The sole impl should be built with StateMap ref at creation (fixed).
	 */

	std::unordered_map<std::string, std::shared_ptr<NumericMutable<double>>>
	createParameterMapForModel (Context & c, const TransitionModel & model) {
		const auto & modelParameters = model.getParameters ();
		const auto nbParameters = modelParameters.size ();
		std::unordered_map<std::string, std::shared_ptr<NumericMutable<double>>> map;
		for (std::size_t i = 0; i < nbParameters; ++i) {
			const auto & param = modelParameters[i];
			map.emplace (model.getParameterNameWithoutNamespace (param.getName ()),
			             NumericMutable<double>::create (c, param.getValue ()));
		}
		return map;
	}

	NodeRefVec
	createDependencyVector (const TransitionModel & model,
	                        const std::function<NodeRef (const std::string &)> & getParameter) {
		const auto & modelParameters = model.getParameters ();
		const auto nbParameters = modelParameters.size ();
		NodeRefVec deps (nbParameters);
		for (std::size_t i = 0; i < nbParameters; ++i) {
			std::string nonNamespacedName =
			    model.getParameterNameWithoutNamespace (modelParameters[i].getName ());
			auto dep = getParameter (nonNamespacedName);
			if (!dep) {
				throw Exception ("createDependencyVector (TransitionModel): model parameter not found: " +
				                 nonNamespacedName);
			}
			deps[i] = std::move (dep);
		}
		return deps;
	}

	// Model node

	std::shared_ptr<ConfiguredModel>
	ConfiguredModel::create (Context & c, NodeRefVec && deps,
	                         std::unique_ptr<TransitionModel> && model) {
		if (!model) {
			throw Exception ("ConfiguredModel(): nullptr TransitionModel");
		}
		// Check dependencies
		const auto nbParameters = model->getParameters ().size ();
		checkDependenciesNotNull (typeid (Self), deps);
		checkDependencyVectorSize (typeid (Self), deps, nbParameters);
		checkDependencyRangeIsValue<double> (typeid (Self), deps, 0, nbParameters);
		return std::make_shared<Self> (std::move (deps), std::move (model));
	}

	ConfiguredModel::ConfiguredModel (NodeRefVec && deps, std::unique_ptr<TransitionModel> && model)
	    : Value<const TransitionModel *> (std::move (deps), model.get ()),
	      model_ (std::move (model)) {}

	ConfiguredModel::~ConfiguredModel () = default;

	const std::string & ConfiguredModel::getParameterName (std::size_t index) {
		return model_->getParameters ()[index].getName ();
	}
	std::size_t ConfiguredModel::getParameterIndex (const std::string & name) {
		return static_cast<std::size_t> (model_->getParameters ().whichParameterHasName (name));
	}

	std::string ConfiguredModel::description () const { return "Model(" + model_->getName () + ")"; }
	std::string ConfiguredModel::debugInfo () const {
		return "nbState=" + std::to_string (model_->getAlphabet ()->getSize ());
	}

	NodeRef ConfiguredModel::recreate (Context & c, NodeRefVec && deps) const {
		auto m = Self::create (c, std::move (deps), std::unique_ptr<TransitionModel>{model_->clone ()});
		m->config = this->config; // Duplicate derivation config
		return m;
	}

	void ConfiguredModel::compute () {
		// Update each internal model bpp::Parameter with the dependency
		auto & modelParameters = model_->getParameters ();
		const auto nbParameters = this->nbDependencies ();
		for (std::size_t i = 0; i < nbParameters; ++i) {
			auto & v = accessValueConstCast<double> (*this->dependency (i));
			auto & p = modelParameters[i];
			if (p.getValue () != v) {
				// TODO improve bpp::Model interface to change values by index.
				model_->setParameterValue (model_->getParameterNameWithoutNamespace (p.getName ()), v);
			}
		}
	}

	// EquilibriumFrequenciesFromModel

	ValueRef<Eigen::RowVectorXd>
	EquilibriumFrequenciesFromModel::create (Context & c, NodeRefVec && deps,
	                                         const Dimension<Eigen::RowVectorXd> & dim) {
		checkDependenciesNotNull (typeid (Self), deps);
		checkDependencyVectorSize (typeid (Self), deps, 1);
		checkNthDependencyIs<ConfiguredModel> (typeid (Self), deps, 0);
		return std::make_shared<Self> (std::move (deps), dim);
	}

	EquilibriumFrequenciesFromModel::EquilibriumFrequenciesFromModel (
	    NodeRefVec && deps, const Dimension<Eigen::RowVectorXd> & dim)
	    : Value<Eigen::RowVectorXd> (std::move (deps)), targetDimension (dim) {}

	std::string EquilibriumFrequenciesFromModel::debugInfo () const {
		using namespace numeric;
		return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension);
	}

	NodeRef EquilibriumFrequenciesFromModel::derive (Context & c, const Node & node) {
		// d(equFreqs)/dn = sum_i d(equFreqs)/dx_i * dx_i/dn
		// In this case x_i are the model parameters, dependencies of the single dep of this node.
		auto dep = this->dependency (0);
		const auto & model = static_cast<const ConfiguredModel &> (*dep);
		auto parameterDim = Dimension<double>{}; // Dimension for model params
		NodeRefVec derivativeSumDeps;
		for (std::size_t i = 0; i < model.nbDependencies (); ++i) {
			// First compute dxi_dn. If this maps to a constant 0, do not compute df_dxi at all (costly).
			auto dxi_dn = model.dependency (i)->derive (c, node);
			if (!dxi_dn->isConstantAnd (NumericalProperty::Zero)) {
				auto buildFWithNewXi = [this, i, &model](Context & c, ValueRef<double> newDep,
				                                         const Dimension<T> & nodeDim) {
					// The sub-graph that will be replicated with shifted inputs is: equFreq -> model
					NodeRefVec newModelDeps = model.dependencies ();
					newModelDeps[i] = std::move (newDep);
					auto newModel = model.recreate (c, std::move (newModelDeps));
					return Self::create (c, {std::move (newModel)}, targetDimension);
				};
				auto df_dxi = generateNumericalDerivative<T, double> (
				    c, model.config, model.dependency (i), parameterDim, targetDimension, buildFWithNewXi);
				derivativeSumDeps.emplace_back (CWiseMul<T, std::tuple<T, double>>::create (
				    c, {std::move (df_dxi), std::move (dxi_dn)}, targetDimension));
			}
		}
		return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps));
	}

	void EquilibriumFrequenciesFromModel::compute () {
		const auto * model = accessValueConstCast<const TransitionModel *> (*this->dependency (0));
		const auto & freqsFromModel = model->getFrequencies ();
		auto & r = this->accessValueMutable ();
		r = Eigen::Map<const T> (freqsFromModel.data (),
		                         static_cast<Eigen::Index> (freqsFromModel.size ()));
	}
} // namespace dataflow
} // namespace bpp
