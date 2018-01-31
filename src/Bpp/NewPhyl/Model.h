//
// File: Model.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-29
// Last modified: 2017-05-29
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

#pragma once
#ifndef BPP_NEWPHYL_MODEL_H
#define BPP_NEWPHYL_MODEL_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowTemplates.h>
#include <Bpp/NewPhyl/LinearAlgebraFwd.h>
#include <Bpp/NewPhyl/PhylogenyTypes.h>
#include <Bpp/NewPhyl/Signed.h>
#include <map>
#include <memory>
#include <string>

namespace bpp {
class TransitionModel;

/** Interface for accessing model parameters by their names.
 * Associates Value<double> nodes to names.
 */
class ModelParameterAccessByName {
public:
	virtual ~ModelParameterAccessByName () = default;
	virtual DF::ValueRef<double> getModelParameter (const std::string & name) const = 0;
};

/** Impl for ModelParameterAccessByName: map of DF::Mutable<double>.
 * Represents the set of parameters for one TransitionModel.
 * Mutable nodes are created for each TransitionModel parameter.
 * They are registered under their non-namespaced names.
 * They can be changed to point to other Mutable<double> objects.
 */
class ModelParameterMap : public ModelParameterAccessByName {
public:
	/** Creates a new ModelParameterMap.
	 * One Mutable<double> object is created for each model parameter.
	 * It is registered in the map with its non-namespaced name.
	 * Mutable<double> nodes are initialized with values from the TransitionModel parameters.
	 * The reference to the model is only used in this constructor, not stored.
	 */
	ModelParameterMap (const TransitionModel & model);

	/// Impl of ModelParameterAccessByName, returns the parameter or throws.
	DF::ValueRef<double> getModelParameter (const std::string & name) const override;

	/** Access Mutable directly, or throws if not found (non-const version).
	 * Can change which Mutable is associated to this name.
	 * This is useful to have multiple models share a subset of parameter values.
	 */
	DF::MutableRef<double> & operator[] (const std::string & name);

	/** Access Mutable directly or throws if not found (const version).
	 *  Mutable cannot be changed, but its value can.
	 */
	const DF::MutableRef<double> & operator[] (const std::string & name) const;

	/// Direct map access (for info, debug, iteration).
	const std::map<std::string, DF::MutableRef<double>> & getMap () const {
		return mutableNodeByName_;
	}

private:
	std::map<std::string, DF::MutableRef<double>> mutableNodeByName_;
};

namespace DF {
	/** Create a dependency vector suitable for a Model class constructor.
	 * The vector is built from model parameter names, and Value<double> nodes in a map-like object.
	 * For each named parameter in the model, a value node of the same node is taken from the object.
	 * Only non-namespaced names are tried.
	 * If no node is found in the map-like object, an exception is thrown.
	 */
	NodeRefVec createDependencyVector (const TransitionModel & model,
	                                   const ModelParameterAccessByName & depsByName);

	/** Data flow node representing a Model configured with parameter values.
	 * This class wraps a bpp::TransitionModel as a data flow node.
	 * It depends on Value<double> nodes (one for each parameter declared in the model).
	 * It provides a dummy value representing the "model configured by its parameters".
	 * This dummy value is then used by other node types to compute equilibrium frequencies,
	 * transition matrices and their derivatives.
	 *
	 * The dummy value is implemented as a pointer to the internal model for simplicity.
	 */
	class Model : public Value<const TransitionModel *> {
	public:
		/// Internal constructor, see Builder<Model>::make for doc.
		Model (NodeRefVec && deps, std::unique_ptr<TransitionModel> && model);

		~Model ();

		// Access some node information TODO namespacing semantics !
		SizeType nbParameters () const noexcept;
		ValueRef<double> getParameter (SizeType index);
		ValueRef<double> getParameter (const std::string & name);
		const std::string & getParameterName (SizeType index);

		std::string description () const final;
		std::string debugInfo () const final;

		// TODO  derivation (customizable)
		bool isDerivable (const Node & node) const final;

		NodeRef rebuild (NodeRefVec && deps) const final;

	private:
		void compute () final;
		std::unique_ptr<TransitionModel> model_;
	};

	/// Always build Model nodes by using makeNode, which uses one of the make functions defined here.
	template <> struct Builder<Model> {
		/** Create a new model node from a dependency vector.
		 * Model parameters are given by a dependency vector of Value<double> nodes.
		 * The number and order of parameters is given by the TransitionModel internal ParameterList.
		 */
		static std::shared_ptr<Model> make (NodeRefVec && deps,
		                                    std::unique_ptr<TransitionModel> && model);

		/** Create a new model for an association from parameter names to Value<double> nodes.
		 * Internally, this builds a dependency vector using createDependencyVector.
		 * It will throw if some parameter names are node found in depsByName.
		 */
		static std::shared_ptr<Model> make (const ModelParameterAccessByName & depsByName,
		                                    std::unique_ptr<TransitionModel> && model);
	};

	// Compute nodes

	// (model) -> Vector of freqs
	class EquilibriumFrequenciesFromModel;
	template <> struct Builder<EquilibriumFrequenciesFromModel> {
    // Dim == nbStates FIXME add a dim type for Equ Freq ?
		static ValueRef<VectorDouble> make (NodeRefVec && deps, const Dimension<VectorDouble> & dim);
	};

	// (model, branch length) -> transition matrix
	class TransitionMatrixFromModel;
	template <> struct Builder<TransitionMatrixFromModel> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, const TransitionMatrixDimension & dim);
	};

	// (model, branch length) -> d(transition matrix)
	class TransitionMatrixFromModelBrlenDerivative;
	template <> struct Builder<TransitionMatrixFromModelBrlenDerivative> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, const TransitionMatrixDimension & dim);
	};

	// (model, branch length) -> d2(transition matrix)
	class TransitionMatrixFromModelBrlenSecondDerivative;
	template <> struct Builder<TransitionMatrixFromModelBrlenSecondDerivative> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, const TransitionMatrixDimension & dim);
	};
} // namespace DF
} // namespace bpp

#endif // BPP_NEWPHYL_MODEL_H
