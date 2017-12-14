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
class SubstitutionModel;

namespace DF {
	/** Data flow node representing a Model configured with parameter values.
	 * This class wraps a bpp::SubstitutionModel as a data flow node.
	 * It depends on Value<double> nodes (one for each parameter declared in the model).
	 * It provides a dummy value representing the "model configured by its parameters".
	 * This dummy value is then used by other node types to compute equilibrium frequencies,
	 * transition matrices and their derivatives.
	 *
	 * The dummy value is implemented as a pointer to the internal model for simplicity.
	 */
	class Model : public Value<const SubstitutionModel *> {
	public:
		/** Create a new model node from a dependency vector.
		 * Model parameters are given by a dependency vector of Value<double> nodes.
		 * The number and order of parameters is given by the SubstitutionModel internal ParameterList.
		 */
		Model (NodeRefVec && deps, std::unique_ptr<SubstitutionModel> && model);

		~Model ();

		SizeType nbParameters () const noexcept;

		// Legacy FIXME
		Model (std::unique_ptr<SubstitutionModel> model);
		ParameterRef<double> getParameter (SizeType index);
		ParameterRef<double> getParameter (const std::string & name);
		const std::string & getParameterName (SizeType index);

		std::string description () const override final;
		std::string debugInfo () const override final;

		// TODO  derivation (customizable)
		bool isDerivable (const Node & node) const override final;

	private:
		void compute () override final;
		std::unique_ptr<SubstitutionModel> model_;
	};

	/** Create a dependency vector suitable for a Model class constructor.
	 * The vector is built from the model internal parameter names, and Value<double> nodes in a map.
	 * For each named parameter in the model, a value node of the same node is taken from the map.
	 * Both namespaced and non-namespaced names are tried.
	 * If no node is found, an exception is thrown.
	 */
	NodeRefVec createDependencyVector (const SubstitutionModel & model,
	                                   const std::map<std::string, ValueRef<double>> & depsByName);

	// Compute nodes

	// (model) -> Vector of freqs
	class EquilibriumFrequenciesFromModel;
	template <> struct Builder<EquilibriumFrequenciesFromModel> {
		static ValueRef<VectorDouble> make (NodeRefVec && deps, SizeType nbStates); // FIXME dimension
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
