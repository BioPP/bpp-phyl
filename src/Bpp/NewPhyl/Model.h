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
#include <Bpp/NewPhyl/NodeSpecification.h>
#include <Bpp/NewPhyl/Signed.h>
#include <Eigen/Core>
#include <memory>
#include <string> // description

namespace bpp {
class SubstitutionModel;

namespace Phyl {
	using TransitionMatrix = Eigen::MatrixXd;
	using FrequencyVector = Eigen::VectorXd;

	// DF Node representing a model with its parameter (wrapper to master code).

	class ModelNode : public DF::Value<const SubstitutionModel *>::Impl {
		// TODO wrap SubstitutionModel in a Pimpl ModelValue class.
	public:
		ModelNode (std::unique_ptr<SubstitutionModel> model);
		~ModelNode ();

		SizeType nbParameters () const noexcept { return this->dependencies ().size (); }
		DF::Parameter<double> getParameter (SizeType index);
		DF::Parameter<double> getParameter (const std::string & name);
		const std::string & getParameterName (SizeType index);

		void compute () final;

		std::string description () const override;

	private:
		std::unique_ptr<SubstitutionModel> model_;
	};

	// Compute nodes

	struct ComputeEquilibriumFrequenciesFromModelOp {
		using ResultType = FrequencyVector;
		using ArgumentTypes = std::tuple<const SubstitutionModel *>;
		static void compute (FrequencyVector & freqs, const SubstitutionModel * model);
		static std::string description ();
		// Should init with freq vector size (nb states)
	};
	using ComputeEquilibriumFrequenciesFromModelNode =
	    DF::GenericFunctionComputation<ComputeEquilibriumFrequenciesFromModelOp>;

	struct ComputeTransitionMatrixFromModelOp {
		using ResultType = TransitionMatrix;
		enum { Model, BrLen };
		using ArgumentTypes = std::tuple<const SubstitutionModel *, double>;
		static void compute (TransitionMatrix & matrix, const SubstitutionModel * model, double brlen);
		static std::string description ();
		// Should init with (nb_state, nb_state)
	};
	using ComputeTransitionMatrixFromModelNode =
	    DF::GenericFunctionComputation<ComputeTransitionMatrixFromModelOp>;

	struct ComputeTransitionMatrixFirstDerivativeFromModelOp {
		using ResultType = TransitionMatrix;
		enum { Model, BrLen };
		using ArgumentTypes = std::tuple<const SubstitutionModel *, double>;
		static void compute (TransitionMatrix & matrix, const SubstitutionModel * model, double brlen);
		static std::string description ();
		// Should init with (nb_state, nb_state)
	};
	using ComputeTransitionMatrixFirstDerivativeFromModelNode =
	    DF::GenericFunctionComputation<ComputeTransitionMatrixFirstDerivativeFromModelOp>;

	struct ComputeTransitionMatrixSecondDerivativeFromModelOp {
		using ResultType = TransitionMatrix;
		enum { Model, BrLen };
		using ArgumentTypes = std::tuple<const SubstitutionModel *, double>;
		static void compute (TransitionMatrix & matrix, const SubstitutionModel * model, double brlen);
		static std::string description ();
		// Should init with (nb_state, nb_state)
	};
	using ComputeTransitionMatrixSecondDerivativeFromModelNode =
	    DF::GenericFunctionComputation<ComputeTransitionMatrixSecondDerivativeFromModelOp>;

	// Specs

	class ModelEquilibriumFrequenciesSpec
	    : public DF::NodeSpecAlwaysGenerate<ComputeEquilibriumFrequenciesFromModelNode> {
	public:
		ModelEquilibriumFrequenciesSpec (DF::Node modelParameter, SizeType nbStates);
		DF::NodeSpecificationVec computeDependencies () const;
		DF::Node buildNode (DF::NodeVec deps) const;

	private:
		DF::Node modelParameter_;
		SizeType nbStates_;
	};

	class ModelTransitionMatrixSpec
	    : public DF::NodeSpecAlwaysGenerate<ComputeTransitionMatrixFromModelNode> {
	public:
		ModelTransitionMatrixSpec (DF::Node modelParameter, DF::Node branchLengthParameter,
		                           SizeType nbStates);
		DF::NodeSpecificationVec computeDependencies () const;
		DF::Node buildNode (DF::NodeVec deps) const;

	private:
		DF::Node modelParameter_;
		DF::Node branchLengthParameter_;
		SizeType nbStates_;
	};

	class ModelTransitionMatrixFirstDerivativeSpec
	    : public DF::NodeSpecAlwaysGenerate<ComputeTransitionMatrixFirstDerivativeFromModelNode> {
	public:
		ModelTransitionMatrixFirstDerivativeSpec (DF::Node modelParameter,
		                                          DF::Node branchLengthParameter, SizeType nbStates);
		DF::NodeSpecificationVec computeDependencies () const;
		DF::Node buildNode (DF::NodeVec deps) const;

	private:
		DF::Node modelParameter_;
		DF::Node branchLengthParameter_;
		SizeType nbStates_;
	};

	class ModelTransitionMatrixSecondDerivativeSpec
	    : public DF::NodeSpecAlwaysGenerate<ComputeTransitionMatrixSecondDerivativeFromModelNode> {
	public:
		ModelTransitionMatrixSecondDerivativeSpec (DF::Node modelParameter,
		                                           DF::Node branchLengthParameter, SizeType nbStates);
		DF::NodeSpecificationVec computeDependencies () const;
		DF::Node buildNode (DF::NodeVec deps) const;

	private:
		DF::Node modelParameter_;
		DF::Node branchLengthParameter_;
		SizeType nbStates_;
	};
}
}

#endif // BPP_NEWPHYL_MODEL_H
