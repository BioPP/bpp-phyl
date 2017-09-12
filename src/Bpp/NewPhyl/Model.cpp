//
// File: Model.cpp
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

#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/Model.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <cassert>

namespace bpp {
namespace Phyl {
	// Model DF Node

	ModelNode::ModelNode (std::unique_ptr<SubstitutionModel> model)
	    : DF::Value<const SubstitutionModel *> (DF::noDependency, model.get ()),
	      model_ (std::move (model)) {
		// TODO support already built paramater nodes
		model_->setNamespace ({}); // Delete namespace prefix
		const auto & parameters = model_->getParameters ();
		for (auto i : index_range (parameters))
			this->appendDependency (DF::createNode<DF::Parameter<double>> (parameters[i].getValue ()));
		this->makeValid (); // Initially valid
	}

	ModelNode::~ModelNode () = default;

	DF::ParameterRef<double> ModelNode::getParameter (SizeType index) {
		return DF::convertRef<DF::Parameter<double>> (this->dependencies ().at (index));
	}
	DF::ParameterRef<double> ModelNode::getParameter (const std::string & name) {
		return getParameter (
		    static_cast<SizeType> (model_->getParameters ().whichParameterHasName (name)));
	}
	const std::string & ModelNode::getParameterName (SizeType index) {
		return model_->getParameters ()[static_cast<std::size_t> (index)].getName ();
	}

	void ModelNode::compute () {
		// Update current model params with ours
		auto & modelParams = model_->getParameters ();
		for (auto i : index_range (this->dependencies ())) {
			auto v = DF::accessValueUnsafe<double> (*this->dependencies ()[i]);
			auto & p = modelParams[static_cast<std::size_t> (i)];
			if (p.getValue () != v)
				model_->setParameterValue (p.getName (), v);
		}
	}

	std::string ModelNode::description () const { return "Model(" + model_->getName () + ")"; }

	// Compute node functions

	void ComputeEquilibriumFrequenciesFromModelOp::compute (FrequencyVector & freqs,
	                                                        const SubstitutionModel * model) {
		auto & freqsFromModel = model->getFrequencies ();
		freqs = Eigen::Map<const FrequencyVector> (freqsFromModel.data (),
		                                           Eigen::Index (freqsFromModel.size ()));
	}

	std::string ComputeEquilibriumFrequenciesFromModelOp::description () {
		return "EquilibriumFreqs";
	}

	namespace {
		/* For now copy matrix cell by cell.
		 * TODO use eigen internally in SubsitutionModel ! (not perf critical for now though)
		 * FIXME if multithreading, internal model state may be a problem !
		 */
		void bppToEigen (const Matrix<double> & bppMatrix, TransitionMatrix & eigenMatrix) {
			assert (eigenMatrix.rows () == bppMatrix.getNumberOfRows ());
			assert (eigenMatrix.cols () == bppMatrix.getNumberOfColumns ());
			for (auto i : range (eigenMatrix.rows ()))
				for (auto j : range (eigenMatrix.cols ()))
					eigenMatrix (i, j) =
					    bppMatrix (static_cast<std::size_t> (i), static_cast<std::size_t> (j));
		}
	}

	void ComputeTransitionMatrixFromModelOp::compute (TransitionMatrix & matrix,
	                                                  const SubstitutionModel * model, double brlen) {
		bppToEigen (model->getPij_t (brlen), matrix);
	}

	std::string ComputeTransitionMatrixFromModelOp::description () { return "TransitionMatrix"; }

	void ComputeTransitionMatrixFirstDerivativeFromModelOp::compute (TransitionMatrix & matrix,
	                                                                 const SubstitutionModel * model,
	                                                                 double brlen) {
		bppToEigen (model->getdPij_dt (brlen), matrix);
	}

	std::string ComputeTransitionMatrixFirstDerivativeFromModelOp::description () {
		return "d(TransitionMatrix)/dt";
	}

	void ComputeTransitionMatrixSecondDerivativeFromModelOp::compute (TransitionMatrix & matrix,
	                                                                  const SubstitutionModel * model,
	                                                                  double brlen) {
		bppToEigen (model->getd2Pij_dt2 (brlen), matrix);
	}

	std::string ComputeTransitionMatrixSecondDerivativeFromModelOp::description () {
		return "d2(TransitionMatrix)/dt2";
	}

	// Specs

	ModelEquilibriumFrequenciesSpec::ModelEquilibriumFrequenciesSpec (DF::NodeRef modelParameter,
	                                                                  SizeType nbStates)
	    : modelParameter_ (std::move (modelParameter)), nbStates_ (nbStates) {}
	DF::NodeSpecificationVec ModelEquilibriumFrequenciesSpec::computeDependencies () const {
		return DF::makeNodeSpecVec (DF::NodeSpecReturnParameter (modelParameter_));
	}
	DF::NodeRef ModelEquilibriumFrequenciesSpec::buildNode (DF::NodeRefVec deps) const {
		return DF::createNode<ComputeEquilibriumFrequenciesFromModelNode> (std::move (deps), nbStates_);
	}

	ModelTransitionMatrixSpec::ModelTransitionMatrixSpec (DF::NodeRef modelParameter,
	                                                      DF::NodeRef branchLengthParameter,
	                                                      SizeType nbStates)
	    : modelParameter_ (std::move (modelParameter)),
	      branchLengthParameter_ (std::move (branchLengthParameter)),
	      nbStates_ (nbStates) {}
	DF::NodeSpecificationVec ModelTransitionMatrixSpec::computeDependencies () const {
		return DF::makeNodeSpecVec (DF::NodeSpecReturnParameter (modelParameter_),
		                            DF::NodeSpecReturnParameter (branchLengthParameter_));
	}
	DF::NodeRef ModelTransitionMatrixSpec::buildNode (DF::NodeRefVec deps) const {
		return DF::createNode<ComputeTransitionMatrixFromModelNode> (std::move (deps), nbStates_,
		                                                             nbStates_);
	}

	ModelTransitionMatrixFirstDerivativeSpec::ModelTransitionMatrixFirstDerivativeSpec (
	    DF::NodeRef modelParameter, DF::NodeRef branchLengthParameter, SizeType nbStates)
	    : modelParameter_ (std::move (modelParameter)),
	      branchLengthParameter_ (std::move (branchLengthParameter)),
	      nbStates_ (nbStates) {}
	DF::NodeSpecificationVec ModelTransitionMatrixFirstDerivativeSpec::computeDependencies () const {
		return DF::makeNodeSpecVec (DF::NodeSpecReturnParameter (modelParameter_),
		                            DF::NodeSpecReturnParameter (branchLengthParameter_));
	}
	DF::NodeRef ModelTransitionMatrixFirstDerivativeSpec::buildNode (DF::NodeRefVec deps) const {
		return DF::createNode<ComputeTransitionMatrixFirstDerivativeFromModelNode> (
		    std::move (deps), nbStates_, nbStates_);
	}

	ModelTransitionMatrixSecondDerivativeSpec::ModelTransitionMatrixSecondDerivativeSpec (
	    DF::NodeRef modelParameter, DF::NodeRef branchLengthParameter, SizeType nbStates)
	    : modelParameter_ (std::move (modelParameter)),
	      branchLengthParameter_ (std::move (branchLengthParameter)),
	      nbStates_ (nbStates) {}
	DF::NodeSpecificationVec ModelTransitionMatrixSecondDerivativeSpec::computeDependencies () const {
		return DF::makeNodeSpecVec (DF::NodeSpecReturnParameter (modelParameter_),
		                            DF::NodeSpecReturnParameter (branchLengthParameter_));
	}
	DF::NodeRef ModelTransitionMatrixSecondDerivativeSpec::buildNode (DF::NodeRefVec deps) const {
		return DF::createNode<ComputeTransitionMatrixSecondDerivativeFromModelNode> (
		    std::move (deps), nbStates_, nbStates_);
	}
}
}
