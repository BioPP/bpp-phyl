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

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/Model.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <iostream>

namespace bpp {
namespace Phyl {
	ModelNode::ModelNode (std::unique_ptr<SubstitutionModel> model)
	    : DF::Value<const SubstitutionModel *>::Impl (model.get ()), model_ (std::move (model)) {

		model_->setNamespace ({}); // Delete namespace prefix
		const auto & parameters = model_->getParameters ();
		for (auto i : index_range (parameters))
			this->appendDependency (DF::Parameter<double>::create (parameters[i].getValue ()));

		this->makeValid (); // Initially valid
	}

	ModelNode::~ModelNode () = default;

	DF::Parameter<double> ModelNode::getParameter (std::size_t index) {
		return DF::Parameter<double>{this->dependencies ().at (index)};
	}
	DF::Parameter<double> ModelNode::getParameter (const std::string & name) {
		return getParameter (model_->getParameters ().whichParameterHasName (name));
	}
	const std::string & ModelNode::getParameterName (std::size_t index) {
		return model_->getParameters ()[index].getName ();
	}

	void ModelNode::compute () {
		// Update current model params with ours
		auto & modelParams = model_->getParameters ();
		for (auto i : index_range (this->dependencies ())) {
			auto v = DF::getValueUnsafe<double> (this->dependencies ()[i]);
			auto & p = modelParams[i];
			if (p.getValue () != v)
				model_->setParameterValue (p.getName (), v);
		}
	}

	std::string ModelNode::description () const { return "Model(" + model_->getName () + ")"; }

	void ComputeEquilibriumFrequenciesFromModelOp::compute (FrequencyVector & freqs,
	                                                        const SubstitutionModel * model) {
		auto & freqsFromModel = model->getFrequencies ();
		freqs = Eigen::Map<const FrequencyVector> (freqsFromModel.data (),
		                                           Eigen::Index (freqsFromModel.size ()));
	}

	void ComputeTransitionMatrixFromModelOp::compute (TransitionMatrix & matrix,
	                                                  const SubstitutionModel * model, double brlen) {
		auto & matrixFromModel = model->getPij_t (brlen);
		for (auto i : range (matrix.rows ()))
			for (auto j : range (matrix.cols ()))
				matrix (i, j) = matrixFromModel (std::size_t (i), std::size_t (j));
	}
}
}
