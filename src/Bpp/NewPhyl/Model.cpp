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

#include <Bpp/NewPhyl/DataFlowInternalTemplates.h>
#include <Bpp/NewPhyl/Model.h>
#include <Bpp/NewPhyl/Range.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <cassert>

namespace bpp {
namespace Phyl {
	namespace DF {
		// Model DF Node

		Model::Model (std::unique_ptr<SubstitutionModel> && model)
		    : Value<const SubstitutionModel *> (noDependency, model.get ()),
		      model_ (std::move (model)) {
			// TODO support already built paramater nodes
			model_->setNamespace ({}); // Delete namespace prefix
			const auto & parameters = model_->getParameters ();
			for (auto i : index_range (parameters))
				this->appendDependency (Parameter<double>::create (parameters[i].getValue ()));
		}

		Model::~Model () = default;

		SizeType Model::nbParameters () const noexcept { return this->dependencies ().size (); }
		ParameterRef<double> Model::getParameter (SizeType index) {
			return convertRef<bpp::DF::Parameter<double>> (this->dependencies ().at (index));
		}
		ParameterRef<double> Model::getParameter (const std::string & name) {
			return getParameter (
			    static_cast<SizeType> (model_->getParameters ().whichParameterHasName (name)));
		}
		const std::string & Model::getParameterName (SizeType index) {
			return model_->getParameters ()[static_cast<std::size_t> (index)].getName ();
		}

		void Model::compute () {
			// Update current model params with ours
			auto & modelParams = model_->getParameters ();
			for (auto i : index_range (this->dependencies ())) {
				auto v = accessValueUncheckedCast<double> (*this->dependencies ()[i]);
				auto & p = modelParams[static_cast<std::size_t> (i)];
				if (p.getValue () != v)
					model_->setParameterValue (p.getName (), v);
			}
		}

		std::string Model::description () const { return "Model(" + model_->getName () + ")"; }
		std::string Model::debugInfo () const {
			return "nbState=" + std::to_string (model_->getAlphabet ()->getSize ());
		}

		std::shared_ptr<Model> Model::create (std::unique_ptr<SubstitutionModel> && model) {
			return std::make_shared<Model> (std::move (model));
		}

		// Compute node functions

		EquilibriumFrequenciesFromModel::EquilibriumFrequenciesFromModel (NodeRefVec && deps,
		                                                                  SizeType nbStates)
		    : Value<VectorDouble> (std::move (deps), nbStates) {
			checkDependencies (*this);
		}
		void EquilibriumFrequenciesFromModel::compute () {
			callWithValues (*this, [](VectorDouble & freqs, const SubstitutionModel * model) {
				auto & freqsFromModel = model->getFrequencies ();
				freqs = Eigen::Map<const VectorDouble> (freqsFromModel.data (),
				                                        static_cast<Eigen::Index> (freqsFromModel.size ()));
			});
		}
		std::string EquilibriumFrequenciesFromModel::debugInfo () const {
			// TODO VectorDoubledebugInfo overload in Matrix.h
			auto & v = this->accessValue ();
			auto s = "nbState=" + std::to_string (v.rows ()) + " props=";
			if (v == VectorDouble::Zero (v.rows ()))
				s += '0';
			if (v == VectorDouble::Ones (v.rows ()))
				s += '1';
			return s;
		}
		std::shared_ptr<EquilibriumFrequenciesFromModel>
		EquilibriumFrequenciesFromModel::create (NodeRefVec && deps, SizeType nbStates) {
			return std::make_shared<EquilibriumFrequenciesFromModel> (std::move (deps), nbStates);
		}
		std::shared_ptr<EquilibriumFrequenciesFromModel>
		EquilibriumFrequenciesFromModel::create (ValueRef<const SubstitutionModel *> model,
		                                         SizeType nbStates) {
			return create (NodeRefVec{std::move (model)}, nbStates);
		}

		namespace {
			/* For now copy matrix cell by cell.
			 * TODO use eigen internally in SubsitutionModel ! (not perf critical for now though)
			 * FIXME if multithreading, internal model state may be a problem !
			 */
			void bppToEigen (const Matrix<double> & bppMatrix, MatrixDouble & eigenMatrix) {
				assert (eigenMatrix.rows () == bppMatrix.getNumberOfRows ());
				assert (eigenMatrix.cols () == bppMatrix.getNumberOfColumns ());
				for (auto i : range (eigenMatrix.rows ()))
					for (auto j : range (eigenMatrix.cols ()))
						eigenMatrix (i, j) =
						    bppMatrix (static_cast<std::size_t> (i), static_cast<std::size_t> (j));
			}
		} // namespace

		TransitionMatrixFromModel::TransitionMatrixFromModel (NodeRefVec && deps, SizeType nbStates)
		    : Value<MatrixDouble> (std::move (deps), nbStates, nbStates) {
			checkDependencies (*this);
		}
		void TransitionMatrixFromModel::compute () {
			callWithValues (*this, [](MatrixDouble & matrix, const SubstitutionModel * model,
			                          double brlen) { bppToEigen (model->getPij_t (brlen), matrix); });
		}
		std::string TransitionMatrixFromModel::debugInfo () const {
			return debugInfoFor (this->accessValue ());
		}
		std::shared_ptr<TransitionMatrixFromModel>
		TransitionMatrixFromModel::create (NodeRefVec && deps, SizeType nbStates) {
			return std::make_shared<TransitionMatrixFromModel> (std::move (deps), nbStates);
		}
		std::shared_ptr<TransitionMatrixFromModel>
		TransitionMatrixFromModel::create (ValueRef<const SubstitutionModel *> model,
		                                   ValueRef<double> brlen, SizeType nbStates) {
			return create (NodeRefVec{std::move (model), std::move (brlen)}, nbStates);
		}

		// TODO restore
		// bppToEigen (model->getdPij_dt (brlen), matrix);
		// bppToEigen (model->getd2Pij_dt2 (brlen), matrix);
	} // namespace DF
} // namespace Phyl
} // namespace bpp
