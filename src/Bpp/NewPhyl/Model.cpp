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

#include <Bpp/Exceptions.h>
#include <Bpp/NewPhyl/DataFlowInternalTemplates.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/NewPhyl/IntegerRange.h>
#include <Bpp/NewPhyl/LinearAlgebra.h>
#include <Bpp/NewPhyl/Model.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <algorithm>
#include <cassert>

namespace bpp {
namespace DF {
	/* Model values derivation:
	 *
	 * We only have analytical derivation for transition matrix with respect to brlen.
	 * A model is defined as "derivable(x)" if its parameters do not depend on "x".
	 *
	 * With a "derivable(x)" model:
	 * - equilibrium frequencies are a constant (derive to 0s).
	 * - transition matrices are only dependent on brlen (which may depend on "x", constant if not).
	 *
	 * derive(x) methods for model value compute nodes assume the model is derivable(x).
	 * derive(x) will not fail is not the case, but the derivative will be wrong.
	 * This is checked by an assert in debug mode.
	 * derive(x) for these nodes should not be called if isDerivable(x) is false.
	 * FIXME check anyway at derive(x), throw exception ?
	 *
	 * A non derivable(x) model is a model whose parameters depend on "x".
	 * It can be derived numerically.
	 */

	// Model DF Node

	Model::Model (NodeRefVec && deps, std::unique_ptr<SubstitutionModel> && model)
	    : Value<const SubstitutionModel *> (std::move (deps), model.get ()),
	      model_ (std::move (model)) {
		checkDependencyPattern (typeid (Model), dependencies (),
		                        ArrayOfValues<double>{nbParameters ()});
	}

	NodeRefVec createDependencyVector (const SubstitutionModel & model,
	                                   const std::map<std::string, ValueRef<double>> & depsByName) {
		NodeRefVec deps;
		const auto & modelParameters = model.getParameters ();
		for (auto i : range (modelParameters.size ())) {
			const auto & namespacedName = modelParameters[i].getName ();
			auto names = {namespacedName, model.getParameterNameWithoutNamespace (namespacedName)};

			auto findFirstFoundName = [&depsByName, &names]() {
				for (const auto & name : names) {
					auto it = depsByName.find (name);
					if (it != depsByName.end ())
						return it;
				}
				return depsByName.end ();
			};
			auto it = findFirstFoundName ();

			if (it != depsByName.end ()) {
				deps.emplace_back (it->second);
			} else {
				std::string msg =
				    "createDependencyVector(model, depsByName): parameter name not found. tested: ";
				for (const auto & name : names) {
					msg += name;
					msg += ", ";
				}
				throw Exception (std::move (msg));
			}
		}
		return deps;
	}

	Model::Model (std::unique_ptr<SubstitutionModel> model)
	    : Value<const SubstitutionModel *> (noDependency, model.get ()), model_ (std::move (model)) {
		// TODO remove
		const auto & parameters = model_->getParameters ();
		for (auto i : range (parameters.size ()))
			this->appendDependency (DF::makeNode<DF::Mutable<double>> (parameters[i].getValue ()));
	}

	Model::~Model () = default;

	SizeType Model::nbParameters () const noexcept { return this->dependencies ().size (); }

	// TODO remove too
	MutableRef<double> Model::getParameter (SizeType index) {
		assert (0 <= index);
		assert (index < this->nbDependencies ());
		return convertRef<DF::Mutable<double>> (this->dependency (index));
	}
	MutableRef<double> Model::getParameter (const std::string & name) {
		return getParameter (
		    static_cast<SizeType> (model_->getParameters ().whichParameterHasName (name)));
	}
	const std::string & Model::getParameterName (SizeType index) {
		return model_->getParameters ()[static_cast<std::size_t> (index)].getName ();
	}

	std::string Model::description () const { return "Model(" + model_->getName () + ")"; }
	std::string Model::debugInfo () const {
		return "nbState=" + std::to_string (model_->getAlphabet ()->getSize ());
	}

	bool Model::isDerivable (const Node & node) const {
		return std::none_of (
		    this->dependencies ().begin (), this->dependencies ().end (),
		    [&node](const NodeRef & dep) { return dep->isTransitivelyDependentOn (node); });
	}

	NodeRef Model::rebuild (NodeRefVec && deps) const {
		return makeNode<Model> (std::move (deps), std::unique_ptr<SubstitutionModel>{model_->clone ()});
	}

	void Model::compute () {
		// Update internal model bpp::Parameter with ours
		auto & modelParams = model_->getParameters ();
		for (auto i : range (nbParameters ())) {
			auto & v = accessValidValueConstCast<double> (this->dependency (i));
			auto & p = modelParams[static_cast<std::size_t> (i)];
			if (p.getValue () != v)
				model_->setParameterValue (model_->getParameterNameWithoutNamespace (p.getName ()), v);
		}
	}

	// Compute node functions

	class EquilibriumFrequenciesFromModel : public Value<VectorDouble> {
	public:
		using Dependencies = FunctionOfValues<const SubstitutionModel *>;
		EquilibriumFrequenciesFromModel (NodeRefVec && deps, SizeType nbStates)
		    : Value<VectorDouble> (std::move (deps), nbStates) {}
		std::string debugInfo () const override final {
			using std::to_string;
			return Value<VectorDouble>::debugInfo () + " nbState=" + to_string (dimensions (*this));
		}
		NodeRef derive (const Node & node) override final {
			assert (isDerivable (node));
			return Builder<Constant<VectorDouble>>::makeZero (dimensions (*this));
		}
		bool isDerivable (const Node & node) const override final {
			return derivableIfAllDepsAre (*this, node);
		}
		NodeRef rebuild (NodeRefVec && deps) const override final {
			return makeNode<EquilibriumFrequenciesFromModel> (std::move (deps), dimensions (*this).size);
		}

	private:
		void compute () override final {
			callWithValues (*this, [](VectorDouble & freqs, const SubstitutionModel * model) {
				auto & freqsFromModel = model->getFrequencies ();
				freqs = Eigen::Map<const VectorDouble> (freqsFromModel.data (),
				                                        static_cast<Eigen::Index> (freqsFromModel.size ()));
			});
		}
	};
	ValueRef<VectorDouble> Builder<EquilibriumFrequenciesFromModel>::make (NodeRefVec && deps,
	                                                                       SizeType nbStates) {
		checkDependencies<EquilibriumFrequenciesFromModel> (deps);
		return std::make_shared<EquilibriumFrequenciesFromModel> (std::move (deps), nbStates);
	}

	namespace {
		/* For now copy matrix cell by cell.
		 * TODO use eigen internally in SubsitutionModel ! (not perf critical for now though)
		 * FIXME if multithreading, internal model state must be removed !
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

	class TransitionMatrixFromModel : public Value<MatrixDouble> {
	public:
		using Dependencies = FunctionOfValues<const SubstitutionModel *, double>;

		TransitionMatrixFromModel (NodeRefVec && deps, const TransitionMatrixDimension & dim)
		    : Value<MatrixDouble> (std::move (deps), dim.rows, dim.cols) {}
		std::string debugInfo () const override final {
			return Value<MatrixDouble>::debugInfo () + " " +
			       to_string (TransitionMatrixDimension (dimensions (*this)));
		}
		NodeRef derive (const Node & node) override final {
			assert (isDerivable (node));
			auto dim = TransitionMatrixDimension (dimensions (*this));
			auto & modelNode = this->dependency (0);
			auto & brlenNode = this->dependency (1);
			auto dTransMat_dBrlen =
			    makeNode<TransitionMatrixFromModelBrlenDerivative> ({modelNode, brlenNode}, dim);
			return makeNode<CWiseMulScalarMatrixDouble> (
			    {brlenNode->derive (node), std::move (dTransMat_dBrlen)}, dim);
		}
		bool isDerivable (const Node & node) const override final {
			return derivableIfAllDepsAre (*this, node);
		}
		NodeRef rebuild (NodeRefVec && deps) const override final {
			return makeNode<TransitionMatrixFromModel> (std::move (deps), dimensions (*this));
		}

	private:
		void compute () override final {
			callWithValues (*this, [](MatrixDouble & matrix, const SubstitutionModel * model,
			                          double brlen) { bppToEigen (model->getPij_t (brlen), matrix); });
		}
	};
	ValueRef<MatrixDouble>
	Builder<TransitionMatrixFromModel>::make (NodeRefVec && deps,
	                                          const TransitionMatrixDimension & dim) {
		checkDependencies<TransitionMatrixFromModel> (deps);
		return std::make_shared<TransitionMatrixFromModel> (std::move (deps), dim);
	}

	class TransitionMatrixFromModelBrlenDerivative : public Value<MatrixDouble> {
	public:
		using Dependencies = FunctionOfValues<const SubstitutionModel *, double>;

		TransitionMatrixFromModelBrlenDerivative (NodeRefVec && deps,
		                                          const TransitionMatrixDimension & dim)
		    : Value<MatrixDouble> (std::move (deps), dim.rows, dim.cols) {}
		std::string debugInfo () const override final {
			return Value<MatrixDouble>::debugInfo () + " " +
			       to_string (TransitionMatrixDimension (dimensions (*this)));
		}
		NodeRef derive (const Node & node) override final {
			assert (isDerivable (node));
			auto dim = TransitionMatrixDimension (dimensions (*this));
			auto & modelNode = this->dependency (0);
			auto & brlenNode = this->dependency (1);
			auto d2TransMat_dBrlen2 =
			    makeNode<TransitionMatrixFromModelBrlenSecondDerivative> ({modelNode, brlenNode}, dim);
			return makeNode<CWiseMulScalarMatrixDouble> (
			    {brlenNode->derive (node), std::move (d2TransMat_dBrlen2)}, dim);
		}
		bool isDerivable (const Node & node) const override final {
			return derivableIfAllDepsAre (*this, node);
		}
		NodeRef rebuild (NodeRefVec && deps) const override final {
			return makeNode<TransitionMatrixFromModelBrlenDerivative> (std::move (deps),
			                                                           dimensions (*this));
		}

	private:
		void compute () override final {
			callWithValues (*this, [](MatrixDouble & matrix, const SubstitutionModel * model,
			                          double brlen) { bppToEigen (model->getdPij_dt (brlen), matrix); });
		}
	};
	ValueRef<MatrixDouble>
	Builder<TransitionMatrixFromModelBrlenDerivative>::make (NodeRefVec && deps,
	                                                         const TransitionMatrixDimension & dim) {
		checkDependencies<TransitionMatrixFromModelBrlenDerivative> (deps);
		return std::make_shared<TransitionMatrixFromModelBrlenDerivative> (std::move (deps), dim);
	}

	class TransitionMatrixFromModelBrlenSecondDerivative : public Value<MatrixDouble> {
	public:
		using Dependencies = FunctionOfValues<const SubstitutionModel *, double>;

		TransitionMatrixFromModelBrlenSecondDerivative (NodeRefVec && deps,
		                                                const TransitionMatrixDimension & dim)
		    : Value<MatrixDouble> (std::move (deps), dim.rows, dim.cols) {}
		std::string debugInfo () const override final {
			return Value<MatrixDouble>::debugInfo () + " " +
			       to_string (TransitionMatrixDimension (dimensions (*this)));
		}
		NodeRef rebuild (NodeRefVec && deps) const override final {
			return makeNode<TransitionMatrixFromModelBrlenSecondDerivative> (std::move (deps),
			                                                                 dimensions (*this));
		}

	private:
		void compute () override final {
			callWithValues (*this,
			                [](MatrixDouble & matrix, const SubstitutionModel * model, double brlen) {
				                bppToEigen (model->getd2Pij_dt2 (brlen), matrix);
			                });
		}
	};
	ValueRef<MatrixDouble> Builder<TransitionMatrixFromModelBrlenSecondDerivative>::make (
	    NodeRefVec && deps, const TransitionMatrixDimension & dim) {
		checkDependencies<TransitionMatrixFromModelBrlenSecondDerivative> (deps);
		return std::make_shared<TransitionMatrixFromModelBrlenSecondDerivative> (std::move (deps), dim);
	}
} // namespace DF
} // namespace bpp
