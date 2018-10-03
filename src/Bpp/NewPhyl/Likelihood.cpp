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
    /* For now copy matrix cell by cell.
     * TODO use eigen internally in SubstitutionModel ! (not perf critical for now though)
     * FIXME if multithreading, internal model state must be removed !
     */
    static void copyBppToEigen (const Matrix<double> & bppMatrix, Eigen::MatrixXd & eigenMatrix) {
      const auto eigenRows = static_cast<Eigen::Index> (bppMatrix.getNumberOfRows ());
      const auto eigenCols = static_cast<Eigen::Index> (bppMatrix.getNumberOfColumns ());
      eigenMatrix.resize (eigenRows, eigenCols);
      for (Eigen::Index i = 0; i < eigenRows; ++i) {
        for (Eigen::Index j = 0; j < eigenCols; ++j) {
          eigenMatrix (i, j) = bppMatrix (static_cast<std::size_t> (i), static_cast<std::size_t> (j));
        }
      }
    }

    /* Helper function for generating numerical derivatives of model computation nodes.
     * Assuming we have a v = f(model, stuff) node, with v of type T.
     * df/dn = sum_i df/dx_i * dx_i/dn + df/dstuff * dstuff/dn.
     * This function returns a NodeRefVec containings the nodes: {df/dx_i * dx_i/dn} for i in order.
     * buildFWithNewModel(newModel) should create the f(newModel, stuff) node.
     */
    template <typename T, typename B>
    static NodeRefVec generateModelDerivativeSumDepsForModelComputations (
      Context & c, ConfiguredModel & model, const Node & derivationNode, const Dimension<T> & targetDimension,
      B buildFWithNewModel) {
      NodeRefVec derivativeSumDeps;
      for (std::size_t i = 0; i < model.nbDependencies (); ++i) {
        // First compute dxi_dn. If this maps to a constant 0, do not compute df_dxi at all (costly).
        auto dxi_dn = model.dependency (i)->derive (c, derivationNode);
        if (!dxi_dn->hasNumericalProperty (NumericalProperty::ConstantZero)) {
          auto buildFWithNewXi = [&c, i, &model, &buildFWithNewModel](ValueRef<double> newDep) {
            // The sub-graph that will be replicated with shifted inputs is: f(model(x_i), stuff)
            NodeRefVec newModelDeps = model.dependencies ();
            newModelDeps[i] = std::move (newDep);
            auto newModel = model.recreate (c, std::move (newModelDeps));
            return buildFWithNewModel (std::move (newModel));
          };
          auto df_dxi = generateNumericalDerivative<T, double> (
            c, model.config, model.dependency (i), Dimension<double>{}, targetDimension, buildFWithNewXi);
          derivativeSumDeps.emplace_back (CWiseMul<T, std::tuple<double, T>>::create (
            c, {std::move (dxi_dn), std::move (df_dxi)}, targetDimension));
        }
      }
      return derivativeSumDeps;
    }

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

    NodeRefVec createDependencyVector (const TransitionModel & model,
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

    std::shared_ptr<ConfiguredModel> ConfiguredModel::create (Context & c, NodeRefVec && deps,
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
      : Value<const TransitionModel *> (std::move (deps), model.get ()), model_ (std::move (model)) {}

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

    NodeRef ConfiguredModel::recreate (Context & c, NodeRefVec && deps) {
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
      : Value<Eigen::RowVectorXd> (std::move (deps)), targetDimension_ (dim) {}

    std::string EquilibriumFrequenciesFromModel::debugInfo () const {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    NodeRef EquilibriumFrequenciesFromModel::derive (Context & c, const Node & node) {
      // d(equFreqs)/dn = sum_i d(equFreqs)/dx_i * dx_i/dn (x_i = model parameters)
      auto modelDep = this->dependency (0);
      auto & model = static_cast<ConfiguredModel &> (*modelDep);
      auto buildFWithNewModel = [this, &c](NodeRef && newModel) {
        return Self::create (c, {std::move (newModel)}, targetDimension_);
      };
      NodeRefVec derivativeSumDeps = generateModelDerivativeSumDepsForModelComputations<T> (
        c, model, node, targetDimension_, buildFWithNewModel);
      return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
    }

    NodeRef EquilibriumFrequenciesFromModel::recreate (Context & c, NodeRefVec && deps) {
      return Self::create (c, std::move (deps), targetDimension_);
    }

    void EquilibriumFrequenciesFromModel::compute () {
      const auto * model = accessValueConstCast<const TransitionModel *> (*this->dependency (0));
      const auto & freqsFromModel = model->getFrequencies ();
      auto & r = this->accessValueMutable ();
      r = Eigen::Map<const T> (freqsFromModel.data (), static_cast<Eigen::Index> (freqsFromModel.size ()));
    }

    // TransitionMatrixFromModel

    ValueRef<Eigen::MatrixXd> TransitionMatrixFromModel::create (Context & c, NodeRefVec && deps,
                                                                 const Dimension<Eigen::MatrixXd> & dim) {
      checkDependenciesNotNull (typeid (Self), deps);
      checkDependencyVectorSize (typeid (Self), deps, 2);
      checkNthDependencyIs<ConfiguredModel> (typeid (Self), deps, 0);
      checkNthDependencyIsValue<double> (typeid (Self), deps, 1);
      return std::make_shared<Self> (std::move (deps), dim);
    }

    TransitionMatrixFromModel::TransitionMatrixFromModel (NodeRefVec && deps,
                                                          const Dimension<Eigen::MatrixXd> & dim)
      : Value<Eigen::MatrixXd> (std::move (deps)), targetDimension_ (dim) {}

    std::string TransitionMatrixFromModel::debugInfo () const {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    NodeRef TransitionMatrixFromModel::derive (Context & c, const Node & node) {
      // dtm/dn = sum_i dtm/dx_i * dx_i/dn + dtm/dbrlen + dbrlen/dn (x_i = model parameters).
      auto modelDep = this->dependency (0);
      auto brlenDep = this->dependency (1);
      // Model part
      auto & model = static_cast<ConfiguredModel &> (*modelDep);
      auto buildFWithNewModel = [this, &c, &brlenDep](NodeRef && newModel) {
        return Self::create (c, {std::move (newModel), brlenDep}, targetDimension_);
      };
      NodeRefVec derivativeSumDeps = generateModelDerivativeSumDepsForModelComputations<T> (
        c, model, node, targetDimension_, buildFWithNewModel);
      // Brlen part, use specific node
      auto dbrlen_dn = brlenDep->derive (c, node);
      if (!dbrlen_dn->hasNumericalProperty (NumericalProperty::ConstantZero)) {
        auto df_dbrlen =
          TransitionMatrixFromModelFirstBrlenDerivative::create (c, {modelDep, brlenDep}, targetDimension_);
        derivativeSumDeps.emplace_back (CWiseMul<T, std::tuple<double, T>>::create (
          c, {std::move (dbrlen_dn), std::move (df_dbrlen)}, targetDimension_));
      }
      return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
    }

    NodeRef TransitionMatrixFromModel::recreate (Context & c, NodeRefVec && deps) {
      return Self::create (c, std::move (deps), targetDimension_);
    }

    void TransitionMatrixFromModel::compute () {
      const auto * model = accessValueConstCast<const TransitionModel *> (*this->dependency (0));
      const auto brlen = accessValueConstCast<double> (*this->dependency (1));
      auto & r = this->accessValueMutable ();
      copyBppToEigen (model->getPij_t (brlen), r);
    }

    // TransitionMatrixFromModelFirstBrlenDerivative

    ValueRef<Eigen::MatrixXd>
    TransitionMatrixFromModelFirstBrlenDerivative::create (Context & c, NodeRefVec && deps,
                                                           const Dimension<Eigen::MatrixXd> & dim) {
      checkDependenciesNotNull (typeid (Self), deps);
      checkDependencyVectorSize (typeid (Self), deps, 2);
      checkNthDependencyIs<ConfiguredModel> (typeid (Self), deps, 0);
      checkNthDependencyIsValue<double> (typeid (Self), deps, 1);
      return std::make_shared<Self> (std::move (deps), dim);
    }

    TransitionMatrixFromModelFirstBrlenDerivative::TransitionMatrixFromModelFirstBrlenDerivative (
      NodeRefVec && deps, const Dimension<Eigen::MatrixXd> & dim)
      : Value<Eigen::MatrixXd> (std::move (deps)), targetDimension_ (dim) {}

    std::string TransitionMatrixFromModelFirstBrlenDerivative::debugInfo () const {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    NodeRef TransitionMatrixFromModelFirstBrlenDerivative::derive (Context & c, const Node & node) {
      // dtm/dn = sum_i dtm/dx_i * dx_i/dn + dtm/dbrlen + dbrlen/dn (x_i = model parameters).
      auto modelDep = this->dependency (0);
      auto brlenDep = this->dependency (1);
      // Model part
      auto & model = static_cast<ConfiguredModel &> (*modelDep);
      auto buildFWithNewModel = [this, &c, &brlenDep](NodeRef && newModel) {
        return Self::create (c, {std::move (newModel), brlenDep}, targetDimension_);
      };
      NodeRefVec derivativeSumDeps = generateModelDerivativeSumDepsForModelComputations<T> (
        c, model, node, targetDimension_, buildFWithNewModel);
      // Brlen part, use specific node
      auto dbrlen_dn = brlenDep->derive (c, node);
      if (!dbrlen_dn->hasNumericalProperty (NumericalProperty::ConstantZero)) {
        auto df_dbrlen =
          TransitionMatrixFromModelSecondBrlenDerivative::create (c, {modelDep, brlenDep}, targetDimension_);
        derivativeSumDeps.emplace_back (CWiseMul<T, std::tuple<double, T>>::create (
          c, {std::move (dbrlen_dn), std::move (df_dbrlen)}, targetDimension_));
      }
      return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
    }

    NodeRef TransitionMatrixFromModelFirstBrlenDerivative::recreate (Context & c, NodeRefVec && deps) {
      return Self::create (c, std::move (deps), targetDimension_);
    }

    void TransitionMatrixFromModelFirstBrlenDerivative::compute () {
      const auto * model = accessValueConstCast<const TransitionModel *> (*this->dependency (0));
      const auto brlen = accessValueConstCast<double> (*this->dependency (1));
      auto & r = this->accessValueMutable ();
      copyBppToEigen (model->getdPij_dt (brlen), r);
    }

    // TransitionMatrixFromModelSecondBrlenDerivative

    ValueRef<Eigen::MatrixXd>
    TransitionMatrixFromModelSecondBrlenDerivative::create (Context & c, NodeRefVec && deps,
                                                            const Dimension<Eigen::MatrixXd> & dim) {
      checkDependenciesNotNull (typeid (Self), deps);
      checkDependencyVectorSize (typeid (Self), deps, 2);
      checkNthDependencyIs<ConfiguredModel> (typeid (Self), deps, 0);
      checkNthDependencyIsValue<double> (typeid (Self), deps, 1);
      return std::make_shared<Self> (std::move (deps), dim);
    }

    TransitionMatrixFromModelSecondBrlenDerivative::TransitionMatrixFromModelSecondBrlenDerivative (
      NodeRefVec && deps, const Dimension<Eigen::MatrixXd> & dim)
      : Value<Eigen::MatrixXd> (std::move (deps)), targetDimension_ (dim) {}

    std::string TransitionMatrixFromModelSecondBrlenDerivative::debugInfo () const {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    NodeRef TransitionMatrixFromModelSecondBrlenDerivative::derive (Context & c, const Node & node) {
      // dtm/dn = sum_i dtm/dx_i * dx_i/dn + dtm/dbrlen + dbrlen/dn (x_i = model parameters).
      auto modelDep = this->dependency (0);
      auto brlenDep = this->dependency (1);
      // Model part
      auto & model = static_cast<ConfiguredModel &> (*modelDep);
      auto buildFWithNewModel = [this, &c, &brlenDep](NodeRef && newModel) {
        return Self::create (c, {std::move (newModel), brlenDep}, targetDimension_);
      };
      NodeRefVec derivativeSumDeps = generateModelDerivativeSumDepsForModelComputations<T> (
        c, model, node, targetDimension_, buildFWithNewModel);
      // Brlen part : no specific node, use numerical derivation.
      auto dbrlen_dn = brlenDep->derive (c, node);
      if (!dbrlen_dn->hasNumericalProperty (NumericalProperty::ConstantZero)) {
        auto buildFWithNewBrlen = [this, &c, &modelDep](ValueRef<double> newBrlen) {
          return Self::create (c, {modelDep, std::move (newBrlen)}, targetDimension_);
        };
        auto df_dbrlen = generateNumericalDerivative<T, double> (
          c, model.config, brlenDep, Dimension<double> (), targetDimension_, buildFWithNewBrlen);
        derivativeSumDeps.emplace_back (CWiseMul<T, std::tuple<double, T>>::create (
          c, {std::move (dbrlen_dn), std::move (df_dbrlen)}, targetDimension_));
      }
      return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
    }

    NodeRef TransitionMatrixFromModelSecondBrlenDerivative::recreate (Context & c, NodeRefVec && deps) {
      return Self::create (c, std::move (deps), targetDimension_);
    }

    void TransitionMatrixFromModelSecondBrlenDerivative::compute () {
      const auto * model = accessValueConstCast<const TransitionModel *> (*this->dependency (0));
      const auto brlen = accessValueConstCast<double> (*this->dependency (1));
      auto & r = this->accessValueMutable ();
      copyBppToEigen (model->getd2Pij_dt2 (brlen), r);
    }
  } // namespace dataflow
} // namespace bpp
