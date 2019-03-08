//
// File: Model.cpp
// Authors:
// Created: mercredi 10 octobre 2018, à 15h 34
// Last modified: mercredi 10 octobre 2018, à 15h 34
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
#include <Bpp/NewPhyl/Model.h>
#include <Bpp/NewPhyl/Parametrizable.h>

using namespace std;

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

    // Model node

    ConfiguredModel::ConfiguredModel (Context& context, NodeRefVec && deps, std::unique_ptr<TransitionModel> && model)
      : Value<const TransitionModel*> (std::move (deps), model.get ()), AbstractParametrizable(model->getNamespace()), context_(context), model_(std::move(model))
    {
      for (const auto& dep:dependencies())
      {
        const auto& param=std::dynamic_pointer_cast<ConfiguredParameter>(dep);
        shareParameter_(param);
      }
    }

    ConfiguredModel::~ConfiguredModel () = default;

    std::string ConfiguredModel::description () const { return "Model(" + model_->getName () + ")"; }

    std::string ConfiguredModel::debugInfo () const {
      return "nbState=" + std::to_string (model_->getAlphabet ()->getSize ());
    }

    // Model node additional arguments = (type of bpp::TransitionModel).
    // Everything else is determined by the node dependencies.
    bool ConfiguredModel::compareAdditionalArguments (const Node & other) const {
      const auto * derived = dynamic_cast<const Self *> (&other);
      if (derived == nullptr) {
        return false;
      } else {
        const auto & thisModel = *model_;
        const auto & otherModel = *derived->model_;
        return typeid (thisModel) == typeid (otherModel);
      }
    }
    
    std::size_t ConfiguredModel::hashAdditionalArguments () const {
      const auto & bppModel = *model_;
      return typeid (bppModel).hash_code ();
    }

    NodeRef ConfiguredModel::recreate (Context & c, NodeRefVec && deps) {
      auto m = ConfiguredParametrizable::createConfigured<Target, Self> (c, std::move (deps), std::unique_ptr<Target>{model_->clone ()});
      m->config = this->config; // Duplicate derivation config
      return m;
    }
    
    // EquilibriumFrequenciesFromModel

    EquilibriumFrequenciesFromModel::EquilibriumFrequenciesFromModel (
      NodeRefVec && deps, const Dimension<Eigen::RowVectorXd> & dim)
      : Value<Eigen::RowVectorXd> (std::move (deps)), targetDimension_ (dim) {}

    std::string EquilibriumFrequenciesFromModel::debugInfo () const {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    // EquilibriumFrequenciesFromModel additional arguments = ().
    bool EquilibriumFrequenciesFromModel::compareAdditionalArguments (const Node & other) const {
      return dynamic_cast<const Self *> (&other) != nullptr;
    }

    NodeRef EquilibriumFrequenciesFromModel::derive (Context & c, const Node & node) {
      // d(equFreqs)/dn = sum_i d(equFreqs)/dx_i * dx_i/dn (x_i = model parameters)
      auto modelDep = this->dependency (0);
      auto & model = static_cast<Dep &> (*modelDep);
      auto buildFWithNewModel = [this, &c](NodeRef && newModel) {
        return ConfiguredParametrizable::createVector<Dep, Self> (c, {std::move (newModel)}, targetDimension_);
      };
      
      NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T> (
        c, model, node, targetDimension_, buildFWithNewModel);
      return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
    }

    NodeRef EquilibriumFrequenciesFromModel::recreate (Context & c, NodeRefVec && deps) {
      return ConfiguredParametrizable::createVector<Dep, Self> (c, std::move (deps), targetDimension_);
    }

    void EquilibriumFrequenciesFromModel::compute () {
      const auto * model = accessValueConstCast<const TransitionModel *> (*this->dependency (0));
      const auto & freqsFromModel = model->getFrequencies ();
      auto & r = this->accessValueMutable ();
      r = Eigen::Map<const T> (freqsFromModel.data (), static_cast<Eigen::Index> (freqsFromModel.size ()));
    }

    // TransitionMatrixFromModel

    TransitionMatrixFromModel::TransitionMatrixFromModel (NodeRefVec && deps,
                                                          const Dimension<Eigen::MatrixXd> & dim)
      : Value<Eigen::MatrixXd> (std::move (deps)), targetDimension_ (dim) {}

    std::string TransitionMatrixFromModel::debugInfo () const {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    // TransitionMatrixFromModel additional arguments = ().
    bool TransitionMatrixFromModel::compareAdditionalArguments (const Node & other) const {
      return dynamic_cast<const Self *> (&other) != nullptr;
    }

    NodeRef TransitionMatrixFromModel::derive (Context & c, const Node & node) {
      // dtm/dn = sum_i dtm/dx_i * dx_i/dn + dtm/dbrlen * dbrlen/dn (x_i = model parameters).
      
      auto modelDep = this->dependency (0);
      auto brlenDep = this->dependency (1);
      // Model part
      auto & model = static_cast<Dep &> (*modelDep);
      auto buildFWithNewModel = [this, &c, &brlenDep](NodeRef && newModel) {
        return ConfiguredParametrizable::createMatrix<Dep, Self> (c, {std::move (newModel), brlenDep}, targetDimension_);
      };
      NodeRefVec derivativeSumDeps = bpp::dataflow::ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T> (
        c, model, node, targetDimension_, buildFWithNewModel);
      // Brlen part, use specific node
      auto dbrlen_dn = brlenDep->derive (c, node);
      if (!dbrlen_dn->hasNumericalProperty (NumericalProperty::ConstantZero)) {
        auto df_dbrlen =
          ConfiguredParametrizable::createMatrix<Dep, TransitionMatrixFromModelFirstBrlenDerivative>(c, {modelDep, brlenDep}, targetDimension_);
        derivativeSumDeps.emplace_back (CWiseMul<T, std::tuple<double, T>>::create (
                                          c, {std::move (dbrlen_dn), std::move (df_dbrlen)}, targetDimension_));
      }
      return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
    }

    NodeRef TransitionMatrixFromModel::recreate (Context & c, NodeRefVec && deps) {
      return ConfiguredParametrizable::createMatrix<Dep, Self> (c, std::move (deps), targetDimension_);
    }

    void TransitionMatrixFromModel::compute () {
      const auto * model = accessValueConstCast<const TransitionModel *> (*this->dependency (0));
      const auto brlen = accessValueConstCast<double> (*this->dependency (1)->dependency(0));
      auto & r = this->accessValueMutable ();
      copyBppToEigen (model->getPij_t (brlen), r);
    }

    // TransitionMatrixFromModelFirstBrlenDerivative

    TransitionMatrixFromModelFirstBrlenDerivative::TransitionMatrixFromModelFirstBrlenDerivative (
      NodeRefVec && deps, const Dimension<Eigen::MatrixXd> & dim)
      : Value<Eigen::MatrixXd> (std::move (deps)), targetDimension_ (dim) {}

    std::string TransitionMatrixFromModelFirstBrlenDerivative::debugInfo () const {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    // TransitionMatrixFromModelFirstBrlenDerivative additional arguments = ().
    bool
    TransitionMatrixFromModelFirstBrlenDerivative::compareAdditionalArguments (const Node & other) const {
      return dynamic_cast<const Self *> (&other) != nullptr;
    }

    NodeRef TransitionMatrixFromModelFirstBrlenDerivative::derive (Context & c, const Node & node) {
      // dtm/dn = sum_i dtm/dx_i * dx_i/dn + dtm/dbrlen + dbrlen/dn (x_i = model parameters).
      auto modelDep = this->dependency (0);
      auto brlenDep = this->dependency (1);
      // Model part
      auto & model = static_cast<Dep &> (*modelDep);
      auto buildFWithNewModel = [this, &c, &brlenDep](NodeRef && newModel) {
        return ConfiguredParametrizable::createMatrix<Dep, Self> (c, {std::move (newModel), brlenDep}, targetDimension_);
      };
      
      NodeRefVec derivativeSumDeps = bpp::dataflow::ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T> (
        c, model, node, targetDimension_, buildFWithNewModel);
      // Brlen part, use specific node
      auto dbrlen_dn = brlenDep->derive (c, node);
      if (!dbrlen_dn->hasNumericalProperty (NumericalProperty::ConstantZero)) {
        auto df_dbrlen =
          ConfiguredParametrizable::createMatrix<Dep, TransitionMatrixFromModelSecondBrlenDerivative> (c, {modelDep, brlenDep}, targetDimension_);
        derivativeSumDeps.emplace_back (CWiseMul<T, std::tuple<double, T>>::create (
                                          c, {std::move (dbrlen_dn), std::move (df_dbrlen)}, targetDimension_));
      }
      return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
    }

    NodeRef TransitionMatrixFromModelFirstBrlenDerivative::recreate (Context & c, NodeRefVec && deps) {
      return ConfiguredParametrizable::createMatrix<Dep, Self> (c, std::move (deps), targetDimension_);
    }

    void TransitionMatrixFromModelFirstBrlenDerivative::compute () {
      const auto * model = accessValueConstCast<const TransitionModel *> (*this->dependency (0));
      const auto brlen = accessValueConstCast<double> (*this->dependency (1)->dependency(0));
      auto & r = this->accessValueMutable ();
      copyBppToEigen (model->getdPij_dt (brlen), r);
    }

    
    ////////////////////////////////////////////////////////
    // TransitionMatrixFromModelSecondBrlenDerivative

    TransitionMatrixFromModelSecondBrlenDerivative::TransitionMatrixFromModelSecondBrlenDerivative (
      NodeRefVec && deps, const Dimension<Eigen::MatrixXd> & dim)
      : Value<Eigen::MatrixXd> (std::move (deps)), targetDimension_ (dim) {}

    std::string TransitionMatrixFromModelSecondBrlenDerivative::debugInfo () const {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    // TransitionMatrixFromModelSecondBrlenDerivative additional arguments = ().
    bool
    TransitionMatrixFromModelSecondBrlenDerivative::compareAdditionalArguments (const Node & other) const {
      return dynamic_cast<const Self *> (&other) != nullptr;
    }

    NodeRef TransitionMatrixFromModelSecondBrlenDerivative::derive (Context & c, const Node & node) {
      // dtm/dn = sum_i dtm/dx_i * dx_i/dn + dtm/dbrlen + dbrlen/dn (x_i = model parameters).
      auto modelDep = this->dependency (0);
      auto brlenDep = this->dependency (1);
      // Model part
      auto & model = static_cast<Dep &> (*modelDep);
      auto buildFWithNewModel = [this, &c, &brlenDep](NodeRef && newModel) {
        return ConfiguredParametrizable::createMatrix<Dep, Self> (c, {std::move (newModel), brlenDep}, targetDimension_);
      };
      NodeRefVec derivativeSumDeps = bpp::dataflow::ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T> (
        c, model, node, targetDimension_, buildFWithNewModel);
      // Brlen part : no specific node, use numerical derivation.
      auto dbrlen_dn = brlenDep->derive (c, node);
      if (!dbrlen_dn->hasNumericalProperty (NumericalProperty::ConstantZero)) {
        auto buildFWithNewBrlen = [this, &c, &modelDep](ValueRef<double> newBrlen) {
          return ConfiguredParametrizable::createMatrix<Dep, Self> (c, {modelDep, std::move (newBrlen)}, targetDimension_);
        };
        auto df_dbrlen = generateNumericalDerivative<T, double> (
          c, model.config, brlenDep, Dimension<double> (), targetDimension_, buildFWithNewBrlen);
        derivativeSumDeps.emplace_back (CWiseMul<T, std::tuple<double, T>>::create (
                                          c, {std::move (dbrlen_dn), std::move (df_dbrlen)}, targetDimension_));
      }
      return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
    }

    NodeRef TransitionMatrixFromModelSecondBrlenDerivative::recreate (Context & c, NodeRefVec && deps) {
      return ConfiguredParametrizable::createMatrix<Dep, Self> (c, std::move (deps), targetDimension_);
    }

    void TransitionMatrixFromModelSecondBrlenDerivative::compute () {
      const auto * model = accessValueConstCast<const TransitionModel *> (*this->dependency (0));
      const auto brlen = accessValueConstCast<double> (*this->dependency (1)->dependency(0));
      auto & r = this->accessValueMutable ();
      copyBppToEigen (model->getd2Pij_dt2 (brlen), r);
    }
    
  } // namespace dataflow
} // namespace bpp
