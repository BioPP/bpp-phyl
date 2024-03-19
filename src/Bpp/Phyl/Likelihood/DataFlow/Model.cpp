// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Exceptions.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Model.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Parametrizable.h>
#include <Bpp/Phyl/Model/MixedTransitionModel.h>


using namespace std;
using namespace bpp;

// Model node

ConfiguredModel::ConfiguredModel(
    Context& context,
    NodeRefVec&& deps,
    shared_ptr<BranchModelInterface>&& model) :
  Value<std::shared_ptr<BranchModelInterface>>(
    std::move(deps),
    model),
  AbstractParametrizable(model->getNamespace()), // , context_(context)
  model_(model)
{
  for (const auto& dep : dependencies())
  {
    shareParameter_(std::dynamic_pointer_cast<ConfiguredParameter>(dep));
  }
}

ConfiguredModel::~ConfiguredModel () = default;

std::string ConfiguredModel::description () const { return "Model(" + model_->getName () + ")"; }

std::string ConfiguredModel::debugInfo () const
{
  return "nbState=" + std::to_string (model_->getAlphabet ()->getSize ());
}

// Model node additional arguments = (type of BranchModel).
// Everything else is determined by the node dependencies.
bool ConfiguredModel::compareAdditionalArguments (const Node_DF& other) const
{
  const auto* derived = dynamic_cast<const Self*>(&other);
  if (derived == nullptr)
  {
    return false;
  }
  else
  {
    const auto& thisModel = *model_;
    const auto& otherModel = *derived->model_;
    return typeid (thisModel) == typeid (otherModel);
  }
}

std::size_t ConfiguredModel::hashAdditionalArguments () const
{
  const auto& bppModel = *model_;
  return typeid (bppModel).hash_code ();
}

NodeRef ConfiguredModel::recreate (Context& c, NodeRefVec&& deps)
{
  auto m = ConfiguredParametrizable::createConfigured<Target, Self>(c, std::move (deps), std::unique_ptr<Target>{dynamic_cast<Target*>(model_->clone ())});
  m->config = this->config; // Duplicate derivation config
  return m;
}

// EquilibriumFrequenciesFromModel

EquilibriumFrequenciesFromModel::EquilibriumFrequenciesFromModel (
    NodeRefVec&& deps, const Dimension<Eigen::RowVectorXd>& dim)
  : Value<Eigen::RowVectorXd>(std::move (deps)), targetDimension_ (dim)
{}

std::string EquilibriumFrequenciesFromModel::debugInfo () const
{
  using namespace numeric;
  return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
}

// EquilibriumFrequenciesFromModel additional arguments = ().
bool EquilibriumFrequenciesFromModel::compareAdditionalArguments (const Node_DF& other) const
{
  return dynamic_cast<const Self*>(&other) != nullptr;
}

NodeRef EquilibriumFrequenciesFromModel::derive (Context& c, const Node_DF& node)
{
  // d(equFreqs)/dn = sum_i d(equFreqs)/dx_i * dx_i/dn (x_i = model parameters)
  auto modelDep = this->dependency (0);
  auto& model = static_cast<Dep&>(*modelDep);
  NodeRef subNode = this->nbDependencies() == 1 ? 0 : this->dependency (1);
  auto buildFWithNewModel = [this, &c, &subNode](NodeRef&& newModel) {
        return ConfiguredParametrizable::createRowVector<Dep, Self>(c, {std::move (newModel), subNode}, targetDimension_);
      };

  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T>(
        c, model, node, targetDimension_, buildFWithNewModel);
  return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
}

NodeRef EquilibriumFrequenciesFromModel::recreate (Context& c, NodeRefVec&& deps)
{
  return ConfiguredParametrizable::createRowVector<Dep, Self>(c, std::move (deps), targetDimension_);
}

void EquilibriumFrequenciesFromModel::compute ()
{
  const Vdouble* freqsFromModel;
  const auto* mixmodel = dynamic_cast<const MixedTransitionModelInterface*>(accessValueConstCast<const BranchModelInterface*>(*this->dependency (0)));
  if (mixmodel && nbDependencies() > 1 && dependency(1))
  {
    auto nMod = accessValueConstCast<size_t>(*this->dependency(1));
    freqsFromModel = &mixmodel->nModel(nMod).getFrequencies();
  }
  else
  {
    const auto pmodel = dynamic_cast<const TransitionModelInterface*>(accessValueConstCast<const BranchModelInterface*>(*this->dependency (0)));
    if (pmodel)
      freqsFromModel = &pmodel->getFrequencies();
    else
      throw Exception("EquilibriumFrequenciesFromModel::compute only possible for Transition Models.");
  }

  auto& r = this->accessValueMutable ();
  r = Eigen::Map<const T>(freqsFromModel->data (), static_cast<Eigen::Index>(freqsFromModel->size ()));
}

// TransitionMatrixFromModel

TransitionMatrixFromModel::TransitionMatrixFromModel (NodeRefVec&& deps,
    const Dimension<Eigen::MatrixXd>& dim)
  : Value<Eigen::MatrixXd>(std::move (deps)), targetDimension_ (dim)
{}

std::string TransitionMatrixFromModel::debugInfo () const
{
  using namespace numeric;
  const auto nDeriv = accessValueConstCast<size_t>(*this->dependency (2));
  auto ret = debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_) + ":nDeriv=" + TextTools::toString(nDeriv);
  return ret;
}

// TransitionMatrixFromModel additional arguments = ().
bool TransitionMatrixFromModel::compareAdditionalArguments (const Node_DF& other) const
{
  return dynamic_cast<const Self*>(&other) != nullptr;
}

NodeRef TransitionMatrixFromModel::derive (Context& c, const Node_DF& node)
{
  // dtm/dn = sum_i dtm/dx_i * dx_i/dn + dtm/dbrlen * dbrlen/dn (x_i = model parameters).
  auto modelDep = this->dependency (0);
  auto brlenDep = this->dependency (1);
  NodeRef derivNode = this->dependency (2);
  NodeRef subNode = nbDependencies() < 4 ? 0 : this->dependency (3);

  const auto nDeriv = accessValueConstCast<size_t>(*derivNode);

  // Model part
  auto& model = static_cast<Dep&>(*modelDep);
  auto buildFWithNewModel = [this, &c, &brlenDep, &derivNode, &subNode](NodeRef&& newModel) {
        return ConfiguredParametrizable::createMatrix<Dep, Self>(c, {std::move (newModel), brlenDep, derivNode, subNode}, targetDimension_);
      };
  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T>(
        c, model, node, targetDimension_, buildFWithNewModel);
  // Brlen part, use specific node
  auto dbrlen_dn = brlenDep->derive (c, node);
  if (!dbrlen_dn->hasNumericalProperty (NumericalProperty::ConstantZero))
  {
    auto nDerivp = NumericConstant<size_t>::create(c, nDeriv + 1);

    auto df_dbrlen =
        ConfiguredParametrizable::createMatrix<Dep, TransitionMatrixFromModel>(c, {modelDep, brlenDep, nDerivp, subNode}, targetDimension_);
    derivativeSumDeps.emplace_back (CWiseMul<T, std::tuple<double, T>>::create (
          c, {std::move (dbrlen_dn), std::move (df_dbrlen)}, targetDimension_));
  }
  return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
}

NodeRef TransitionMatrixFromModel::recreate (Context& c, NodeRefVec&& deps)
{
  return ConfiguredParametrizable::createMatrix<Dep, Self>(c, std::move (deps), targetDimension_);
}

void TransitionMatrixFromModel::compute ()
{
  const auto brlen = accessValueConstCast<double>(*this->dependency (1)->dependency(0));
  const auto nDeriv = accessValueConstCast<size_t>(*this->dependency (2));

  auto& r = this->accessValueMutable ();

  const auto* model1 = accessValueConstCast<const BranchModelInterface*>(*this->dependency (0));

  const auto* mixmodel = dynamic_cast<const MixedTransitionModelInterface*>(accessValueConstCast<const BranchModelInterface*>(*this->dependency (0)));

#ifdef DEBUG
  std::cerr << "=== TransitionMatrixFromModel::compute === " << this << std::endl;
  std::cerr << "brlen= " << brlen << std::endl;
  std::cerr << "nDeriv=" << nDeriv << std::endl;
  std::cerr << "model= " << model1 << " : " << model1->getName() << std::endl;
  model1->getParameters().printParameters(std::cerr);
  std::cerr << "=== end TransitionMatrixFromModel::compute === " << this << std::endl << std::endl;
#endif

  if (mixmodel && nbDependencies() >= 4 && this->dependency(3)) // in case there is a submodel
  {
    auto nMod = accessValueConstCast<size_t>(*this->dependency (3));

    switch (nDeriv)
    {
    case 0:
      copyBppToEigen (mixmodel->nModel(nMod).getPij_t (brlen), r);
      break;
    case 1:
      copyBppToEigen (mixmodel->nModel(nMod).getdPij_dt (brlen), r);
      break;
    case 2:
      copyBppToEigen (mixmodel->nModel(nMod).getd2Pij_dt2 (brlen), r);
      break;
    default:
      throw Exception("TransitionMatrixFromModel likelihood derivate " + TextTools::toString(nDeriv) + " not defined.");
    }
  }
  else
  {
    const auto model = dynamic_cast<const TransitionModelInterface*>(model1);

    if (!model)
      throw Exception("TransitionMatrixFromModel::compute only possible for Transition Models.");

    switch (nDeriv)
    {
    case 0:
      copyBppToEigen (model->getPij_t (brlen), r);
      break;
    case 1:
      copyBppToEigen (model->getdPij_dt (brlen), r);
      break;
    case 2:
      copyBppToEigen (model->getd2Pij_dt2 (brlen), r);
      break;
    default:
      throw Exception("TransitionMatrixFromModel likelihood derivate " + TextTools::toString(nDeriv) + " not defined.");
    }
  }
}

////////////////////////////////////////////////////////////
// TransitionFunctionFromModel

TransitionFunctionFromModel::TransitionFunctionFromModel (NodeRefVec&& deps,
    const Dimension<T>& dim)
  : Value<TransitionFunction>(std::move (deps)), targetDimension_ (dim) {}

std::string TransitionFunctionFromModel::debugInfo () const
{
  using namespace numeric;
  const auto nDeriv = accessValueConstCast<size_t>(*this->dependency (2));
  auto ret = debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_) + ":nDeriv=" + TextTools::toString(nDeriv);
  return ret;
}

// TransitionFunctionFromModel additional arguments = ().
bool TransitionFunctionFromModel::compareAdditionalArguments (const Node_DF& other) const
{
  return dynamic_cast<const Self*>(&other) != nullptr;
}

std::shared_ptr<TransitionFunctionFromModel> TransitionFunctionFromModel::create (Context& c, NodeRefVec&& deps, const Dimension<T>& dim)
{
  checkDependenciesNotNull (typeid (Self), deps);
  checkDependencyVectorSize (typeid (Self), deps, 3);
  checkNthDependencyIs<ConfiguredModel>(typeid (Self), deps, 0);
  checkNthDependencyIs<ConfiguredParameter>(typeid (Self), deps, 1);
//  checkNthDependencyIs<size_t> (typeid (Self), deps, 2);

  return cachedAs<TransitionFunctionFromModel>(c, std::make_shared<TransitionFunctionFromModel>(std::move (deps), dim));
}

NodeRef TransitionFunctionFromModel::derive (Context& c, const Node_DF& node)
{
  // df(v)/dn = sum_i df(v)/dx_i * dx_i/dn + df(v)/dbrlen * dbrlen/dn + df(v)/dv * dv/dn (x_i = model parameters, v = variable of f).
  if (&node == &NodeX)
    throw Exception("TransitionFunctionFromModel::derive : Jacobian not implemented. Ask developpers.");

  auto modelDep = this->dependency (0);
  auto brlenDep = this->dependency (1);
  NodeRef derivNode = this->dependency (2);
  const auto nDeriv = accessValueConstCast<size_t>(*this->dependency (2));

  // Model part
  auto& model = static_cast<Dep&>(*modelDep);
  auto buildFWithNewModel = [this, &c, &brlenDep, &derivNode](NodeRef&& newModel) {
        return TransitionFunctionFromModel::create(c, {std::move (newModel), brlenDep, derivNode}, targetDimension_);
      };

  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T>(c, model, node, targetDimension_, buildFWithNewModel);

  // Brlen part, use specific node
  auto dbrlen_dn = brlenDep->derive (c, node);
  if (!dbrlen_dn->hasNumericalProperty (NumericalProperty::ConstantZero))
  {
    auto nDerivp = NumericConstant<size_t>::create(c, nDeriv + 1);
    auto df_dbrlen = TransitionFunctionFromModel::create(c, {modelDep, brlenDep, nDerivp}, targetDimension_);
    derivativeSumDeps.emplace_back (CWiseMul<T, std::tuple<double, T>>::create (
          c, {std::move (dbrlen_dn), std::move (df_dbrlen)}, targetDimension_));
  }

  // Vector part, use specific node


  return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), targetDimension_);
}

NodeRef TransitionFunctionFromModel::recreate (Context& c, NodeRefVec&& deps)
{
  return Self::create(c, std::move (deps), targetDimension_);
}

void TransitionFunctionFromModel::compute ()
{
  const auto brlen = accessValueConstCast<double>(*this->dependency (1)->dependency(0));
  const auto* model = accessValueConstCast<const BranchModelInterface*>(*this->dependency(0));
  const auto nDeriv = accessValueConstCast<size_t>(*this->dependency (2));

  auto& r = this->accessValueMutable ();

  auto dimin = VectorDimension(Eigen::Dynamic); // ttargetDimension_;
  auto dimout = Dimension<VectorLik>(Eigen::Dynamic, 1); // ttargetDimension_;

  r = [model, brlen, nDeriv, dimin, dimout](const VectorLik& values)
      {
        switch (nDeriv)
        {
        case 0:
          return numeric::convert(model->Lik_t(numeric::convert<Eigen::Dynamic, 1>(values, dimin), brlen), dimout);
        case 1:
          return numeric::convert(model->dLik_dt(numeric::convert<Eigen::Dynamic, 1>(values, dimin), brlen), dimout);
        case 2:
          return numeric::convert(model->d2Lik_dt2(numeric::convert<Eigen::Dynamic, 1>(values, dimin), brlen), dimout);
//          return numeric::convert<VectorLik, Eigen::VectorXd>(model->d2Lik_dt2(numeric::convert<Eigen::VectorXd, VectorLik>(values, dim), brlen), dim);
        default:
          throw Exception("TransitionFunctionFromModel likelihood derivate " + TextTools::toString(nDeriv) + " not defined.");
        }
      };
}


////////////////////////////////////////
// ProbabilitiesFromMixedModel

ProbabilitiesFromMixedModel::ProbabilitiesFromMixedModel (
    NodeRefVec&& deps, const Dimension<Eigen::RowVectorXd>& dim)
  : Value<Eigen::RowVectorXd>(std::move (deps)), nbClass_ (dim) {}


std::string ProbabilitiesFromMixedModel::debugInfo () const
{
  using namespace numeric;
  return debug (this->accessValueConst ()) + " nbClass=" + to_string (nbClass_);
}

// ProbabilitiesFromMixedModel additional arguments = ().
bool ProbabilitiesFromMixedModel::compareAdditionalArguments (const Node_DF& other) const
{
  return dynamic_cast<const Self*>(&other) != nullptr;
}

NodeRef ProbabilitiesFromMixedModel::derive (Context& c, const Node_DF& node)
{
  // d(Prob)/dn = sum_i d(Prob)/dx_i * dx_i/dn (x_i = distrib parameters)
  auto modelDep = this->dependency (0);
  auto& model = static_cast<ConfiguredModel&>(*modelDep);
  auto buildPWithNewModel = [this, &c](NodeRef&& newModel) {
        return ConfiguredParametrizable::createRowVector<Dep, Self>(c, {std::move (newModel)}, nbClass_);
      };

  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T >(
        c, model, node, nbClass_, buildPWithNewModel);
  return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), nbClass_);
}

NodeRef ProbabilitiesFromMixedModel::recreate (Context& c, NodeRefVec&& deps)
{
  return ConfiguredParametrizable::createRowVector<Dep, Self>(c, std::move (deps), nbClass_);
}

void ProbabilitiesFromMixedModel::compute ()
{
  const auto* mixmodel = dynamic_cast<const MixedTransitionModelInterface*>(accessValueConstCast<const BranchModelInterface*>(*this->dependency (0)));
  const auto& probasFromModel = mixmodel->getProbabilities ();
  auto& r = this->accessValueMutable ();
  r = Eigen::Map<const T>(probasFromModel.data(), static_cast<Eigen::Index>(probasFromModel.size ()));
}

std::shared_ptr<ProbabilitiesFromMixedModel> ProbabilitiesFromMixedModel::create(Context& c, NodeRefVec&& deps)
{
  checkDependenciesNotNull (typeid (Self), deps);
  checkDependencyVectorSize (typeid (Self), deps, 1);
  checkNthDependencyIs<ConfiguredModel>(typeid (Self), deps, 0);
  auto& model = static_cast<ConfiguredModel&>(*deps[0]);
  auto& deps0 = *deps[0];
  const auto mixmodel = dynamic_pointer_cast<const MixedTransitionModelInterface>(model.targetValue());
  if (!mixmodel)
    failureDependencyTypeMismatch(typeid(Self), 0, typeid(MixedTransitionModelInterface), typeid(deps0));

  size_t nbCat = mixmodel->getNumberOfModels();
  return cachedAs<ProbabilitiesFromMixedModel>(c, make_shared<ProbabilitiesFromMixedModel>(std::move(deps), RowVectorDimension(Eigen::Index(nbCat))));
}


////////////////////////////////////////////////////
// ProbabilityFromMixedModel

ProbabilityFromMixedModel::ProbabilityFromMixedModel (
    NodeRefVec&& deps, size_t nCat)
  : Value<double>(std::move (deps)), nCat_ (nCat) {}


std::string ProbabilityFromMixedModel::debugInfo () const
{
  using namespace numeric;
  return "proba=" + TextTools::toString(accessValueConst()) + ":nCat=" + TextTools::toString(nCat_);
}

// ProbabilityFromMixedModel additional arguments = ().
bool ProbabilityFromMixedModel::compareAdditionalArguments (const Node_DF& other) const
{
  const auto* derived = dynamic_cast<const Self*>(&other);
  return derived != nullptr && nCat_ == derived->nCat_;
}

std::shared_ptr<ProbabilityFromMixedModel> ProbabilityFromMixedModel::create(Context& c, NodeRefVec&& deps, size_t nCat)
{
  checkDependenciesNotNull (typeid (Self), deps);
  checkDependencyVectorSize (typeid (Self), deps, 1);
  checkNthDependencyIs<ConfiguredModel>(typeid (Self), deps, 0);
  auto& model = static_cast<ConfiguredModel&>(*deps[0]);
  const auto& deps0 = *deps[0];
  const auto mixmodel = dynamic_pointer_cast<const MixedTransitionModelInterface>(model.targetValue());
  if (!mixmodel)
    failureDependencyTypeMismatch (typeid(Self), 0, typeid(MixedTransitionModelInterface), typeid(deps0));

  return cachedAs<ProbabilityFromMixedModel>(c, std::make_shared<ProbabilityFromMixedModel>(std::move (deps), nCat));
}

NodeRef ProbabilityFromMixedModel::derive (Context& c, const Node_DF& node)
{
  // d(Prob)/dn = sum_i d(Prob)/dx_i * dx_i/dn (x_i = distrib parameters)
  auto modelDep = this->dependency (0);
  auto& model = static_cast<Dep&>(*modelDep);
  auto buildPWithNewModel = [this, &c](NodeRef&& newModel) {
        return this->create (c, {std::move (newModel)}, nCat_);
      };

  NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T >(
        c, model, node, 1, buildPWithNewModel);
  return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), 1);
}

NodeRef ProbabilityFromMixedModel::recreate (Context& c, NodeRefVec&& deps)
{
  return ProbabilityFromMixedModel::create (c, {std::move (deps)}, nCat_);
}

void ProbabilityFromMixedModel::compute ()
{
  const auto* mixmodel = dynamic_cast<const MixedTransitionModelInterface*>(accessValueConstCast<const BranchModelInterface*>(*this->dependency (0)));
  this->accessValueMutable () = mixmodel->getNProbability(nCat_);
}
