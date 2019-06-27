//
// File: DiscreteDistribution.cpp
// Authors:
// Created: jeudi 25 octobre 2018, à 17h 23
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
#include <Bpp/NewPhyl/DiscreteDistribution.h>

using namespace std;

namespace bpp {
  namespace dataflow {

    ConfiguredDistribution::ConfiguredDistribution (Context& context, NodeRefVec && deps, std::unique_ptr<DiscreteDistribution> && distrib)
      : Value<const DiscreteDistribution*> (std::move (deps), distrib.get ()), AbstractParametrizable(distrib->getNamespace()), context_(context), distrib_(std::move(distrib))
    {
      for (const auto& dep:dependencies())
      {
        const auto& param=std::dynamic_pointer_cast<ConfiguredParameter>(dep);
        shareParameter_(param);
      }
    }

    ConfiguredDistribution::~ConfiguredDistribution () = default;

    std::string ConfiguredDistribution::description () const { return "Distribution(" + distrib_->getName () + ")"; }

    std::string ConfiguredDistribution::debugInfo () const {
      return "nbClass=" + std::to_string (distrib_->getNumberOfCategories ());
    }

    // Model node additional arguments = (type of bpp::TransitionModel).
    // Everything else is determined by the node dependencies.
    bool ConfiguredDistribution::compareAdditionalArguments (const Node & other) const {
      const auto * derived = dynamic_cast<const Self *> (&other);
      if (derived == nullptr) {
        return false;
      } else {
        const auto & thisDistrib = *distrib_;
        const auto & otherDistrib = *derived->distrib_;
        return typeid (thisDistrib) == typeid (otherDistrib);
      }
    }
    
    std::size_t ConfiguredDistribution::hashAdditionalArguments () const {
      const auto & bppDistrib = *distrib_;
      return typeid (bppDistrib).hash_code ();
    }

    NodeRef ConfiguredDistribution::recreate (Context & c, NodeRefVec && deps) {
      auto m = ConfiguredParametrizable::createConfigured<Target, Self> (c, std::move (deps), std::unique_ptr<DiscreteDistribution>{distrib_->clone ()});
      m->config = this->config; // Duplicate derivation config
      return m;
    }

    //////////////////////////////////////////////
    // ProbabilitiesFromDiscreteDistribution

    ProbabilitiesFromDiscreteDistribution::ProbabilitiesFromDiscreteDistribution (
      NodeRefVec && deps, const Dimension<Eigen::RowVectorXd> & dim)
      : Value<Eigen::RowVectorXd> (std::move (deps)), nbClass_ (dim) {}

    
    std::string ProbabilitiesFromDiscreteDistribution::debugInfo () const {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " nbClass=" + to_string (nbClass_);
    }

    // ProbabilitiesFromDiscreteDistribution additional arguments = ().
    bool ProbabilitiesFromDiscreteDistribution::compareAdditionalArguments (const Node & other) const {
      return dynamic_cast<const Self *> (&other) != nullptr;
    }

    NodeRef ProbabilitiesFromDiscreteDistribution::derive (Context & c, const Node & node) {
      // d(Prob)/dn = sum_i d(Prob)/dx_i * dx_i/dn (x_i = distrib parameters)
      auto distribDep = this->dependency (0);
      auto & distrib = static_cast<Dep &> (*distribDep);
      auto buildPWithNewDistrib = [this, &c](NodeRef && newDistrib) {
        return ConfiguredParametrizable::createVector<Dep, Self> (c, {std::move (newDistrib)}, nbClass_);
      };
      
      NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<ConfiguredDistribution, T > (
        c, distrib, node, nbClass_, buildPWithNewDistrib);
      return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), nbClass_);
    }

    NodeRef ProbabilitiesFromDiscreteDistribution::recreate (Context & c, NodeRefVec && deps) {
      return ConfiguredParametrizable::createVector<Dep, Self> (c, std::move (deps), nbClass_);
    }

    void ProbabilitiesFromDiscreteDistribution::compute () {
      const auto * distrib = accessValueConstCast<const DiscreteDistribution *> (*this->dependency (0));
      const auto & probasFromDistrib = distrib->getProbabilities ();
      auto & r = this->accessValueMutable ();
      r = Eigen::Map<const T> (probasFromDistrib.data(), static_cast<Eigen::Index> (probasFromDistrib.size ()));
    }


    std::shared_ptr<ProbabilitiesFromDiscreteDistribution> ProbabilitiesFromDiscreteDistribution::create(Context & c, NodeRefVec && deps)
    {
      checkDependenciesNotNull (typeid (Self), deps);
      checkDependencyVectorSize (typeid (Self), deps, 1);
      checkNthDependencyIs<ConfiguredDistribution> (typeid (Self), deps, 0);
      size_t nbCat=accessValueConstCast<DiscreteDistribution*> (*deps[0])->getNumberOfCategories();
      return cachedAs<ProbabilitiesFromDiscreteDistribution> (c, std::make_shared<ProbabilitiesFromDiscreteDistribution> (std::move(deps), rowVectorDimension(Eigen::Index(nbCat))));
    }
      

    ////////////////////////////////////////////////////
    // ProbabilityFromDiscreteDistribution

    ProbabilityFromDiscreteDistribution::ProbabilityFromDiscreteDistribution (
      NodeRefVec && deps, uint nCat)
      : Value<double> (std::move (deps)), nCat_ (nCat) {}

    
    std::string ProbabilityFromDiscreteDistribution::debugInfo () const {
      using namespace numeric;
      return "proba=" + TextTools::toString(accessValueConst()) + ":nCat=" + TextTools::toString(nCat_);
    }

    // ProbabilityFromDiscreteDistribution additional arguments = ().
    bool ProbabilityFromDiscreteDistribution::compareAdditionalArguments (const Node & other) const {
      const auto * derived = dynamic_cast<const Self *> (&other);
      return derived != nullptr && nCat_ == derived->nCat_;
    }

    std::shared_ptr<ProbabilityFromDiscreteDistribution> ProbabilityFromDiscreteDistribution::create (Context & c, NodeRefVec && deps, uint nCat) {
      checkDependenciesNotNull (typeid (Self), deps);
      checkDependencyVectorSize (typeid (Self), deps, 1);
      checkNthDependencyIs<ConfiguredDistribution> (typeid (Self), deps, 0);
      return cachedAs<ProbabilityFromDiscreteDistribution> (c, std::make_shared<ProbabilityFromDiscreteDistribution> (std::move (deps), nCat));
    }

    NodeRef ProbabilityFromDiscreteDistribution::derive (Context & c, const Node & node) {
      // d(Prob)/dn = sum_i d(Prob)/dx_i * dx_i/dn (x_i = distrib parameters)
      auto distribDep = this->dependency (0);
      auto & distrib = static_cast<Dep &> (*distribDep);
      auto buildPWithNewDistrib = [this, &c](NodeRef && newDistrib) {
        return this->create (c, {std::move (newDistrib)}, nCat_);
      };
      
      NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T > (
        c, distrib, node, 1, buildPWithNewDistrib);
      return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), 1);
    }

    NodeRef ProbabilityFromDiscreteDistribution::recreate (Context & c, NodeRefVec && deps) {
      return ProbabilityFromDiscreteDistribution::create (c, {std::move (deps)}, nCat_);
    }

    void ProbabilityFromDiscreteDistribution::compute () {
      const auto * distrib = accessValueConstCast<const DiscreteDistribution *> (*this->dependency (0));
      this->accessValueMutable () = distrib->getProbability((size_t)nCat_);
    }

    ////////////////////////////////////////////////////
    // CategoryFromDiscreteDistribution

    CategoryFromDiscreteDistribution::CategoryFromDiscreteDistribution (
      NodeRefVec && deps, uint nCat)
      : Value<double> (std::move (deps)), nCat_ (nCat) {}

    
    std::string CategoryFromDiscreteDistribution::debugInfo () const {
      using namespace numeric;
      return "rate=" + TextTools::toString(accessValueConst()) + ":nCat=" + TextTools::toString(nCat_);
    }

    // CategoryFromDiscreteDistribution additional arguments = ().
    bool CategoryFromDiscreteDistribution::compareAdditionalArguments (const Node & other) const {
      const auto * derived = dynamic_cast<const Self *> (&other);
      return derived != nullptr && nCat_ == derived->nCat_;
    }

    std::shared_ptr<CategoryFromDiscreteDistribution> CategoryFromDiscreteDistribution::create (Context & c, NodeRefVec && deps, uint nCat) {
      checkDependenciesNotNull (typeid (Self), deps);
      checkDependencyVectorSize (typeid (Self), deps, 1);
      checkNthDependencyIs<ConfiguredDistribution> (typeid (Self), deps, 0);
      return cachedAs<CategoryFromDiscreteDistribution> (c, std::make_shared<CategoryFromDiscreteDistribution> (std::move (deps), nCat));
    }

    NodeRef CategoryFromDiscreteDistribution::derive (Context & c, const Node & node) {
      // d(Prob)/dn = sum_i d(Prob)/dx_i * dx_i/dn (x_i = distrib parameters)
      auto distribDep = this->dependency (0);
      auto & distrib = static_cast<Dep &> (*distribDep);
      auto buildPWithNewDistrib = [this, &c](NodeRef && newDistrib) {
        return this->create (c, {std::move (newDistrib)}, nCat_);
      };
      
      NodeRefVec derivativeSumDeps = ConfiguredParametrizable::generateDerivativeSumDepsForComputations<Dep, T > (
        c, distrib, node, 1, buildPWithNewDistrib);
      return CWiseAdd<T, ReductionOf<T>>::create (c, std::move (derivativeSumDeps), 1);
    }

    NodeRef CategoryFromDiscreteDistribution::recreate (Context & c, NodeRefVec && deps) {
      return CategoryFromDiscreteDistribution::create (c, {std::move (deps)}, nCat_);
    }

    void CategoryFromDiscreteDistribution::compute () {
      const auto * distrib = accessValueConstCast<const DiscreteDistribution *> (*this->dependency (0));
      double categoryFromDistrib = distrib->getCategory(nCat_);
      this->accessValueMutable () = categoryFromDistrib;
    }

  } // namespace dataflow
} // namespace bpp
