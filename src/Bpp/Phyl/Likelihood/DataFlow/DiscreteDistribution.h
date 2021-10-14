//
// File: DiscreteDistribution.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: jeudi 25 octobre 2018, ÃÂ  17h 17
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_DISCRETEDISTRIBUTION_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_DISCRETEDISTRIBUTION_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlow.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Parametrizable.h>
#include <functional>
#include <unordered_map>

#include "Definitions.h"

namespace bpp
{
/* Likelihood discrete distribution.
 */

/** @brief Data flow node representing a DiscreteDistribution
 * configured with parameter values.
 *
 * This class wraps a bpp::DiscreteDistribution as a data flow node.
 * It depends on Value<double> nodes (one for each parameter
 * declared in the distribution).
 * It provides a dummy value representing the "distribution
 * configured by its parameters".
 *
 * The dummy value is implemented as a pointer to the internal
 * distribution for simplicity.
 */

class ConfiguredDistribution : public Value<const DiscreteDistribution*>,
  public AbstractParametrizable
{
  // private:
  //   Context& context_;

public:
  using Self = ConfiguredDistribution;
  using Target = DiscreteDistribution;

  ConfiguredDistribution (Context& context, NodeRefVec&& deps, std::unique_ptr<DiscreteDistribution>&& distrib);
  ~ConfiguredDistribution ();

  ConfiguredDistribution* clone() const
  {
    throw bpp::Exception("ConfiguredDistribution clone should not happen.");
  }

  std::string description () const final;
  std::string debugInfo () const final;
  std::string color() const final
  {
    return "blue";
  }

  bool compareAdditionalArguments (const Node_DF& other) const;

  std::size_t hashAdditionalArguments () const;

  /// Configuration for numerical derivation of computation nodes using this Model.
  NumericalDerivativeConfiguration config;

  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

  const ConfiguredParameter& getConfiguredParameter(const std::string& name)
  {
    return static_cast<const ConfiguredParameter&>(getParameter(name));
  }

private:
  void compute ()
  {
    distrib_->matchParametersValues(getParameters());
  }

  std::unique_ptr<DiscreteDistribution> distrib_;
};

/** Probabilities = f(DiscreteDistribution).
 * Probabilities: RowVector(nbClass).
 * DiscreteDistribution: ConfiguredDistribution.
 *
 * Node construction should be done with the create static method.
 */

class ProbabilitiesFromDiscreteDistribution : public Value<Eigen::RowVectorXd>
{
public:
  using Self = ProbabilitiesFromDiscreteDistribution;
  using Dep = ConfiguredDistribution;
  using T = Eigen::RowVectorXd;

  ProbabilitiesFromDiscreteDistribution (NodeRefVec&& deps, const Dimension<T>& dim);

  std::string debugInfo () const final;

  std::string color() const final
  {
    return "blue";
  }

  bool compareAdditionalArguments (const Node_DF& other) const final;

  NodeRef derive (Context& c, const Node_DF& node) final;
  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

private:
  void compute () final;

  Dimension<T> nbClass_;

public:
  static std::shared_ptr<Self> create (Context& c, NodeRefVec&& deps);
};

/** Proba = f(DiscreteDistribution, Category).
 * Proba: Double.
 * DiscreteDistribution: ConfiguredDistribution.
 * Category: number of the category
 *
 * Node construction should be done with the create static method.
 */

class ProbabilityFromDiscreteDistribution : public Value<double>
{
private:
  uint nCat_;

public:
  using Self = ProbabilityFromDiscreteDistribution;
  using Dep = ConfiguredDistribution;
  using T = double;

  ProbabilityFromDiscreteDistribution (NodeRefVec&& deps, uint nCat_);

  std::string debugInfo () const final;

  bool compareAdditionalArguments (const Node_DF& other) const final;

  NodeRef derive (Context& c, const Node_DF& node) final;
  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

  std::string color() const final
  {
    return "blue";
  }

private:
  void compute () final;

public:
  static std::shared_ptr<Self> create (Context& c, NodeRefVec&& deps, uint nCat);
};

/** Rate = f(DiscreteDistribution, Category).
 * Rate: Double.
 * DiscreteDistribution: ConfiguredDistribution.
 * Category: number of the category
 *
 * Node construction should be done with the create static method.
 */

class CategoryFromDiscreteDistribution : public Value<double>
{
private:
  uint nCat_;

public:
  using Self = CategoryFromDiscreteDistribution;
  using Dep = ConfiguredDistribution;
  using T = double;

  CategoryFromDiscreteDistribution (NodeRefVec&& deps, uint nCat_);

  std::string debugInfo () const final;

  bool compareAdditionalArguments (const Node_DF& other) const final;

  NodeRef derive (Context& c, const Node_DF& node) final;
  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

  std::string color() const final
  {
    return "blue";
  }

private:
  void compute () final;

public:
  static std::shared_ptr<Self> create (Context& c, NodeRefVec&& deps, uint nCat);
};
} // namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_DISCRETEDISTRIBUTION_H
