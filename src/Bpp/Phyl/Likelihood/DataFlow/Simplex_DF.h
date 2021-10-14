//
// File: Simplex_DF.h
// Authors:
//   Laurent Gueguen (2017)
// Created: jeudi 11 octobre 2018, ÃÂ  06h 57
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

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_SIMPLEX_DF_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_SIMPLEX_DF_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Prob/Simplex.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlow.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowCWiseComputing.h>
#include <functional>
#include <unordered_map>

#include "Definitions.h"

namespace bpp
{
class Simplex;

/** @brief Data flow node representing a Frequencies Set
 * configured with parameter values.
 *
 * This class wraps a bpp::Simplex as a data flow node.
 *
 * It depends on Value<double> nodes (one for each parameter
 * declared in the freq set).
 * It provides a dummy value representing the "frequencies set
 * configured by its parameters".
 *
 * The dummy value is implemented as a pointer to the internal
 * frequencies set for simplicity.
 */

class ConfiguredSimplex : public Value<const Simplex*>,
  public AbstractParametrizable
{
  // private:

  //   const Context& context_;

public:
  using Self = ConfiguredSimplex;
  using Target = Simplex;

  ConfiguredSimplex (const Context& context, NodeRefVec&& deps, std::unique_ptr<Simplex>&& simplex);
  ~ConfiguredSimplex ();

  ConfiguredSimplex* clone() const
  {
    throw bpp::Exception("ConfiguredSimplex clone should not happen.");
  }

  std::string description () const final;
  std::string debugInfo () const final;
  std::string color() const final
  {
    return "blue";
  }

  bool compareAdditionalArguments (const Node_DF& other) const;

  std::size_t hashAdditionalArguments () const;

  /// Configuration for numerical derivation of computation nodes using this Simplex.
  NumericalDerivativeConfiguration config;

  NodeRef recreate (Context& c, NodeRefVec&& deps) final;

  const ConfiguredParameter& getConfiguredParameter(const std::string& name)
  {
    return static_cast<const ConfiguredParameter&>(getParameter(name));
  }

private:
  void compute ()
  {
    simplex_->matchParametersValues(getParameters());
  }


  std::unique_ptr<Simplex> simplex_;
};

/** Frequencies = f(Simplex).
 * Frequencies: RowVector(nbState).
 * Frequenciesset: ConfiguredSimplex.
 *
 * Node construction should be done with the create static method.
 */

class FrequenciesFromSimplex : public Value<Eigen::RowVectorXd>
{
public:
  using Self = FrequenciesFromSimplex;
  using T = Eigen::RowVectorXd;

  // static ValueRef<T> create (Context & c, NodeRefVec && deps, const Dimension<T> & dim);
  FrequenciesFromSimplex (NodeRefVec&& deps, const Dimension<T>& dim);

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

  Dimension<T> targetDimension_;

public:
  static std::shared_ptr<Self> create (Context& c, NodeRefVec&& deps);
};
} // namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_SIMPLEX_DF_H
