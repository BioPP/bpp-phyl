//
// File: DiscreteDistribution.h
// Authors: Laurent Guéguen
// Created: jeudi 25 octobre 2018, à 17h 17
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

#pragma once
#ifndef BPP_NEWPHYL_DISCRETE_DISTRIBUTION_H
#define BPP_NEWPHYL_DISCRETE_DISTRIBUTION_H 

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/NewPhyl/Parametrizable.h>
#include <Bpp/Exceptions.h>
#include <functional>
#include <unordered_map>

namespace bpp {
  /* Likelihood discrete distribution.
   */
  class DiscreteDistribution;

  namespace dataflow {

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
    
    class ConfiguredDistribution : public Value<const DiscreteDistribution*>
    {
    public:
      using Self = ConfiguredDistribution;
      using Target = DiscreteDistribution;
      
      ConfiguredDistribution (NodeRefVec && deps, std::unique_ptr<DiscreteDistribution> && distrib);
      ~ConfiguredDistribution ();

      std::string description () const final;
      std::string debugInfo () const final;

      /// Return the index of parameter with the given non namespaced name (or throw).
      std::size_t getParameterIndex (const std::string & name);
      /// Return the non namespaced name for parameter at the given index.
      const std::string & getParameterName (std::size_t index);

      bool compareAdditionalArguments (const Node & other) const;
      
      std::size_t hashAdditionalArguments () const;
      
      /// Configuration for numerical derivation of computation nodes using this Model.
      NumericalDerivativeConfiguration config;

      NodeRef recreate (Context & c, NodeRefVec && deps) final;
      
    private:
      void compute ();

      std::unique_ptr<DiscreteDistribution> distrib_;
      
    };

    /** Probabilities = f(DiscreteDistribution).
     * Probabilities: RowVector(nbClass).
     * DiscreteDistribution: ConfiguredDistribution.
     *
     * Node construction should be done with the create static method.
     */

    class ProbabilitiesFromDiscreteDistribution : public Value<Eigen::RowVectorXd> {
    public:
      using Self = ProbabilitiesFromDiscreteDistribution;
      using Dep = ConfiguredDistribution;
      using T = Eigen::RowVectorXd;

      ProbabilitiesFromDiscreteDistribution (NodeRefVec && deps, const Dimension<T> & dim);

      std::string debugInfo () const final;

      bool compareAdditionalArguments (const Node & other) const final;

      NodeRef derive (Context & c, const Node & node) final;
      NodeRef recreate (Context & c, NodeRefVec && deps) final;

    private:
      void compute () final;

      Dimension<T> nbClass_;
    };

    /** Categories = f(DiscreteDistribution).
     * Categories: RowVector(nbClass).
     * DiscreteDistribution: ConfiguredDistribution.
     *
     * Node construction should be done with the create static method.
     */

    class CategoriesFromDiscreteDistribution : public Value<Eigen::RowVectorXd> {
    public:
      using Self = CategoriesFromDiscreteDistribution;
      using Dep = ConfiguredDistribution;
      using T = Eigen::RowVectorXd;

      CategoriesFromDiscreteDistribution (NodeRefVec && deps, const Dimension<T> & dim);

      std::string debugInfo () const final;

      bool compareAdditionalArguments (const Node & other) const final;

      NodeRef derive (Context & c, const Node & node) final;
      NodeRef recreate (Context & c, NodeRefVec && deps) final;

    private:
      void compute () final;

      Dimension<T> nbClass_;
    };

  } // namespace dataflow
} // namespace bpp

#endif // BPP_NEWPHYL_DISCRETE_DISTRIBUTION_H
