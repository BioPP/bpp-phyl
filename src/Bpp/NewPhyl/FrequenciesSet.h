//
// File: FrequenciesSet.h
// Authors:
//   Laurent Gueguen (2017)
// Created: jeudi 11 octobre 2018, à 06h 57
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
#ifndef BPP_NEWPHYL_FREQUENCIES_SET_H
#define BPP_NEWPHYL_FREQUENCIES_SET_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowNumeric.h>
#include <Bpp/Exceptions.h>
#include <functional>
#include <unordered_map>

namespace bpp {

  class FrequenciesSet;

  namespace dataflow {

    /** @brief Data flow node representing a Frequencies Set
     * configured with parameter values.
     *
     * This class wraps a bpp::FrequenciesSet as a data flow node.
     *
     * It depends on Value<double> nodes (one for each parameter
     * declared in the freq set).
     * It provides a dummy value representing the "frequencies set
     * configured by its parameters".
     *
     * The dummy value is implemented as a pointer to the internal
     * frequencies set for simplicity.
     */
    
    class ConfiguredFrequenciesSet : public Value<const FrequenciesSet*>// ,ConfiguredParametrizable
    {
    public:
      using Self = ConfiguredFrequenciesSet;

      /** Create a new model node from a dependency vector.
       * Model parameters are given by a dependency vector of Value<double> nodes.
       * The number and order of parameters is given by the FrequenciesSet internal ParameterList.
       */
      static std::shared_ptr<ConfiguredFrequenciesSet> create (Context & c, NodeRefVec && deps,
                                                               std::unique_ptr<FrequenciesSet> && freqset);

      ConfiguredFrequenciesSet (NodeRefVec && deps, std::unique_ptr<FrequenciesSet> && freqset);
      ~ConfiguredFrequenciesSet ();
      
      std::string description () const final;
      std::string debugInfo () const final;

      /// Return the index of parameter with the given non namespaced name (or throw).
      std::size_t getParameterIndex (const std::string & name);
      /// Return the non namespaced name for parameter at the given index.
      const std::string & getParameterName (std::size_t index);

      bool compareAdditionalArguments (const Node & other) const;
      
      std::size_t hashAdditionalArguments () const;
      
      /// Configuration for numerical derivation of computation nodes using this FrequenciesSet.
      NumericalDerivativeConfiguration config;

      NodeRef recreate (Context & c, NodeRefVec && deps) final;

    private:
      void compute ();

      std::unique_ptr<FrequenciesSet> freqset_;
    };

  } // namespace dataflow
} // namespace bpp

#endif // BPP_NEWPHYL_FREQUENCIES_SET_H
