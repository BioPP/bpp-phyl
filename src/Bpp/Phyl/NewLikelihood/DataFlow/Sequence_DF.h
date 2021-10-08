//
// File: Sequence.h
// Authors:
//   Laurent Gueguen (2017)
// Created: vendredi 6 mars 2020, à 14h 16
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
#ifndef BPP_SEQUENCE_DF_H
#define BPP_SEQUENCE_DF_H


#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowNumeric.h>
#include <Bpp/Exceptions.h>
#include <functional>
#include <unordered_map>
#include "Definitions.h"

namespace bpp
{
/** @brief Data flow node representing a Sequence as a
 * Value<Eigen::MatrixXd> with a name.
 *
 */

class Sequence_DF : public Value<MatrixLik>
{
private:
  std::string name_;

public:
  using Self = Sequence_DF;
  using T = MatrixLik;

  static ValueRef<T> create (Context& c, T&& value, const std::string& name)
  {
    return cachedAs<Self>(c, std::make_shared<Self>(std::move(value), name));
  }

  Sequence_DF (T&& value, const std::string& name) :
    Value<T>(NodeRefVec{}, std::move(value)),
    name_(name)
  {
    this->makeValid (); // Always valid
  }

  const std::string& getName() const
  {
    return name_;
  }

  std::string debugInfo () const final
  {
    using namespace numeric;
    return name_ + " " + debug(this->accessValueConst ());
  }

  std::string description () const final
  {
    using namespace numeric;
    return Node_DF::description() + "\n" + name_ + "\n" + debug (this->accessValueConst ());
  }

  bool compareAdditionalArguments (const Node_DF& other) const override
  {
    const auto* derived = dynamic_cast<const Self*>(&other);
    return derived != nullptr && name_ == derived->name_ && this->accessValueConst () == derived->accessValueConst ();
  }

  std::string color () const override
  {
    return "grey";
  }

  NodeRef recreate (Context& c, NodeRefVec&& deps) final
  {
    checkRecreateWithoutDependencies (typeid (Self), deps);
    return this->shared_from_this ();
  }

  std::size_t hashAdditionalArguments () const override
  {
    using namespace numeric;
    size_t seed = hash (this->accessValueConst ());
    combineHash<std::string>(seed, name_);
    return seed;
  }

  NodeRef derive (Context& c, const Node_DF& node) final
  {
    const auto dim = Dimension<T>(this->accessValueConst ());
    if (&node == this)
    {
      return ConstantOne<T>::create (c, dim);
    }
    return ConstantZero<T>::create (c, dim);
  }

private:
  void compute () final
  {
    // Constant is valid from construction
    failureComputeWasCalled (typeid (*this));
  }
};
} // namespace bpp


#endif// BPP_SEQUENCE_DF_H
