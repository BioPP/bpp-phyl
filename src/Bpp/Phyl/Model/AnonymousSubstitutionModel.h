// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_ANONYMOUSSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_ANONYMOUSSUBSTITUTIONMODEL_H


#include "AbstractSubstitutionModel.h"

namespace bpp
{
/*
 * @brief Substitution Model with no name nor predefined generator
 * sets, and provides a non-const setGenerator() function to fill
 * the generator directly.
 *
 * Directly inherits from AbstractSubstitutionModel, hence uses the
 * computation methods developped in it.
 */
class AnonymousSubstitutionModel :
  public AbstractSubstitutionModel
{
public:
  AnonymousSubstitutionModel(
      std::shared_ptr<const Alphabet> alpha,
      std::shared_ptr<const StateMapInterface> stateMap) :
    AbstractParameterAliasable("Anonymous"),
    AbstractSubstitutionModel(alpha, stateMap, "Anonymous")
  {}

  virtual ~AnonymousSubstitutionModel() {}

  AnonymousSubstitutionModel* clone() const override
  {
    return new AnonymousSubstitutionModel(*this);
  }

public:
  std::string getName() const override { return "Anonymous"; }

  Matrix<double>& setGenerator() { return generator_; }

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ANONYMOUSSUBSTITUTIONMODEL_H
