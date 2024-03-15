// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_WRAPPEDMODEL_H
#define BPP_PHYL_MODEL_WRAPPEDMODEL_H

#include <Bpp/Seq/Container/SequencedValuesContainer.h>

#include "SubstitutionModel.h"

namespace bpp
{
/**
 * @brief Wrapping model interface
 */
class WrappedModelInterface :
  public virtual BranchModelInterface
{
public:
  WrappedModelInterface() {}
  virtual ~WrappedModelInterface() {}

  virtual const BranchModelInterface& model() const = 0;

public:
};

class WrappedTransitionModelInterface :
  public virtual WrappedModelInterface,
  public virtual TransitionModelInterface
{
public:
  WrappedTransitionModelInterface() {}
  virtual ~WrappedTransitionModelInterface() {}

  virtual const TransitionModelInterface& transitionModel() const = 0;

};


class WrappedSubstitutionModelInterface :
  public virtual WrappedModelInterface,
  public virtual SubstitutionModelInterface
{
public:
  WrappedSubstitutionModelInterface() {}

  virtual ~WrappedSubstitutionModelInterface() {}

  virtual const SubstitutionModelInterface& substitutionModel() const = 0;

};

} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_WRAPPEDMODEL_H
