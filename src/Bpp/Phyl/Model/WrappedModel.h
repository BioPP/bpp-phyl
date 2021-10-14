//
// File: WrappedModel.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: mardi 26 septembre 2017, ÃÂ  16h 18
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

#ifndef BPP_PHYL_MODEL_WRAPPEDMODEL_H
#define BPP_PHYL_MODEL_WRAPPEDMODEL_H

#include <Bpp/Seq/Container/SequencedValuesContainer.h>

#include "SubstitutionModel.h"

namespace bpp
{
/**
 * @brief Wrapping model class
 *
 *
 */

class WrappedModel :
  public virtual BranchModel
{
public:
  WrappedModel() {}
  virtual ~WrappedModel() {}

  virtual const BranchModel& getModel() const = 0;

protected:
  virtual BranchModel& getModel() = 0;

public:
};

class WrappedTransitionModel :
  public virtual WrappedModel,
  public virtual TransitionModel
{
public:
  WrappedTransitionModel() {}
  virtual ~WrappedTransitionModel() {}

  virtual const TransitionModel& getTransitionModel() const = 0;

protected:
  virtual TransitionModel& getTransitionModel() = 0;

public:
};


class WrappedSubstitutionModel :
  public virtual WrappedModel,
  public virtual SubstitutionModel
{
public:
  WrappedSubstitutionModel() {}

  virtual ~WrappedSubstitutionModel() {}

  virtual const SubstitutionModel& getSubstitutionModel() const = 0;

protected:
  virtual SubstitutionModel& getSubstitutionModel() = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_WRAPPEDMODEL_H
