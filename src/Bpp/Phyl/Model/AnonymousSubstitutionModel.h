//
// File: AnonymousSubstitutionModel.h
// Authors:
//   Laurent Gueguen
// Created: samedi 17 juin 2017, à 07h 55
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
 *
 *
 */


class AnonymousSubstitutionModel :
  public AbstractSubstitutionModel
{
public:
  AnonymousSubstitutionModel(const Alphabet* alpha, std::shared_ptr<const StateMap> stateMap) :
    AbstractParameterAliasable("Anonymous"),
    AbstractSubstitutionModel(alpha, stateMap, "Anonymous")
  {}

  virtual ~AnonymousSubstitutionModel() {}

  AnonymousSubstitutionModel* clone() const { return new AnonymousSubstitutionModel(*this); }

public:
  std::string getName() const { return "Anonymous"; }

  size_t getNumberOfStates() const { return size_;}

  Matrix<double>& setGenerator() { return generator_; }

  void updateMatrices()
  {
    AbstractSubstitutionModel::updateMatrices();
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ANONYMOUSSUBSTITUTIONMODEL_H
