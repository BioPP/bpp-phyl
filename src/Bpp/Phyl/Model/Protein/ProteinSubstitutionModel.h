//
// File: ProteinSubstitutionModel.h
// Authors:
//   Julien Dutheil
// Created: 2004-01-21 13:59:18
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

#ifndef BPP_PHYL_MODEL_PROTEIN_PROTEINSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_PROTEIN_PROTEINSUBSTITUTIONMODEL_H


#include "../AbstractSubstitutionModel.h"
#include "../SubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>

namespace bpp
{
/**
 * @brief Specialized interface for protein substitution model.
 */
class ProteinSubstitutionModel :
  public virtual SubstitutionModel
{
public:
  virtual ~ProteinSubstitutionModel() {}

  ProteinSubstitutionModel* clone() const = 0;

public:
  size_t getNumberOfStates() const { return 20; }

  const ProteicAlphabet* getAlphabet() const = 0;
};


/**
 * @brief Specialized interface for protein reversible substitution model.
 */
class ProteinReversibleSubstitutionModel :
  public virtual ProteinSubstitutionModel,
  public virtual ReversibleSubstitutionModel
{
public:
  virtual ~ProteinReversibleSubstitutionModel() {}

  ProteinReversibleSubstitutionModel* clone() const = 0;
};


/**
 * @brief Specialisation abstract class for protein substitution model.
 */
class AbstractProteinSubstitutionModel :
  public AbstractSubstitutionModel,
  public virtual ProteinSubstitutionModel
{
public:
  AbstractProteinSubstitutionModel(const ProteicAlphabet* alpha, std::shared_ptr<const StateMap> stateMap, const std::string& prefix) :
    AbstractParameterAliasable(prefix),
    AbstractSubstitutionModel(alpha, stateMap, prefix) {}

  virtual ~AbstractProteinSubstitutionModel() {}

  AbstractProteinSubstitutionModel* clone() const = 0;

public:
  const ProteicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const ProteicAlphabet*>(alphabet_);
  }
};

/**
 * @brief Specialisation abstract class for reversible protein substitution model.
 */
class AbstractReversibleProteinSubstitutionModel :
  public AbstractReversibleSubstitutionModel,
  public virtual ProteinSubstitutionModel
{
public:
  AbstractReversibleProteinSubstitutionModel(const ProteicAlphabet* alpha, std::shared_ptr<const StateMap> stateMap, const std::string& prefix) :
    AbstractParameterAliasable(prefix),
    AbstractReversibleSubstitutionModel(alpha, stateMap, prefix) {}

  virtual ~AbstractReversibleProteinSubstitutionModel() {}

  AbstractReversibleProteinSubstitutionModel* clone() const = 0;

public:
  const ProteicAlphabet* getAlphabet() const
  {
    return dynamic_cast<const ProteicAlphabet*>(alphabet_);
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_PROTEINSUBSTITUTIONMODEL_H
