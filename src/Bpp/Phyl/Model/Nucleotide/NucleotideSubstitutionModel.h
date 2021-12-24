//
// File: NucleotideSubstitutionModel.h
// Authors:
//   Julien Dutheil
// Created: 2003-05-27 11:03:53
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

#ifndef BPP_PHYL_MODEL_NUCLEOTIDE_NUCLEOTIDESUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_NUCLEOTIDE_NUCLEOTIDESUBSTITUTIONMODEL_H


#include "../AbstractSubstitutionModel.h"
#include "../SubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>

namespace bpp
{
/**
 * @brief Specialisation interface for nucleotide substitution model.
 */
class NucleotideSubstitutionModel :
  public virtual SubstitutionModel
{
public:
  virtual ~NucleotideSubstitutionModel() {}

  NucleotideSubstitutionModel* clone() const override = 0;

public:
  size_t getNumberOfStates() const override { return 4; }

  const NucleicAlphabet* getAlphabet() const override = 0;
};


/**
 * @brief Specialisation interface for rversible nucleotide substitution model.
 */
class NucleotideReversibleSubstitutionModel :
  public virtual NucleotideSubstitutionModel,
  public virtual ReversibleSubstitutionModel
{
public:
  virtual ~NucleotideReversibleSubstitutionModel() {}

  NucleotideReversibleSubstitutionModel* clone() const override = 0;

  size_t getNumberOfStates() const override
  {
    return NucleotideSubstitutionModel::getNumberOfStates();
  }
};


/**
 * @brief Specialisation abstract class for nucleotide substitution model.
 */
class AbstractNucleotideSubstitutionModel :
  public AbstractSubstitutionModel,
  public virtual NucleotideSubstitutionModel
{
public:
  AbstractNucleotideSubstitutionModel(const NucleicAlphabet* alpha, std::shared_ptr<const StateMap> stateMap, const std::string& prefix) :
    AbstractParameterAliasable(prefix),
    AbstractSubstitutionModel(alpha, stateMap, prefix) {}

  virtual ~AbstractNucleotideSubstitutionModel() {}

  AbstractNucleotideSubstitutionModel* clone() const override = 0;

public:
  const NucleicAlphabet* getAlphabet() const override
  {
    return dynamic_cast<const NucleicAlphabet*>(alphabet_);
  }

  size_t getNumberOfStates() const override
  {
    return NucleotideSubstitutionModel::getNumberOfStates();
  }
};

/**
 * @brief Specialisation abstract class for reversible nucleotide substitution model.
 */
class AbstractReversibleNucleotideSubstitutionModel :
  public AbstractReversibleSubstitutionModel,
  public virtual NucleotideReversibleSubstitutionModel
{
public:
  AbstractReversibleNucleotideSubstitutionModel(const NucleicAlphabet* alpha, std::shared_ptr<const StateMap> stateMap, const std::string& prefix) :
    AbstractParameterAliasable(prefix),
    AbstractReversibleSubstitutionModel(alpha, stateMap, prefix) {}

  virtual ~AbstractReversibleNucleotideSubstitutionModel() {}

  AbstractReversibleNucleotideSubstitutionModel* clone() const override = 0;

public:
  const NucleicAlphabet* getAlphabet() const override
  {
    return dynamic_cast<const NucleicAlphabet*>(alphabet_);
  }

  size_t getNumberOfStates() const override
  {
    return NucleotideSubstitutionModel::getNumberOfStates();
  }

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_NUCLEOTIDESUBSTITUTIONMODEL_H
