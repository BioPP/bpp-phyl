// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
class NucleotideSubstitutionModelInterface :
  public virtual SubstitutionModelInterface
{
public:
  virtual ~NucleotideSubstitutionModelInterface() {}

  NucleotideSubstitutionModelInterface* clone() const override = 0;

public:

  virtual std::shared_ptr<const NucleicAlphabet> getNucleicAlphabet() const = 0;
};


/**
 * @brief Specialisation interface for rversible nucleotide substitution model.
 */
class NucleotideReversibleSubstitutionModelInterface :
  public virtual NucleotideSubstitutionModelInterface,
  public virtual ReversibleSubstitutionModelInterface
{
public:
  virtual ~NucleotideReversibleSubstitutionModelInterface() {}

  NucleotideReversibleSubstitutionModelInterface* clone() const override = 0;

};


/**
 * @brief Specialisation abstract class for nucleotide substitution model.
 */
class AbstractNucleotideSubstitutionModel :
  public AbstractSubstitutionModel,
  public virtual NucleotideSubstitutionModelInterface
{
public:
  AbstractNucleotideSubstitutionModel(
      std::shared_ptr<const NucleicAlphabet> alpha,
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& prefix) :
    AbstractSubstitutionModel(alpha, stateMap, prefix) {}

  virtual ~AbstractNucleotideSubstitutionModel() {}
  
  AbstractNucleotideSubstitutionModel* clone() const override = 0;

public:
  std::shared_ptr<const NucleicAlphabet> getNucleicAlphabet() const override
  {
    return std::dynamic_pointer_cast<const NucleicAlphabet>(alphabet_);
  }

};

/**
 * @brief Specialisation abstract class for reversible nucleotide substitution model.
 */
class AbstractReversibleNucleotideSubstitutionModel :
  public AbstractReversibleSubstitutionModel,
  public virtual NucleotideReversibleSubstitutionModelInterface
{
public:
  AbstractReversibleNucleotideSubstitutionModel(
      std::shared_ptr<const NucleicAlphabet> alpha,
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& prefix) :
    AbstractReversibleSubstitutionModel(alpha, stateMap, prefix) {}

  virtual ~AbstractReversibleNucleotideSubstitutionModel() {}

  AbstractReversibleNucleotideSubstitutionModel* clone() const override = 0;

public:
  std::shared_ptr<const NucleicAlphabet> getNucleicAlphabet() const override
  {
    return std::dynamic_pointer_cast<const NucleicAlphabet>(alphabet_);
  }

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_NUCLEOTIDESUBSTITUTIONMODEL_H
