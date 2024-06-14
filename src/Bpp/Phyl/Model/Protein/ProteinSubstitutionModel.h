// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
class ProteinSubstitutionModelInterface :
  public virtual SubstitutionModelInterface
{
public:
  virtual ~ProteinSubstitutionModelInterface() {}

  virtual ProteinSubstitutionModelInterface* clone() const override = 0;

  virtual std::shared_ptr<const ProteicAlphabet> getProteicAlphabet() const = 0;
};


/**
 * @brief Specialized interface for protein reversible substitution model.
 */
class ProteinReversibleSubstitutionModelInterface :
  public virtual ProteinSubstitutionModelInterface,
  public virtual ReversibleSubstitutionModelInterface
{
public:
  virtual ~ProteinReversibleSubstitutionModelInterface() {}

  ProteinReversibleSubstitutionModelInterface* clone() const override = 0;
};


/**
 * @brief Specialisation abstract class for protein substitution model.
 */
class AbstractProteinSubstitutionModel :
  public AbstractSubstitutionModel,
  public virtual ProteinSubstitutionModelInterface
{
public:
  AbstractProteinSubstitutionModel(
      std::shared_ptr<const ProteicAlphabet> alpha,
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& prefix) :
    AbstractSubstitutionModel(alpha, stateMap, prefix) {}

  virtual ~AbstractProteinSubstitutionModel() {}

  AbstractProteinSubstitutionModel* clone() const override = 0;

public:
  std::shared_ptr<const ProteicAlphabet> getProteicAlphabet() const override
  {
    return std::dynamic_pointer_cast<const ProteicAlphabet>(alphabet_);
  }
};

/**
 * @brief Specialisation abstract class for reversible protein substitution model.
 */
class AbstractReversibleProteinSubstitutionModel :
  public AbstractReversibleSubstitutionModel,
  public virtual ProteinReversibleSubstitutionModelInterface
{
public:
  AbstractReversibleProteinSubstitutionModel(
      std::shared_ptr<const ProteicAlphabet> alpha,
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& prefix) :
    AbstractReversibleSubstitutionModel(alpha, stateMap, prefix) {}

  virtual ~AbstractReversibleProteinSubstitutionModel() {}

  AbstractReversibleProteinSubstitutionModel* clone() const override = 0;

public:
  std::shared_ptr<const ProteicAlphabet> getProteicAlphabet() const override
  {
    return std::dynamic_pointer_cast<const ProteicAlphabet>(alphabet_);
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_PROTEINSUBSTITUTIONMODEL_H
