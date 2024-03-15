// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_CODONSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_CODONSUBSTITUTIONMODEL_H



// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include "../SubstitutionModel.h"
#include "../FrequencySet/CodonFrequencySet.h"

namespace bpp
{
/**
 * @brief Interface for codon models
 * @author Laurent Guéguen
 *
 * This class aims at defining methods needed for inheriting codon.
 */
class CoreCodonSubstitutionModelInterface :
  public virtual ParameterAliasable
{
public:
  CoreCodonSubstitutionModelInterface() {}
  virtual ~CoreCodonSubstitutionModelInterface() {}

  virtual CoreCodonSubstitutionModelInterface* clone() const override = 0;

public:
  /**
   * @brief Returns the multiplicative rate specific to two codons
   * specified by their number. The respective generator rate is this
   * rate multiplied by the rate defined by the model defined on
   * nucleotides.
   */
  virtual double getCodonsMulRate(size_t, size_t) const = 0;

  virtual const CodonFrequencySetInterface& codonFrequencySet() const = 0;
  
  virtual bool hasCodonFrequencySet() const = 0;
  
  virtual void setFreq(std::map<int, double>& frequencies) = 0;
};


class CodonSubstitutionModelInterface :
  public virtual CoreCodonSubstitutionModelInterface,
  public virtual SubstitutionModelInterface
{
public:
  CodonSubstitutionModelInterface() {}
  virtual ~CodonSubstitutionModelInterface() {}

  virtual CodonSubstitutionModelInterface* clone() const override = 0;

  virtual std::shared_ptr<const GeneticCode> getGeneticCode() const = 0;
};


/**
 * @brief Interface for reversible codon models
 * @author Julien Dutheil
 */
class CodonReversibleSubstitutionModelInterface :
  public virtual CodonSubstitutionModelInterface,
  public virtual ReversibleSubstitutionModelInterface
{
public:
  CodonReversibleSubstitutionModelInterface() {}
  virtual ~CodonReversibleSubstitutionModelInterface() {}

  virtual CodonReversibleSubstitutionModelInterface* clone() const override = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_CODONSUBSTITUTIONMODEL_H
