//
// File: CodonSubstitutionModel.h
// Authors:
//   Laurent Gueguen
// Created: 2003-12-24 11:03:53
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
