//
// File: AbstractKroneckerCodonSubstitutionModel.h
// Authors:
//   Laurent Gueguen
// Created: mardi 26 juillet 2016, ÃÂ  21h 15
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

#ifndef BPP_PHYL_MODEL_CODON_ABSTRACTKRONECKERCODONSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_ABSTRACTKRONECKERCODONSUBSTITUTIONMODEL_H


#include "../AbstractKroneckerWordSubstitutionModel.h"
#include "../Nucleotide/NucleotideSubstitutionModel.h"
#include "CodonSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

// From the STL:
#include <memory>
#include <set>

namespace bpp
{
/**
 * @brief Abstract class for substitution models on codons allowing
 * multiple substitutions.
 *
 * Rates of multiple substitutions equal the product of single
 * substitutions involved, before removing stop codons.
 *
 * @author Laurent GuÃÂ©guen
 *
 * Objects of this class are built from either one (repeated three
 * times) or three different substitution models of NucleicAlphabets.
 * No model is directly accessible. </p>
 *
 * The parameters of this codon are the same as the ones of the models
 * used. Their names have a new prefix, "i_" where i stands for the
 * the phase (1,2 or 3) in the codon.
 */

class AbstractKroneckerCodonSubstitutionModel :
  public virtual CodonSubstitutionModel,
  public AbstractKroneckerWordSubstitutionModel
{
private:
  const GeneticCode* gCode_;

public:
  /**
   * @brief Build a new AbstractKroneckerCodonSubstitutionModel object from
   * a pointer to a NucleotideSubstitutionModel.
   *
   * @param gCode a pointer toward a genetic code. The codon alphabet from the genetic code will be used by the model class.
   * @param pmod pointer to the NucleotideSubstitutionModel to use in
   *        the three positions. It is owned by the instance.
   * @param st string of the Namespace
   */

  AbstractKroneckerCodonSubstitutionModel(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod,
    const std::string& st);

  /**
   * @brief Build a new AbstractKroneckerCodonSubstitutionModel object from
   * a pointer to a NucleotideSubstitutionModel.
   *
   * @param gCode a pointer toward a genetic code. The codon alphabet from the genetic code will be used by the model class.
   * @param pmod pointer to the NucleotideSubstitutionModel to use in
   *        the three positions. It is owned by the instance.
   * @param vPos a vector of sets of simultaneously changing
   *   positions.
   * @param st string of the Namespace
   */

  AbstractKroneckerCodonSubstitutionModel(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod,
    const std::vector<std::set< size_t> >& vPos,
    const std::string& st);

  /**
   * @brief Build a new AbstractKroneckerCodonSubstitutionModel object
   * from three pointers to NucleotideSubstitutionModels.
   *
   * @param gCode a pointer toward a genetic code. This model instance will own the underlying GeneticCode object and delete it when required.
   *   The codon alphabet from the genetic code will be used by the model class.
   * @param pmod1, pmod2, pmod3 are pointers to the
   *   NucleotideSubstitutionModel to use in the three positions.
   *   All the models must be different objects to avoid redundant
   *   parameters.  They are owned by the instance.
   * @param st string of the Namespace
   */

  AbstractKroneckerCodonSubstitutionModel(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod1,
    NucleotideSubstitutionModel* pmod2,
    NucleotideSubstitutionModel* pmod3,
    const std::string& st);

  /**
   * @brief Build a new AbstractKroneckerCodonSubstitutionModel object
   * from three pointers to NucleotideSubstitutionModels.
   *
   * @param gCode a pointer toward a genetic code. This model instance will own the underlying GeneticCode object and delete it when required.
   *   The codon alphabet from the genetic code will be used by the model class.
   * @param pmod1, pmod2, pmod3 are pointers to the
   *   NucleotideSubstitutionModel to use in the three positions.
   *   All the models must be different objects to avoid redundant
   *   parameters.  They are owned by the instance.
   * @param vPos a vector of sets of simultaneously changing
   *   positions.
   * @param st string of the Namespace
   */

  AbstractKroneckerCodonSubstitutionModel(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod1,
    NucleotideSubstitutionModel* pmod2,
    NucleotideSubstitutionModel* pmod3,
    const std::vector<std::set< size_t> >& vPos,
    const std::string& st);

  virtual ~AbstractKroneckerCodonSubstitutionModel() {}

  AbstractKroneckerCodonSubstitutionModel(const AbstractKroneckerCodonSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    AbstractKroneckerWordSubstitutionModel(model),
    gCode_(model.gCode_)
  {}

  AbstractKroneckerCodonSubstitutionModel& operator=(const AbstractKroneckerCodonSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractKroneckerWordSubstitutionModel::operator=(model);
    gCode_ = model.gCode_;
    return *this;
  }

  AbstractKroneckerCodonSubstitutionModel* clone() const = 0;

protected:
  /**
   * @brief Method inherited from AbstractWordSubstitutionModel
   *
   * This method sets the rates to/from stop codons to zero and
   * performs the multiplication by the specific codon-codon rate.
   */
  void completeMatrices();

public:
  const GeneticCode* getGeneticCode() const { return gCode_; }

  /**
   * @brief Method inherited from CodonSubstitutionModel
   *
   * Here this methods returns 1;
   *
   **/
  virtual double getCodonsMulRate(size_t i, size_t j) const { return 1.; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTKRONECKERCODONSUBSTITUTIONMODEL_H
