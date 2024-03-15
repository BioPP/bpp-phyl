// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
  public virtual CodonSubstitutionModelInterface,
  public virtual AbstractKroneckerWordSubstitutionModel
{
private:
  std::shared_ptr<const GeneticCode> gCode_;

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
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
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
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    const std::vector<std::set<size_t> >& vPos,
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
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
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
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
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

  AbstractKroneckerCodonSubstitutionModel* clone() const override = 0;

protected:
  /**
   * @brief Method inherited from AbstractWordSubstitutionModel
   *
   * This method sets the rates to/from stop codons to zero and
   * performs the multiplication by the specific codon-codon rate.
   */
  void completeMatrices_() override;

public:
  std::shared_ptr<const GeneticCode> getGeneticCode() const override { return gCode_; }

  /**
   * @brief Method inherited from CodonSubstitutionModel
   *
   * Here this methods returns 1;
   *
   **/
  virtual double getCodonsMulRate(size_t i, size_t j) const override { return 1.; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTKRONECKERCODONSUBSTITUTIONMODEL_H
