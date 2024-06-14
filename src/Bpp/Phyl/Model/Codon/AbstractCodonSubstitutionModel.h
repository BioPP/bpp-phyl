// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_ABSTRACTCODONSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_ABSTRACTCODONSUBSTITUTIONMODEL_H


#include "../AbstractWordSubstitutionModel.h"
#include "../Nucleotide/NucleotideSubstitutionModel.h"
#include "CodonSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

// From the STL:
#include <memory>

namespace bpp
{
/**
 * @brief Abstract class for substitution models on codons.
 * @author Laurent Guéguen
 *
 * Objects of this class are built from either one (repeated three
 * times) or three different substitution models of NucleicAlphabets.
 * No model is directly accessible. </p>
 *
 * Only substitutions with one letter changed are accepted. </p>
 *
 * There is one substitution per codon per unit of time
 * on the equilibrium frequency, and each position has its specific rate.
 *
 * The parameters of this codon are the same as the ones of the models
 * used. Their names have a new prefix, "i_" where i stands for the
 * the phase (1,2 or 3) in the codon.
 */
class AbstractCodonSubstitutionModel :
  public virtual CodonSubstitutionModelInterface,
  public AbstractWordSubstitutionModel
{
private:
  /**
   * @brief boolean for the parametrization of the position relative
   * rates. Default : false.
   */
  bool hasParametrizedRates_;
  std::shared_ptr<const GeneticCode> gCode_;

public:
  /**
   * @brief Build a new AbstractCodonSubstitutionModel object from
   * a pointer to a NucleotideSubstitutionModel.
   *
   * @param gCode a pointer toward a genetic code. The codon alphabet from the genetic code will be used by the model class.
   * @param pmod pointer to the NucleotideSubstitutionModel to use in
   *        the three positions. It is owned by the instance.
   * @param st string of the Namespace
   * @param paramRates boolean concerning the presence of position
   * relative rates (default: false)
   */
  AbstractCodonSubstitutionModel(
      std::shared_ptr<const GeneticCode> gCode,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
      const std::string& st,
      bool paramRates = false);

  /**
   * @brief Build a new AbstractCodonSubstitutionModel object
   * from three pointers to NucleotideSubstitutionModels.
   *
   * @param gCode a pointer toward a genetic code. This model instance will own the underlying GeneticCode object and delete it when required.
   *   The codon alphabet from the genetic code will be used by the model class.
   * @param pmod1, pmod2, pmod3 are pointers to the
   *   NucleotideSubstitutionModel to use in the three positions.
   *   All the models must be different objects to avoid redundant
   *   parameters.  They are owned by the instance.
   * @param st string of the Namespace
   * @param paramRates boolean concerning the presence of position
   * relative rates (default: false)
   */
  AbstractCodonSubstitutionModel(
      std::shared_ptr<const GeneticCode> gCode,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
      const std::string& st,
      bool paramRates = false);

  virtual ~AbstractCodonSubstitutionModel() {}

  AbstractCodonSubstitutionModel(const AbstractCodonSubstitutionModel& model) :
    AbstractWordSubstitutionModel(model),
    hasParametrizedRates_(model.hasParametrizedRates_),
    gCode_(model.gCode_)
  {}

  AbstractCodonSubstitutionModel& operator=(const AbstractCodonSubstitutionModel& model)
  {
    AbstractWordSubstitutionModel::operator=(model);
    hasParametrizedRates_ = model.hasParametrizedRates_;
    gCode_ = model.gCode_;
    return *this;
  }

  AbstractCodonSubstitutionModel* clone() const override = 0;

  void setNamespace(const std::string& prefix) override
  {
    AbstractWordSubstitutionModel::setNamespace(prefix);
  }

protected:
  /**
   * @brief Method inherited from AbstractWordSubstitutionModel
   *
   * This method sets the rates to/from stop codons to zero and
   * performs the multiplication by the specific codon-codon rate.
   */
  void completeMatrices_() override;

  void updateMatrices_() override;

public:
  std::shared_ptr<const GeneticCode> getGeneticCode() const override { return gCode_; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_ABSTRACTCODONSUBSTITUTIONMODEL_H
