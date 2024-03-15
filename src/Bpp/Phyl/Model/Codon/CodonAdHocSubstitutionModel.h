// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_CODONADHOCSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_CODONADHOCSUBSTITUTIONMODEL_H


#include "AbstractCodonSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Class for substitution models of codons with
 * several layers of codon models
 *
 * @author Laurent Guéguen
 *
 * Objects of this class are built from three substitution models of
 * NucleicAlphabets. No model is directly accessible. </p>
 *
 * Only substitutions with one letter changed are accepted. </p>
 */
class CodonAdHocSubstitutionModel :
  public AbstractCodonSubstitutionModel
{
private:
  std::vector< std::unique_ptr<CoreCodonSubstitutionModelInterface> > vModel_;

  std::string name_;

  /**
   * @brief optional FrequencySet if model is defined through a
   * FrequencySet.
   */
  std::unique_ptr<CodonFrequencySetInterface> freqSet_;

public:
  /**
   * @brief Build a new CodonAdHocSubstitutionModel object from
   * a pointer to NucleotideSubstitutionModel.
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod  pointer to the NucleotideSubstitutionModel to use
   * in the three positions.
   * The instance will then own this substitution model.
   * @param vpmodel vector of codon models. They will be owned by the model.
   * @param name the name of the model
   */
  CodonAdHocSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    std::vector<std::unique_ptr<CoreCodonSubstitutionModelInterface>>& vpmodel,
    const std::string& name);

  /**
   * @brief Build a new CodonAdHocSubstitutionModel object
   * from three pointers to NucleotideSubstitutionModels.
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod1, pmod2, pmod3 pointers to the
   *   NucleotideSubstitutionModels to use in the three positions.
   * @param vpmodel vector of codon models. They will be owned by the
   * model.
   * @param name the name of the model
   */
  CodonAdHocSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    std::vector<std::unique_ptr<CoreCodonSubstitutionModelInterface>>& vpmodel,
    const std::string& name);

  CodonAdHocSubstitutionModel(const CodonAdHocSubstitutionModel& model);

  CodonAdHocSubstitutionModel& operator=(const CodonAdHocSubstitutionModel& model);

  virtual ~CodonAdHocSubstitutionModel() {}

  CodonAdHocSubstitutionModel* clone() const override
  {
    return new CodonAdHocSubstitutionModel(*this);
  }

public:
  void fireParameterChanged(const ParameterList& parameterlist) override;

  std::string getName() const override
  {
    return name_;
  }

  void setNamespace(const std::string& prefix) override
  {
    AbstractCodonSubstitutionModel::setNamespace(prefix);
    for (auto& model : vModel_)
    {
      model->setNamespace(prefix);
    }
  }

  size_t getNumberOfModels() const
  {
    return vModel_.size();
  }

  const CoreCodonSubstitutionModelInterface& layerModel(size_t i) const
  {
    return *vModel_[i];
  }

  double getCodonsMulRate(size_t i, size_t j) const override;

  void setFreq(std::map<int, double>& frequencies) override;

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    return *freqSet_;
  }

  bool hasCodonFrequencySet() const override
  {
    return (freqSet_ != nullptr);
  }

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_CODONADHOCSUBSTITUTIONMODEL_H
