//
// File: CodonAdHocSubstitutionModel.h
// Authors:
//   Laurent Gueguen
// Created: lundi 30 octobre 2017, ÃÂ  06h 31
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
