// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_CODON_CODONSAMEAARATESUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_CODONSAMEAARATESUBSTITUTIONMODEL_H


#include "../FrequencySet/CodonFrequencySet.h"
#include "../Protein/ProteinSubstitutionModel.h"
#include "CodonSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Class for modelling of non-synonymous rates in
 *   codon models, such that the substitution rates between amino
 *   acids are similar to the ones in an amino acid rate matrix (from
 *   a shared_ptr model).
 *
 * @author Laurent Guéguen
 *
 * From the generator @f$A@f$ between amino-acids, the generator
 *  @f$Q@f$ between codons, and the codon frequencies @f$\pi@f$, the
 *  non-synonymous rate between codons @f$i@f$ and @f$j@f$ is
 *  multiplied by @f$x_{a_i a_j}@f$ (where @f$a_i@f$ is for amino acid
 *  encoded by codon @f$i@f$) such that:
 *
 * \f[
 * \sum_{l;a_l=a_i} \pi_l A_{a_i,a_j} = \sum_{l;a_l=a_i} \sum_{s;a_s=a_j} \pi_l Q_{ls} x_{a_i,a_j}
 * \f]
 *
 * Which means, with @f$ \phi_a \sum_{l;a_l=a} \pi_l @f$
 *
 * \f[
 * x_{a_i,a_j} = \frac{ \phi_{a_i} A_{a_i,a_j}}{{\sum_{l;a_l=a_i} \sum_{s;a_s=a_j} \pi_l Q_{ls}}}
 * \f]
 */
class CodonSameAARateSubstitutionModel :
  public virtual CodonSubstitutionModelInterface,
  public AbstractSubstitutionModel
{
private:
  /**
   * @brief Protein Model which will be used to get similar AA
   * rates.
   */
  std::unique_ptr<ProteinSubstitutionModelInterface> pAAmodel_;

  /**
   * @brief Codon Model which will be copied. Its possible
   * parameters are not copied in this object.
   *
   */
  std::unique_ptr<CodonSubstitutionModelInterface> pCodonModel_;

  /**
   * @brief Protein Model which will be used to get similar AA
   * rates. Its possible parameters are not copied in this object.
   *
   * May be null, if pi is the equilibrium frequency of pCodonModel_.
   */
  std::unique_ptr<CodonFrequencySetInterface> pFreq_;

  std::shared_ptr<const GeneticCode> pgencode_;

  /**
   * @brief 20 x 20 Matrix of the denominator of the multiplicators
   */
  RowMatrix<double> X_;

  /**
   * @brief 20 Vdouble for computation
   */
  Vdouble phi_;

public:
  /**
   * @brief Build a new CodonSameAARateSubstitutionModel object from
   *  a pointer to NucleotideSubstitutionModel.
   *
   * @param pAAmodel shared_ptr to an amino_acid generator
   * @param pCodonModel the codon substitution model
   *
   * @param pFreq the Codon Frequency set, may be null, in which
   * case the equilibrium frequencies of the model are used.
   *
   * @param pgencode the genetic code
   */
  CodonSameAARateSubstitutionModel(
      std::unique_ptr<ProteinSubstitutionModelInterface> pAAmodel,
      std::unique_ptr<CodonSubstitutionModelInterface> pCodonModel,
      std::unique_ptr<CodonFrequencySetInterface> pFreq,
      std::shared_ptr<const GeneticCode> pgencode);

  CodonSameAARateSubstitutionModel(
      const CodonSameAARateSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    AbstractSubstitutionModel(model),
    pAAmodel_ (model.pAAmodel_->clone()),
    pCodonModel_ (model.pCodonModel_->clone()),
    pFreq_ (model.pFreq_ ? model.pFreq_->clone() : nullptr),
    pgencode_ (model.pgencode_),
    X_ (model.X_),
    phi_ (model.phi_)
  {
    compute_();
    updateMatrices_();
  }

  CodonSameAARateSubstitutionModel& operator=(
      const CodonSameAARateSubstitutionModel& model)
  {
    AbstractSubstitutionModel::operator=(model);

    pAAmodel_.reset(model.pAAmodel_->clone());
    pCodonModel_.reset(model.pCodonModel_->clone());
    pFreq_.reset(model.pFreq_ ? model.pFreq_->clone() : nullptr);

    pgencode_ = model.pgencode_;

    X_ = model.X_;
    phi_  = model.phi_;

    return *this;
  }

  CodonSameAARateSubstitutionModel* clone() const override
  {
    return new CodonSameAARateSubstitutionModel(*this);
  }

  virtual ~CodonSameAARateSubstitutionModel() {}

public:
  std::string getName() const override
  {
    return "SameAARate";
  }

  const ProteinSubstitutionModelInterface& proteinModel() const
  {
    return *pAAmodel_;
  }

  const CodonSubstitutionModelInterface& codonModel() const
  {
    return *pCodonModel_;
  }

  void fireParameterChanged(const ParameterList& parameters) override;

  void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix);

    pAAmodel_->setNamespace(prefix + pAAmodel_->getNamespace());
    pCodonModel_->setNamespace(prefix + pCodonModel_->getNamespace());

    if (pFreq_)
      pFreq_->setNamespace(prefix + pFreq_->getNamespace());
  }

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    return pCodonModel_->codonFrequencySet();
  }

  bool hasCodonFrequencySet() const override
  {
    return pCodonModel_->hasCodonFrequencySet();
  }

  std::shared_ptr<const GeneticCode> getGeneticCode() const override
  {
    return pCodonModel_->getGeneticCode();
  }

  void setFreq(std::map<int, double>& frequencies) override
  {
    dynamic_cast<CoreCodonSubstitutionModelInterface&>(*pCodonModel_).setFreq(frequencies);
    matchParametersValues(pCodonModel_->getParameters());
  }

  double getCodonsMulRate(size_t i, size_t j) const override
  {
    return 1.;
  }

private:
  void compute_();
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_CODONSAMEAARATESUBSTITUTIONMODEL_H
