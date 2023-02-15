//
// File: CodonSameAARateSubstitutionModel.h
// Authors:
//   Laurent Gueguen
// Created: vendredi 30 octobre 2020, ÃÂ  17h 44
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
