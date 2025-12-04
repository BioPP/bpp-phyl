// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_POMO_H
#define BPP_PHYL_MODEL_POMO_H


#include "FrequencySet/FrequencySet.h"
#include "AbstractSubstitutionModel.h"

/**
 * @brief POMO-like Model, based on AllelicAlphabet.
 *
 *
 *
 * In the alphabet from which the AllelicAlphabet is built, this model
 * needs a given mutation model with generator (aka SubstitutionModel
 * in Bio++), and allele selection coefficients from a finesse vector
 * (aka FrequencySet in Bio++).
 *
 * If we denote @f$F@f$ this fitness, the drift term towards the
 * fixation of allele @f$i@f$ is proportional @f$F_i@f$, and follows
 * the Moran law.
 *
 * We set the generator @f$Q@f$ of the model on set of states
 * @f$\big\{\{na_i,(N-n)a_j\}, \forall i<j \& \forall n \in [1;N-1]
 * \big\} \cup \{Na_i, \forall i \}@f$,\\ with @f$\mu_{ab} = \rho_{ab}
 * \pi_a@f$ the reversible mutation rates, with @f$\pi_{a}@f$ the
 * equilibrium distribution of this mutation model, and
 * @f$\rho_{ab}=\rho_{ba}@f$, and @f$\phi@f$ the fitnesses of the
 * states:
 *
 *@f$
 * Q_{\{ua_i,(N-u)a_j\} \rightarrow \{va_i,(N-v)a_j\}} = \left\{
 *\begin{array}{cl}
 * N \mu_{a_ia_j} = N \rho_{a_ia_j}\pi_{a_j} & \text{ if } u=N, v=N-1\\
 * N \mu_{a_ja_i} = N \rho_{a_ja_i}\pi_{a_i} & \text{ if } u=0, v=1\\
 * \frac{n(N-n)}{n\phi_{a_i}+(N-n)\phi_{a_j}}\phi_{a_i} & \text{ if } u=n, v=n+1, 0<n<N \\
 * \frac{n(N-n)}{n\phi_{a_i}+(N-n)\phi_{a_j}}\phi_{a_j} & \text{ if } u=n, v=n-1, 0<n<N \\
 * 0 & \text{ if } |u-v|>1
 * \end{array} \right. @f$
 *
 * (with @f$\{Na_i,0a_j\}:=Na_i@f$)
 *
 * The equilibrium distribution @f$\psi@f$ is:
 *
 * @f$
 * \psi_x = \frac 1K \times \left\{
 * \begin{array}{cl}
 * \pi_{a_i} \phi_{a_i}^{N-1} & \text{ if } x = Na_i\\
 * \pi_{a_i} \mu_{a_ia_j} \phi_{a_i}^{n-1}\phi_{a_j}^{N-n-1}(n\phi_{a_i}+(N-n)\phi_{a_j})\frac{N}{n(N-n)}
 * & \text{ if } x = \{na_i,(N-n)a_j\} \\
 * \end{array}
 * \right.
 * @f$
 *
 * with @f$ K = \sum_{i} \pi_{a_i} \phi_{a_i}^{N-1} + \sum_{i \neq j}
 * \pi_{a_i} \pi_{a_j} \rho_{a_ia_j} \sum_{n=1}^{N-1} \phi_{a_i}^{N-n}
 * \phi_{a_j}^{n-1}\frac{N}{n} @f$
 *
 * The generator should be normalized to ensure one substitution per
 * unit of time. See document online:
 *
 * https://www.overleaf.com/project/639760f5aacdcfbde3ced551
 *
 * See:
 *
 * De Maio, N., Schrempf, D., and Kosiol, C. (2015). PoMo: An allele
 * frequency-based approach for 279 species tree estimation. Systematic
 * Biology, 64(6):1018–1031.
 *
 * Borges R, Szöllősi GJ, Kosiol C. Quantifying GC-Biased Gene Conversion
 * in Great Ape Genomes Using Polymorphism-Aware Models.
 * Genetics. 2019 Aug;212(4):1321-1336.
 *
 * Borges, R. and Kosiol, C. (2020). Consistency and identifiability
 * of the polymorphism-aware phylo-genetic models. Journal of
 * Theoretical Biology, 486:110074.
 */
#include <Bpp/Seq/Alphabet/AllelicAlphabet.h>

namespace bpp
{
class POMO :
  public AbstractSubstitutionModel
{
private:
  unsigned int nbAlleles_;

  std::unique_ptr<SubstitutionModelInterface> pmodel_;
  std::unique_ptr<FrequencySetInterface> pfitness_;

public:
  /**
   * @brief Build a POMO instance
   */
  POMO(std::shared_ptr<const AllelicAlphabet> allAlph,
      std::unique_ptr<SubstitutionModelInterface> pmodel,
      std::unique_ptr<FrequencySetInterface> pfitness);

  POMO(const POMO& model) :
    AbstractParameterAliasable(model),
    AbstractSubstitutionModel(model),
    nbAlleles_(model.nbAlleles_),
    pmodel_(model.pmodel_->clone()),
    pfitness_(model.pfitness_ ? model.pfitness_->clone() : nullptr)
  {}

  POMO& operator=(const POMO& model)
  {
    AbstractParameterAliasable::operator=(model);
    nbAlleles_ = model.nbAlleles_;
    pmodel_.reset(model.pmodel_->clone());
    pfitness_.reset(model.pfitness_ ? model.pfitness_->clone() : nullptr);
    return *this;
  }

  POMO* clone() const override
  {
    return new POMO(*this);
  }

  void fireParameterChanged(const ParameterList& parameters) override;

  void setFreq(std::map<int, double>& frequencies) override;

  void setNamespace(const std::string& prefix) override
  {
    AbstractParameterAliasable::setNamespace(prefix);
    pmodel_->setNamespace(prefix + pmodel_->getNamespace());
    if (pfitness_)
      pfitness_->setNamespace(prefix + pfitness_->getNamespace());
  }

  std::shared_ptr<const AllelicAlphabet> getAllelicAlphabet() const
  {
    return std::dynamic_pointer_cast<const AllelicAlphabet>(alphabet_);
  }

  const AllelicAlphabet& allelicAlphabet() const
  {
    return dynamic_cast<const AllelicAlphabet&>(*alphabet_);
  }

  unsigned int getNbAlleles() const
  {
    return nbAlleles_;
  }

  bool hasFitness() const
  {
    return pfitness_ != nullptr;
  }

  const FrequencySetInterface& fitness() const
  {
    return *pfitness_;
  }

  const SubstitutionModelInterface& mutationModel() const
  {
    return *pmodel_;
  }

  std::string getName() const override
  {
    return "POMO";
  }

protected:
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_POMO_H
