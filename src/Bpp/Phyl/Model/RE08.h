// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_RE08_H
#define BPP_PHYL_MODEL_RE08_H


#include "AbstractSubstitutionModel.h"
#include "Codon/CodonSubstitutionModel.h"
#include "Nucleotide/NucleotideSubstitutionModel.h"
#include "Protein/ProteinSubstitutionModel.h"
#include "SubstitutionModel.h"

namespace bpp
{
/**
 * @brief The Rivas-Eddy substitution model with gap characters.
 *
 * This model expends any reversible substitution model with gaps as an additional state.
 * Although the conditionnal subtitution process is reversible, the insertion/deletion process
 * needs not be. The model hence adds two parameters for insertion and deletions, @f$\lambda@f$ and @f$\mu@f$.
 * If we note @f$Q@f$ the (simple) transition matrix (= Markov generator) and @f$Q^\epsilon@f$ the extended one, we have:
 * @f[
 * Q^\epsilon =
 * \left(
 * \begin{array}{ccc|c}
 * & & & \mu \\
 * & \rule{0cm}{1cm}\rule{1cm}{0cm}Q-\mu\delta_{ij} & & \vdots \\
 * & & & \mu \\
 * \hline
 * \lambda\pi_1 & \ldots & \lambda\pi_n & -\lambda \\
 * \end{array}
 * \right)
 * @f]
 * where @f$n@f$ is the number of states of the simple model (in most case equal to the size of the alphabet) and @f$(\pi_1,\ldots,\pi_n)@f$ is the vector of equilibrium frequencies of the conditional model.
 * @f$\delta_{ij}@f$ is 1 if i=j, 0 otherwise.
 * Note that in the original paper @f$Q@f$ is noted as @f$R@f$, and @f$Q_t@f$ is used for the probability matrix, which is referred here as @f$P^\epsilon(t)@f$ for consistency with the documentation of other models.
 *
 * The extended Markov model is reversible, and the equilibrium frequencies are
 * @f[
 * \pi^\epsilon = \left( \pi \cdot \frac{\lambda}{\lambda + \mu}, \frac{\mu}{\lambda + \mu}\right).
 * @f]
 * The corresponding exchangeability matrix is:
 * @f[
 * S^\epsilon =
 * \left(
 * \begin{array}{ccc|c}
 * & & & \lambda + \mu \\
 * & \rule{0cm}{1cm}\rule{1cm}{0cm}(S - \frac{\mu\delta_{ij}}{\pi_i})\frac{\lambda+\mu}{\lambda} & & \vdots \\
 * & & & \lambda + \mu \\
 * \hline
 * \lambda + \mu & \ldots & \lambda + \mu & - (\lambda + \mu) \\
 * \end{array}
 * \right)
 * @f]
 * The eigen values and vectors are computed numerically, but the transition probabilities are computed analytically from the simple substitution model, together with the first and second order derivatives according to time.
 *
 * Please note that the original Rivas and Eddy method uses this substitution model with a modification of the Felsenstein algorithm.
 *
 * Reference:
 * - Rivas E and Eddy SR (2008), _Probabilistic Phylogenetic Inference with Insertions and Deletions_, 4(9):e1000172, in _PLoS Computational Biology_.
 */
class RE08 :
  public AbstractReversibleSubstitutionModel
{
private:
  std::unique_ptr<ReversibleSubstitutionModelInterface> simpleModel_;
  RowMatrix<double> simpleGenerator_;
  RowMatrix<double> simpleExchangeabilities_;
  mutable double exp_;
  mutable RowMatrix<double> p_;
  double lambda_;
  double mu_;
  std::string nestedPrefix_;

public:
  /**
   * @brief Build a new Rivas & Eddy model from a standard substitution model.
   *
   * The alphabet and number of states for the extended model will be derived from the simple one.
   *
   * @param simpleModel The simple model to use to build the extended one.
   * THE RE08 class will own the simple one, meaning that it will be destroyed together with the RE08 instance, and cloned when cloning the RE08 instance.
   * To prevent the original simple model to be destroyed, you should make a copy of it before creating the RE08 instance.
   * @param lambda Insertion rate.
   * @param mu     Deletion rate.
   */
  RE08(
      std::unique_ptr<ReversibleSubstitutionModelInterface> simpleModel,
      double lambda = 0.1,
      double mu = 0.1);

  RE08(const RE08& model) :
    AbstractParameterAliasable(model),
    AbstractReversibleSubstitutionModel(model),
    simpleModel_(model.simpleModel_->clone()),
    simpleGenerator_(model.simpleGenerator_),
    simpleExchangeabilities_(model.simpleExchangeabilities_),
    exp_(model.exp_),
    p_(model.p_),
    lambda_(model.lambda_),
    mu_(model.mu_),
    nestedPrefix_(model.nestedPrefix_)
  {}

  RE08& operator=(const RE08& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractReversibleSubstitutionModel::operator=(model);
    simpleModel_.reset(model.simpleModel_->clone());
    simpleGenerator_         = model.simpleGenerator_;
    simpleExchangeabilities_ = model.simpleExchangeabilities_;
    exp_                     = model.exp_;
    p_                       = model.p_;
    lambda_                  = model.lambda_;
    mu_                      = model.mu_;
    nestedPrefix_            = model.nestedPrefix_;
    return *this;
  }

  virtual ~RE08() {}

  RE08* clone() const override { return new RE08(*this); }

public:
  double Pij_t    (size_t i, size_t j, double d) const override;
  double dPij_dt  (size_t i, size_t j, double d) const override;
  double d2Pij_dt2(size_t i, size_t j, double d) const override;
  const Matrix<double>& getPij_t    (double d) const override;
  const Matrix<double>& getdPij_dt  (double d) const override;
  const Matrix<double>& getd2Pij_dt2(double d) const override;

  std::string getName() const override { return "RE08"; }

  /**
   * @brief This method is forwarded to the simple model.
   *
   * @param data The data to be passed to the simple model (gaps will be ignored).
   * @param pseudoCount A (typically small) value to add to each count to avoid 0 estimates.
   */
  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override
  {
    simpleModel_->setFreqFromData(data, pseudoCount);
  }

  void setFreq(std::map<int, double>& frequencies) override
  {
    simpleModel_->setFreq(frequencies);
  }

  void fireParameterChanged(const ParameterList& parameters) override
  {
    simpleModel_->matchParametersValues(parameters);
    lambda_ = getParameter_(0).getValue();
    mu_     = getParameter_(1).getValue();
    updateMatrices_();
  }

  double getInitValue(size_t i, int state) const override;

  void setNamespace(const std::string& prefix) override;

  const ReversibleSubstitutionModelInterface& nestedModel() const
  { 
    return *simpleModel_; 
  }

protected:
  
  void updateMatrices_() override;
  
  ReversibleSubstitutionModelInterface& nestedModel_() { return *simpleModel_; }
};


/**
 * @brief This is a wrapper class of RE08 for nucleotide substitution models.
 */
class RE08Nucleotide :
  public RE08,
  public virtual NucleotideSubstitutionModelInterface
{
public:
  /**
   * @brief Build a new Rivas & Eddy model from a nucleotide
   * substitution model.
   *
   * The alphabet and number of states for the extended model will
   * be derived from the simple one.
   *
   * @param nucleotideModel The simple model to use to build the
   * extended one. THE RE08 class will own the simple one, meaning
   * that it will be destroyed together with the RE08 instance, and
   * cloned when cloning the RE08 instance. To prevent the original
   * simple model to be destroyed, you should make a copy of it
   * before creating the RE08 instance.
   * @param lambda Insertion rate.
   * @param mu     Deletion rate.
   */
  RE08Nucleotide(
      std::unique_ptr<NucleotideReversibleSubstitutionModelInterface> nucleotideModel,
      double lambda = 0.1,
      double mu = 0.1) :
    AbstractParameterAliasable("RE08."),
    RE08(std::move(nucleotideModel), lambda, mu) {}

  virtual ~RE08Nucleotide() {}

  RE08Nucleotide* clone() const override { return new RE08Nucleotide(*this); }

public:
  std::shared_ptr<const NucleicAlphabet> getNucleicAlphabet() const override
  {
    return std::dynamic_pointer_cast<const NucleicAlphabet>(alphabet_);
  }

  const NucleotideReversibleSubstitutionModelInterface& nestedModel() const
  {
    return dynamic_cast<const NucleotideReversibleSubstitutionModelInterface&>(RE08::nestedModel());
  }

  size_t getNumberOfStates() const override { return RE08::getNumberOfStates(); }
};

/**
 * @brief This is a wrapper class of RE08 for protein substitution models.
 */
class RE08Protein :
  public RE08,
  public virtual ProteinSubstitutionModelInterface
{
public:
  /**
   * @brief Build a new Rivas & Eddy model from a protein substitution model.
   *
   * The alphabet and number of states for the extended model will be derived from the simple one.
   *
   * @param proteinModel The simple model to use to build the extended one.
   * THE RE08 class will own the simple one, meaning that it will be destroyed together with the RE08 instance, and cloned when cloning the RE08 instance.
   * To prevent the original simple model to be destroyed, you should make a copy of it before creating the RE08 instance.
   * @param lambda Insertion rate.
   * @param mu     Deletion rate.
   */
  RE08Protein(
      std::unique_ptr<ProteinReversibleSubstitutionModelInterface> proteinModel,
      double lambda = 0.1,
      double mu = 0.1) :
    AbstractParameterAliasable("RE08."),
    RE08(std::move(proteinModel), lambda, mu) {}

  virtual ~RE08Protein() {}

  RE08Protein* clone() const override{ return new RE08Protein(*this); }

public:
  std::shared_ptr<const ProteicAlphabet> getProteicAlphabet() const override
  {
    return std::dynamic_pointer_cast<const ProteicAlphabet>(alphabet_);
  }
  
  const ProteinReversibleSubstitutionModelInterface& nestedModel() const
  {
    return dynamic_cast<const ProteinReversibleSubstitutionModelInterface&>(RE08::nestedModel());
  }

  size_t getNumberOfStates() const override { return RE08::getNumberOfStates(); }
};


/**
 * @brief This is a wrapper class of RE08 for codon substitution models.
 */
class RE08Codon :
  public RE08,
  public virtual CodonSubstitutionModelInterface
{
public:
  /**
   * @brief Build a new Rivas & Eddy model from a codon substitution model.
   *
   * The alphabet and number of states for the extended model will be derived from the simple one.
   *
   * @param codonModel The simple model to use to build the extended one.
   * THE RE08 class will own the simple one, meaning that it will be destroyed together with the RE08 instance, and cloned when cloning the RE08 instance.
   * To prevent the original simple model to be destroyed, you should make a copy of it before creating the RE08 instance.
   * @param lambda Insertion rate.
   * @param mu     Deletion rate.
   */
  RE08Codon(
      std::unique_ptr<CodonReversibleSubstitutionModelInterface> codonModel,
      double lambda = 0.1,
      double mu = 0.1) :
    AbstractParameterAliasable("RE08."),
    RE08(std::move(codonModel), lambda, mu) {}

  virtual ~RE08Codon() {}

  RE08Codon* clone() const override { return new RE08Codon(*this); }

public:
  const CodonReversibleSubstitutionModelInterface& nestedModel() const
  {
    return dynamic_cast<const CodonReversibleSubstitutionModelInterface&>(RE08::nestedModel());
  }

  std::shared_ptr<const GeneticCode> getGeneticCode () const override
  {
    return nestedModel().getGeneticCode();
  }

  double getCodonsMulRate(size_t i, size_t j) const override
  {
    return nestedModel().getCodonsMulRate(i, j);
  }

  size_t getNumberOfStates() const override { return RE08::getNumberOfStates(); }

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    return dynamic_cast<const CodonFrequencySetInterface&>(RE08::nestedModel().frequencySet());
  }

  bool hasCodonFrequencySet() const override
  {
    return true;
  }

  void setFreq(std::map<int, double>& frequencies) override
  {
    RE08::nestedModel_().setFreq(frequencies);
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_RE08_H
