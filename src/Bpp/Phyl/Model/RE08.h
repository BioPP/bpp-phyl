//
// File: RE08.h
// Created by: Julien Dutheil
// Created on: Mon Dec 29 10:15 2008
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use, 
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info". 

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability. 

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or 
  data to be ensured and,  more generally, to use and operate it in the 
  same conditions as regards security. 

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _RE08_H_
#define _RE08_H_

#include "SubstitutionModel.h"
#include "AbstractSubstitutionModel.h"
#include "Nucleotide/NucleotideSubstitutionModel.h"
#include "Protein/ProteinSubstitutionModel.h"
#include "Codon/CodonSubstitutionModel.h"

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
  class RE08:
    public AbstractReversibleSubstitutionModel
  {
  private:
    std::auto_ptr<ReversibleSubstitutionModel> simpleModel_;
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
    RE08(ReversibleSubstitutionModel* simpleModel, double lambda = 0.1, double mu = 0.1);

    RE08(const RE08& model):
      AbstractParameterAliasable(model),
      //AbstractSubstitutionModel(model),
      AbstractReversibleSubstitutionModel(model),
      simpleModel_(dynamic_cast<ReversibleSubstitutionModel*>(model.simpleModel_->clone())),
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
      AbstractSubstitutionModel::operator=(model);
      AbstractReversibleSubstitutionModel::operator=(model);
      simpleModel_.reset(dynamic_cast<ReversibleSubstitutionModel*>(model.simpleModel_->clone()));
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

    RE08* clone() const { return new RE08(*this); }

  public:
	
    double Pij_t    (size_t i, size_t j, double d) const;
    double dPij_dt  (size_t i, size_t j, double d) const;
    double d2Pij_dt2(size_t i, size_t j, double d) const;
    const Matrix<double>& getPij_t    (double d) const;
    const Matrix<double>& getdPij_dt  (double d) const;
    const Matrix<double>& getd2Pij_dt2(double d) const;

    std::string getName() const { return "RE08"; }

    /**
     * @brief This method is forwarded to the simple model.
     *
     * @param data The data to be passed to the simple model (gaps will be ignored).
     * @param pseudoCount A (typically small) value to add to each count to avoid 0 estimates.
     */
    void setFreqFromData(const SequenceContainer& data, double pseudoCount = 0)
    {
      simpleModel_->setFreqFromData(data, pseudoCount);
    }
	
    void fireParameterChanged(const ParameterList& parameters)
    {
      AbstractParameterAliasable::fireParameterChanged(parameters);      
      simpleModel_->matchParametersValues(parameters);
      lambda_ = getParameter_(0).getValue();
      mu_     = getParameter_(1).getValue();
      updateMatrices();
    }

    size_t getNumberOfStates() const { return size_; }

    double getInitValue(size_t i, int state) const throw (IndexOutOfBoundsException, BadIntException);
  
    void setNamespace(const std::string& prefix);

    const ReversibleSubstitutionModel* getNestedModel() const { return simpleModel_.get(); }

  protected:

    void updateMatrices();
  };

 
  
  /**
   * @brief This is a wrapper class of RE08 for nucleotide substitution models.
   */
  class RE08Nucleotide:
    public RE08,
    public virtual NucleotideSubstitutionModel
  {
  public:
    /**
     * @brief Build a new Rivas & Eddy model from a nuclotide substitution model.
     * 
     * The alphabet and number of states for the extended model will be derived from the simple one.
     *
     * @param nucleotideModel The simple model to use to build the extended one.
     * THE RE08 class will own the simple one, meaning that it will be destroyed together with the RE08 instance, and cloned when cloning the RE08 instance.
     * To prevent the original simple model to be destroyed, you
     * should make a copy of it before creating the RE08 instance.
     *
     * @param lambda Insertion rate.
     * @param mu     Deletion rate.
     */

    RE08Nucleotide(NucleotideReversibleSubstitutionModel* nucleotideModel, double lambda = 0.1, double mu = 0.1):
      AbstractParameterAliasable("RE08."),
      RE08(nucleotideModel, lambda, mu) {};

    virtual ~RE08Nucleotide() {}

    RE08Nucleotide* clone() const { return new RE08Nucleotide(*this); }

  public:
    const NucleicAlphabet* getAlphabet() const { return dynamic_cast<const NucleicAlphabet*>(alphabet_); }
    const NucleotideReversibleSubstitutionModel* getNestedModel() const {
      return dynamic_cast<const NucleotideReversibleSubstitutionModel*>(RE08::getNestedModel());
    }
    
    size_t getNumberOfStates() const { return RE08::getNumberOfStates(); }

  };
 
  /**
   * @brief This is a wrapper class of RE08 for protein substitution models.
   */
  class RE08Protein:
    public RE08,
    public virtual ProteinSubstitutionModel
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
    RE08Protein(ProteinReversibleSubstitutionModel* proteinModel, double lambda = 0.1, double mu = 0.1):
      AbstractParameterAliasable("RE08."),
      RE08(proteinModel, lambda, mu) {};

    virtual ~RE08Protein() {}

    RE08Protein* clone() const { return new RE08Protein(*this); }

  public:
    const ProteicAlphabet* getAlphabet() const { return dynamic_cast<const ProteicAlphabet*>(alphabet_); }
    const ProteinReversibleSubstitutionModel* getNestedModel() const {
      return dynamic_cast<const ProteinReversibleSubstitutionModel*>(RE08::getNestedModel());
    }

    size_t getNumberOfStates() const { return RE08::getNumberOfStates(); }
  };



  /**
   * @brief This is a wrapper class of RE08 for codon substitution models.
   */
  class RE08Codon:
    public RE08,
    public virtual CodonSubstitutionModel
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
    RE08Codon(CodonReversibleSubstitutionModel* codonModel, double lambda = 0.1, double mu = 0.1):
      AbstractParameterAliasable("RE08."),
      RE08(codonModel, lambda, mu) {};

    virtual ~RE08Codon() {}

    RE08Codon* clone() const { return new RE08Codon(*this); }

  public:
    const CodonReversibleSubstitutionModel* getNestedModel() const {
      return dynamic_cast<const CodonReversibleSubstitutionModel*>(RE08::getNestedModel());
    }

    const GeneticCode* getGeneticCode () const {
      return getNestedModel()->getGeneticCode();
    }

    double	getCodonsMulRate (size_t i, size_t j) const { 
      return getNestedModel()->getCodonsMulRate(i, j);
    }

    size_t getNumberOfStates() const { return RE08::getNumberOfStates(); }
  };


} //end of namespace bpp.

#endif	//_RE08_H_

