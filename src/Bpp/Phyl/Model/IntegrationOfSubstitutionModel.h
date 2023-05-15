//
// File: IntegrationOfSubstitutionModel.h
// Authors:
//   Laurent Gueguen
// Created: samedi 24 octobre 2015, Ã  18h 28
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_MODEL_INTEGRATIONOFSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_INTEGRATIONOFSUBSTITUTIONMODEL_H


#include "AbstractFromSubstitutionModelTransitionModel.h"

namespace bpp
{
/**
 * @brief From a model, compute the expected transition probabilities:
 *
 * if the branch length @f$l@f$ is a random variable that follows a
 * distribution of density @f$\delta@f$, then transition probability
 * @f$P(l)@f$ is also a random variable which @f$m^\text{th}@f$ moment is:
 * @f$E(P^m)=\int_l e^{mlQ} \delta(l) dl@f$
 *
 * A convenient modeling of the variability of the branch length would
 * be, with @f$\forall i, \zeta_i > 0@f$ and @f$\zeta_1 + ... + \zeta_n = 1@f$,
 * @f$\delta_t = \frac 1n \sum_i
 * \Gamma(k,\frac{n \zeta_i t}{k})@f$, with @f$k \in \mb N^*@f$
 * with expectation @f$E(T)= t@f$ and variance
 * @f$var(T)=\frac{t^2}{k} \sum_i \zeta_i^2 @f$.
 *
 * Then, from the fact that @f$E(P)@f$ is like the moment-generating
 * function for a matrix:
 * 
 * 
 * @f$E(P) = \frac 1n \sum_i (I-\frac{\zeta_i n t}{k}Q)^{-k}@f$
 * 
 *
 * It has the same parameters as the SubModel, plus \c "thetas" of the
 * Simplex encoding the @f$\zeta@f$s.
 *
 * Its construction depends on value of integer \c "k" (which is not a
 * parameter) and the number of @f$\zeta@f$, \c "n".
 *
 */

class IntegrationOfSubstitutionModel :
  public AbstractFromSubstitutionModelTransitionModel
{
private:
  Simplex svar_;
  
  /***
   * @brief vector of Zetas.
   *
   **/

  Vdouble vZeta_;

  /****
   * @brief k
   *
   **/
  
  uint k_;
  
  /****
   * @brief Store intermediate matrices for derivatives
   *  Matrix (I-zeta_i/{k * Q * t )^{-k}
   **/
  
  mutable std::vector<RowMatrix<double>>  vMatInv_;

  mutable double time_; // compare time_ with t to ensure the right time is used

  mutable RowMatrix<double> tmp_, tmp2_, tmp3_;

  void computeInv_(double t) const;

public:
  IntegrationOfSubstitutionModel(std::unique_ptr<SubstitutionModelInterface>& originalModel, uint k, uint n);

  IntegrationOfSubstitutionModel(const IntegrationOfSubstitutionModel& fmsm) :
    AbstractParameterAliasable(fmsm),
    AbstractWrappedModel(fmsm),
    AbstractWrappedTransitionModel(fmsm),
    AbstractFromSubstitutionModelTransitionModel(fmsm),
    svar_(fmsm.svar_),
    vZeta_(fmsm.vZeta_),
    k_(fmsm.k_),
    vMatInv_(),
    time_(-1),
    tmp_(model().getNumberOfStates(), model().getNumberOfStates()),
    tmp2_(model().getNumberOfStates(), model().getNumberOfStates()),
    tmp3_(model().getNumberOfStates(), model().getNumberOfStates())
  {
    for (size_t i=0; i<fmsm.vMatInv_.size();i++)
      vMatInv_.push_back(RowMatrix<double>(model().getNumberOfStates(), model().getNumberOfStates()));
  }

  IntegrationOfSubstitutionModel& operator=(const IntegrationOfSubstitutionModel& fmsm)
  {
    AbstractFromSubstitutionModelTransitionModel::operator=(fmsm);
    svar_ = fmsm.svar_;
    vZeta_ = fmsm.vZeta_;
    k_ = fmsm.k_;

    vMatInv_.resize(fmsm.vMatInv_.size());
    for (size_t i=0; i<fmsm.vMatInv_.size();i++)
      vMatInv_.push_back(RowMatrix<double>(model().getNumberOfStates(), model().getNumberOfStates()));

    time_=-1;
    tmp_.resize(model().getNumberOfStates(), model().getNumberOfStates());
    tmp2_.resize(model().getNumberOfStates(), model().getNumberOfStates());
    tmp3_.resize(model().getNumberOfStates(), model().getNumberOfStates());
    return *this;
  }

  uint k() const { return k_;}

  size_t numberOfGammas() const { return vZeta_.size();}

  ~IntegrationOfSubstitutionModel() {}

  IntegrationOfSubstitutionModel* clone() const { return new IntegrationOfSubstitutionModel(*this); }

public:
  double Pij_t    (size_t i, size_t j, double t) const;
  double dPij_dt  (size_t i, size_t j, double t) const;
  double d2Pij_dt2(size_t i, size_t j, double t) const;

  const Matrix<double>& getPij_t(double t) const;

  const Matrix<double>& getdPij_dt(double t) const;

  const Matrix<double>& getd2Pij_dt2(double t) const;

  double getInitValue(size_t i, int state) const override
  {
    return transitionModel().getInitValue(i, state);
  }

  double getRate() const override
  {
    return transitionModel().getRate();
  }

  void setRate(double rate) override
  {
    return transitionModel_().setRate(rate);
  }

  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override
  {
    std::map<int, double> freqs;
    SequenceContainerTools::getFrequencies(data, freqs, pseudoCount);
    // Re-compute generator and eigen values:
    transitionModel_().setFreq(freqs);
  }

  double freq(size_t i) const override { return transitionModel().freq(i); }

  void setFreq(std::map<int, double>& frequencies) override
  {
    transitionModel_().setFreq(frequencies);
  }

  const Vdouble& getFrequencies() const override
  {
    return transitionModel().getFrequencies();
  }

  bool computeFrequencies() const override
  {
    return transitionModel().computeFrequencies();
  }

  /**
   * @return Set if equilibrium frequencies should be computed from
   * the generator
   */
  void computeFrequencies(bool yn) override
  {
    transitionModel_().computeFrequencies(yn);
  }

  /**
   * @}
   */

protected:

  Vdouble& getFrequencies_() override
  {
    return transitionModel_().getFrequencies_();
  }

  
public:
  
  /*
   * @brief The integrative model has same equilibrium frequencies as
   * SubModel.
   *
   */
    
  std::string getName() const
  {
    return "Integrate." + model().getName();
  }

  void fireParameterChanged(const ParameterList& parameters);

  /*
   * @}
   *
   */
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ONECHANGETRANSITIONMODEL_H
