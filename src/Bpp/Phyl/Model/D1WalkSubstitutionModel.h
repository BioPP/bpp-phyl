// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_D1WALKSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_D1WALKSUBSTITUTIONMODEL_H

#include <Bpp/Seq/Alphabet/IntegerAlphabet.h>

#include "AbstractSubstitutionModel.h"


namespace bpp
{
/**
 * @brief The D1Walk substitution model for Integer alphabet [0;N-1].
 *        In this model, substitutions are possible between adjacent
 *        states only. The model is defined through an equilibrium
 *        distribution, and rates towards a state are proportional to
 *        its frequency at equilibrium.
 *
 *  The generator is the dot product between a reversible matrix \f[ S
 *  \f] and an equilibrium vector \f[ \pi \f] with
 *
 * \f[
 * S = \begin{pmatrix}
 *  \cdots & 1    & 0 & \ldots & \ldots  & \ldots & \ldots \\
 *  1  & \cdots & 1 & 0      & \\ldots  & \ldots & \ldots \\
 * 0 & 1  & \cdots & 1 & 0 & \ldots & \ldots \\
 * \ldots & 0 & 1  & \cdots & 1 & 0 & \ldots \\
 * \vdots & \ddots & \ddots  & \ddots & \ddots & \ddots & \vdots \\
 * \ldots &\ldots & 0 & 1  & \cdots  & 1 & 0  \\
 * \ldots & \ldots &\ldots & 0 & 1  & \cdots  & 1  \\
 * \ldots & \ldots & \ldots &\ldots & 0 & 1  & \cdots \\
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = \left(\pi_1, \pi_2, \dots, \pi_N)
 * \f]
 *
 * with a normalization factor \f$P\f$:
 *
 * \f[
 * Q = \frac{1}{P} S . \pi = \frac{1}{P}\begin{pmatrix}
 *  -\pi_2 & \pi_2    & 0 & \ldots & \ldots  & \ldots & \ldots \\
 *  \pi_1  & -(\pi_1+\pi_3) & \pi_3 & 0      & \\ldots  & \ldots & \ldots \\
 * 0 & \pi_2  & -(\pi_2+\pi_4) & \pi_4 & 0 & \ldots & \ldots \\
 * \ldots & 0 & \pi_3  & -(\pi_3+\pi_5) & \pi_5 & 0 & \ldots \\
 * \vdots & \ddots & \ddots  & \ddots & \ddots & \ddots & \vdots \\
 * \ldots &\ldots & 0 & \pi_{N-3}  & -(\pi_{N-3}+\pi_{N-1})\cdots  & \pi_{N-1} & 0  \\
 * \ldots & \ldots &\ldots & 0 & \pi_{N-2}  & -(\pi_{N-2}+\pi_N)\cdots  & \pi_{N}  \\
 * \ldots & \ldots & \ldots &\ldots & 0 & \pi_{N-1}  & -\pi_{N-1} \\
 * \end{pmatrix}
 * \f]
 *
 * with \f$P = 2 \sum_i^{N-1} \pi_i . \pi_{i+1} \f$
 *
 * This model hence include N parameters for frequencies, and N-1 free
 * parameters that are denoted \c "theta1", \c "theta2", ...
 *
 * The parametrization depends on the method used.
 * Default method is 1 (ie global ratio).
 * @see Simplex
 *
 * In addition, an optional rate_ factor defines the mean rate of the
 * model.
 *
 * The generator of this model is diagonalized numericaly.
 * See AbstractSubstitutionModel for details of how the probabilities are computed.
 *
 */

class D1WalkSubstitutionModel :
  public AbstractReversibleSubstitutionModel
{
private:
  /**
   * @brief The Equilibrium Frequency Set
   *
   **/

  std::shared_ptr<FullFrequencySet> freqSet_;

public:
  /**
   * @brief Build a D1Walk model on a given IntegerAlphabet.
   *
   * @param alpha the IntegerAlphabet.
   * @param method the method used to compute the equilibrium
   * frequencies from the theta parameters.
   */


  D1WalkSubstitutionModel(std::shared_ptr<const IntegerAlphabet> alpha, unsigned short method = 1);

  D1WalkSubstitutionModel(const D1WalkSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    AbstractReversibleSubstitutionModel(model),
    freqSet_(model.freqSet_->clone())
  {}

  D1WalkSubstitutionModel& operator=(const D1WalkSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractReversibleSubstitutionModel::operator=(model);
    freqSet_ = std::shared_ptr<FullFrequencySet>(model.freqSet_->clone());
    return *this;
  }

  virtual ~D1WalkSubstitutionModel() {}

  D1WalkSubstitutionModel* clone() const override { return new D1WalkSubstitutionModel(*this); }

public:
  std::string getName() const override
  {
    return "D1Walk";
  }

  void fireParameterChanged(const ParameterList& parameters) override
  {
    freqSet_->matchParametersValues(parameters);
    AbstractReversibleSubstitutionModel::fireParameterChanged(parameters);
  }

  const FrequencySetInterface& frequencySet() const override
  {
    return *freqSet_;
  }

  void setFreq(std::map<int, double>& freq) override;

protected:
  /**
   * In the case of the model of Jukes & Cantor, this method is useless since
   * the generator is fixed! No matrice can be changed... This method is only
   * used in the constructor of the class.
   */
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_D1WALKSUBSTITUTIONMODEL_H
