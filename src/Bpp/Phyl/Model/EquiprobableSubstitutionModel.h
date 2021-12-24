//
// File: EquiprobableSubstitutionModel.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: mercredi 13 juillet 2016, ÃÂ  08h 49
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

#ifndef BPP_PHYL_MODEL_EQUIPROBABLESUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_EQUIPROBABLESUBSTITUTIONMODEL_H

#include <Bpp/Seq/Alphabet/Alphabet.h>

#include "AbstractSubstitutionModel.h"
#include "FrequencySet/FrequencySet.h"

namespace bpp
{
/**
 * @brief The EquiprobableSubstitutionModel substitution model for any kind of
 * alphabet. Jukes-Cantor models are specific case of this model,
 * applied to nucleotides or proteins.
 *
 * If the number of states is \f$N\f$.
 *
 * All rates equal:
 * \f[
 * \begin{pmatrix}
 * \ddots & r      & \ldots & r \\
 * r      & \ddots & \ddots & \vdots \\
 * \vdots & \ddots & \ddots & r \\
 * r      & \ldots & r      & \ddots \\
 * \end{pmatrix}
 * \f]
 * \f[
 * \pi = diag\left(\frac{1}{N}, \ldots, \frac{1}{N}\right)
 * \f]
 * Normalization: \f$r\f$ is set so that \f$\sum_i Q_{i,i}\pi_i = -1\f$:
 * \f[
 * S = \begin{pmatrix}
 * -N           & \frac{N}{N-1} & \ldots        & \frac{N}{N-1} \\
 * \frac{N}{N-1} &           -N & \ddots        & \vdots \\
 * \vdots        & \ddots        & \ddots        & \frac{N}{N-1} \\
 * \frac{N}{N-1} & \ldots        & \frac{N}{N-1} & -N \\
 * \end{pmatrix}
 * \f]
 * The normalized generator is obtained by taking the dot product of \f$S\f$ and \f$pi\f$:
 * \f[
 * Q = S . \pi = \begin{pmatrix}
 * -1 & \frac{1}{N-1}     & \ldots & \frac{1}{N-1} \\
 * \frac{1}{N-1} & -1     & \ddots & \vdots \\
 * \vdots       & \ddots & \ddots & \frac{1}{N-1} \\
 * \frac{1}{N-1} & \ldots & \frac{1}{N-1} & -1 \\
 * \end{pmatrix}
 * \f]
 *
 * The eigen values are \f$\left(0, -\frac{N}{N-1}, \ldots, -\frac{N}{N-1}\right)\f$,
 * the left eigen vectors are, by row:
 * \f[
 * U = \begin{pmatrix}
 *   \frac{1}{N} &        \ldots &  \frac{1}{N} &   \frac{1}{N} &  \frac{1}{N} \\
 *  -\frac{1}{N} &        \ldots & -\frac{1}{N} &  \frac{N-1}{N} & -\frac{1}{N} \\
 *         \vdots &        \ddots & \frac{N-1}{N} & -\frac{1}{N}  & -\frac{1}{N} \\
 *  -\frac{1}{N} &        \ddots &        \ddots &         \vdots &        \vdots \\
 *  \frac{N-1}{N} & -\frac{1}{N} &        \ldots &  -\frac{1}{N} & -\frac{1}{N} \\
 * \end{pmatrix}
 * \f]
 * and the right eigen vectors are, by column:
 * \f[
 * U^-1 = \begin{pmatrix}
 *      1 &      0 &  \ldots &      0 &  1 \\
 * \vdots & \vdots &  \ddots & \ddots &  0 \\
 *      1 &      0 &       1 & \ddots & \vdots \\
 *      1 &      1 &       0 & \ldots &  0 \\
 *      1 &     -1 &      -1 & \ldots & -1 \\
 * \end{pmatrix}
 * \f]
 *
 * In addition, a rate_ factor defines the mean rate of the model.
 *
 * The probabilities of changes are computed analytically using the formulas:
 * \f[
 * P_{i,j}(t) = \begin{cases}
 * \frac{1}{N} + \frac{N-1}{N}e^{- rate\_ * \frac{N}{N-1}t}& \text{if $i=j$}, \\
 * \frac{1}{N} - \frac{1}{N}e^{- rate\_ * \frac{N}{N-1}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 *
 * First and second order derivatives are also computed analytically using the formulas:
 * \f[
 * \frac{\partial P_{i,j}(t)}{\partial t} = rate\_ * \begin{cases}
 * -e^{- rate\_ * \frac{N}{N-1}t}           & \text{if $i=j$}, \\
 * \frac{1}{N-1}e^{- rate\_ * \frac{N}{N-1}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 * \f[
 * \frac{\partial^2 P_{i,j}(t)}{\partial t^2} =  rate\_^2 * \begin{cases}
 * \frac{N}{N-1}e^{- rate\_ * \frac{N}{N-1}t}  & \text{if $i=j$}, \\
 * -\frac{N}{(N-1)^2}e^{- rate\_ * \frac{N}{N-1}t} & \text{otherwise}.
 * \end{cases}
 * \f]
 *
 */
class EquiprobableSubstitutionModel :
  public AbstractReversibleSubstitutionModel
{
private:
  mutable double exp_;
  mutable RowMatrix<double> p_;
  std::shared_ptr<FrequencySet> freqSet_;

public:
  /**
   * @brief Build a simple equiprobable model, with original equilibrium frequencies.
   *
   * @param alpha An alphabet.
   */
  EquiprobableSubstitutionModel(const Alphabet* alpha);

  /**
   * @brief Build an equiprobable model with special equilibrium frequencies.
   *
   * @param alpha An alphabet.
   * @param freqSet A pointer toward a frequencies set, which will be owned by this instance.
   * @param initFreqs Tell if the frequency set should be initialized with the original JTT92 values.
   * Otherwise, the values of the set will be used.
   */
  EquiprobableSubstitutionModel(const Alphabet* alpha, FrequencySet* freqSet, bool initFreqs = false);

  EquiprobableSubstitutionModel(const EquiprobableSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    AbstractReversibleSubstitutionModel(model),
    exp_(model.exp_),
    p_(model.p_),
    freqSet_(std::shared_ptr<FrequencySet>(model.freqSet_->clone()))
  {}

  EquiprobableSubstitutionModel& operator=(const EquiprobableSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    AbstractReversibleSubstitutionModel::operator=(model);
    exp_ = model.exp_;
    p_   = model.p_;
    freqSet_ = std::shared_ptr<FrequencySet>(model.freqSet_->clone());
    return *this;
  }

  virtual ~EquiprobableSubstitutionModel() {}

  EquiprobableSubstitutionModel* clone() const { return new EquiprobableSubstitutionModel(*this); }

public:
  double Pij_t    (size_t i, size_t j, double d) const;
  double dPij_dt  (size_t i, size_t j, double d) const;
  double d2Pij_dt2(size_t i, size_t j, double d) const;
  const Matrix<double>& getPij_t    (double d) const;
  const Matrix<double>& getdPij_dt  (double d) const;
  const Matrix<double>& getd2Pij_dt2(double d) const;

  std::string getName() const
  {
    if (freqSet_->getNamespace().find("+F.") != std::string::npos)
      return "Equi+F";
    else
      return "Equi";
  }

  void fireParameterChanged(const ParameterList& parameters)
  {
    freqSet_->matchParametersValues(parameters);
    freq_ = freqSet_->getFrequencies();
    AbstractReversibleSubstitutionModel::fireParameterChanged(parameters);
  }

  void setFrequencySet(const FrequencySet& freqSet)
  {
    freqSet_ = std::shared_ptr<FrequencySet>(freqSet.clone());
    resetParameters_();
    addParameters_(freqSet_->getParameters());
  }

  const std::shared_ptr<FrequencySet> getFrequencySet() const { return freqSet_; }

  void setFreq(std::map<int, double>& freq);

protected:
  /**
   * In the case of the model of Jukes & Cantor, this method is useless since
   * the generator is fixed! No matrice can be changed... This method is only
   * used in the constructor of the class.
   */
  void updateMatrices();
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_EQUIPROBABLESUBSTITUTIONMODEL_H
