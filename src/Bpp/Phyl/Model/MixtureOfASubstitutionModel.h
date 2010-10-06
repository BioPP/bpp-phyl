//
// File: MixtureOfASubstitutionModel.h
// Created by: David Fournier, Laurent Gueguen
//

/*
   Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _MIXTUREOFASUBSTITUTIONMODEL_H_
#define _MIXTUREOFASUBSTITUTIONMODEL_H_

#include "MixedSubstitutionModel.h"
#include <Bpp/Numeric/Prob.all>
#include <Bpp/Numeric/VectorTools.h>

#include <vector>
#include <string>
#include <map>
#include <cstring> // C lib for string copy

namespace bpp
{
/**
 * @brief Substitution models defined as a mixture of "simple"
 * substitution models.
 * @author Laurent Guéguen
 *
 * All the models are of the same type (for example T92 or GY94),
 * and their parameter values can follow discrete distributions.
 *
 * In this kind of model, there is no generator or transition
 * probabilities.
 *
 * There is a map with connection from parameter names to discrete
 * distributions, and then a related vector of "simple" substitution
 * models for all the combinations of parameter values.
 *
 * For example:
 * HKY85(kappa=Gamma(n=3,alpha=2,beta=5),
 *       theta=TruncExponential(n=4,lambda=0.2,tp=1),
 *       theta1=0.4,
 *       theta2=TruncExponential(n=5,lambda=0.6,tp=1))
 *
 * defines 3*4*5=60 different HKY85 models.
 *
 * If a distribution parameter does not respect the constraints of
 * this parameter, there is an Exception at the creation of the
 * wrong mdel, if any.
 *
 * When used through a MixedTreeLikelihood objets, all the models have
 * a specific probability, defined through the probabilities of the
 * several parameter distributions. The computing of the likelihoods
 * and probabilities are the expectation of the "simple" models
 * values.
 *
 */
class MixtureOfASubstitutionModel :
  public MixedSubstitutionModel
{
private:
  std::map<std::string, DiscreteDistribution*> distributionMap_;

  std::vector<SubstitutionModel*> modelsContainer_;
  std::vector<double> probas_;
  
public:
  MixtureOfASubstitutionModel(const Alphabet* alpha,
                         SubstitutionModel* model,
                         std::map<std::string, DiscreteDistribution*> parametersDistributionsList) throw(Exception);

  MixtureOfASubstitutionModel(const MixtureOfASubstitutionModel&);
  
  MixtureOfASubstitutionModel& operator=(const MixtureOfASubstitutionModel&);

  ~MixtureOfASubstitutionModel();

  MixtureOfASubstitutionModel* clone() const { return new MixtureOfASubstitutionModel(*this); }

public:
  /**
   * @brief Returns a specific model from the mixture
   */
  const SubstitutionModel* getNModel(unsigned int i) const
  {
    return modelsContainer_[i];
  }

  SubstitutionModel* getNModel(unsigned int i)
  {
    return modelsContainer_[i];
  }

  /**
   * @brief Returns the  probability of a specific model from the mixture
   */
  
  double getNProbability(unsigned int i) const
  {
    return probas_[i];
  }
  
  const std::vector<double>& getProbabilities() const
  {
    return probas_;
  }
  
  unsigned int getNumberOfModels() const
  {
    return modelsContainer_.size();
  }

  std::string getName() const { return "MixtureOfASubstitutionModel"; }

  void updateMatrices();
  
  unsigned int getNumberOfStates() const;

  double Pij_t(unsigned int i, unsigned int j, double t) const;
  double dPij_dt(unsigned int i, unsigned int j, double t) const;
  double d2Pij_dt2(unsigned int i, unsigned int j, double t) const;
  const Matrix<double>& getPij_t(double t) const;
  const Matrix<double>& getdPij_dt(double t) const;
  const Matrix<double>& getd2Pij_dt2(double t) const;
  const Vdouble& getFrequencies();
  double freq(unsigned int i) const;

};
} // end of namespace bpp.

#endif  // _MIXTUREOFASUBSTITUTIONMODEL_H_
