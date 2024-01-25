//
// File: MixtureOfATransitionModel.h
// Authors:
//   David Fournier, Laurent Gueguen
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

#ifndef BPP_PHYL_MODEL_MIXTUREOFATRANSITIONMODEL_H
#define BPP_PHYL_MODEL_MIXTUREOFATRANSITIONMODEL_H

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/VectorTools.h>
#include <cstring> // C lib for string copy
#include <map>
#include <string>
#include <vector>

#include "AbstractMixedTransitionModel.h"

namespace bpp
{
/**
 * @brief Transition models defined as a mixture of nested
 * substitution models.
 * @author Laurent Guéguen
 *
 * All the nested models are of the same type (for example T92 or
 * GY94), and their parameter values can follow discrete
 * distributions.
 *
 * In this kind of model, there is no generator.
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
 * defines 3*4*5=60 different HKY85 nested models with rate one.
 *
 * Optionnal arguments are used to homogeneize the rates of the nested
 * models. Default values sets all the rates to 1, and given values
 * are the two letters (from_ & to_) between which the substitution
 * rates are the same in all nested models.
 *
 * For example:
 * HKY85(kappa=Gamma(n=3,alpha=2,beta=5),
 *       theta=TruncExponential(n=4,lambda=0.2,tp=1),
 *       theta1=0.4,
 *       theta2=TruncExponential(n=5,lambda=0.6,tp=1),
 *       from=A, to=C)
 *
 * defines 3*4*5=60 different HKY85 nested models with the same A->C
 * substitution rate.
 *
 * If a distribution parameter does not respect the constraints of
 * this parameter, there is an Exception at the creation of the
 * wrong model, if any.
 *
 * When used through a MixedTreeLikelihood objets, all the models have
 * a specific probability, defined through the probabilities of the
 * several parameter distributions. The computing of the likelihoods
 * and probabilities are the expectation of the "simple" models
 * values.
 *
 */
class MixtureOfATransitionModel :
  public AbstractMixedTransitionModel
{
private:
  std::map<std::string, std::unique_ptr<DiscreteDistributionInterface>> distributionMap_;

protected:
  int from_, to_;

public:
  MixtureOfATransitionModel(
    std::shared_ptr<const Alphabet> alpha,
    std::unique_ptr<TransitionModelInterface> model,
    std::map<std::string, std::unique_ptr<DiscreteDistributionInterface>>& parametersDistributionsList,
    int ffrom = -1, int tto = -1);

  MixtureOfATransitionModel(const MixtureOfATransitionModel&);

  MixtureOfATransitionModel& operator=(const MixtureOfATransitionModel&);

  virtual ~MixtureOfATransitionModel();

  MixtureOfATransitionModel* clone() const override { return new MixtureOfATransitionModel(*this); }

public:
  std::string getName() const override { return "MixedModel"; }

  /**
   * @brief retrieve a pointer to the submodel with the given name.
   */
  const TransitionModelInterface& model(const std::string& name) const override
  {
    size_t nbmod = getNumberOfModels();

    for (size_t i = 0; i < nbmod; ++i)
    {
      auto& model = nModel(i);
      if (model.getName() == name)
        return model;
    }

    throw NullPointerException("MixtureOfATransitionModel::model(). No model with the specified name.");
  }

  const TransitionModelInterface& model(size_t i) const
  {
    return AbstractMixedTransitionModel::nModel(i);
  }
  
  /**
   * @brief Returns the vector of numbers of the submodels in the
   * mixture that match a description of the parameters numbers.
   *
   * @param desc is the description of the class indexes of the mixed
   * parameters. Syntax is like: kappa_1,gamma_3,delta_2
   */
  Vuint getSubmodelNumbers(const std::string& desc) const override;

  /**
   * @brief sets the eq frequencies of the first nested model, and
   * adapts the parameters at best to it (surely there is a better way
   * to manage this).
   */
  void setFreq(std::map<int, double>&) override;

  /**
   * @brief Tells whether a DiscreteDistribution is associated with a given
   * parameter name.
   * @param parName name of the parameter
   */
  bool hasDistribution(std::string& parName) const
  {
    return (distributionMap_.find(parName) != distributionMap_.end());
  }

  /**
   * @brief returns the DiscreteDistribution associated with a given
   * parameter name.
   * @param parName name of the parameter
   */
  const DiscreteDistributionInterface& distribution(std::string& parName) const
  {
    if (distributionMap_.find(parName) != distributionMap_.end())
      return *distributionMap_.find(parName)->second;
    else
      throw Exception("MixtureOfATransitionModel::distribution(). No distribution with name: '" + parName + "'.");
  }


  /**
   *@brief Numbers of the states between which the substitution rates
   * of all the submodels must be equal. If they are set to -1, this
   * constraint does not exist among the submodels.
   *
   */
  int from() const { return from_; }

  int to() const { return to_; }

protected:

  void updateMatrices_() override;
  
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_MIXTUREOFATRANSITIONMODEL_H
