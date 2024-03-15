// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_MIXTUREOFTRANSITIONMODELS_H
#define BPP_PHYL_MODEL_MIXTUREOFTRANSITIONMODELS_H

#include <Bpp/Numeric/VectorTools.h>
#include <cstring> // C lib for string copy
#include <map>
#include <string>
#include <vector>

#include "AbstractMixedTransitionModel.h"

namespace bpp
{
/**
 * @brief Transition models defined as a mixture of several
 * substitution models.
 * @author Laurent Guéguen
 *
 * All the models can be of different types (for example T92 or
 * GY94), and each model has a specific probability and rate.
 *
 *
 * The probabilities and rates of the models are independent
 * parameters, handled directly, under the constraint that the
 * expectation of the rates on the distribution of the models must
 * equal one.
 *
 * If there are @f$n@f$ models, @f$p_i@f$ is the probability of
 * model i (@f$\sum_{i=1}^{n} p_i = 1@f$) and the probabilities
 * are defined by relative probabilities parameters @f$rp_i@f$
 * (called "relprobai") with:
 * @f[
 * 1 <= i < n, p_i = (1-rp_1)*(1-rp_2)...(1-rp_{i-1})*rp_{i}
 * @f]
 * @f[
 * p_n = (1-rp_1)*(1-rp_2)...(1-rp_{n-1})
 * @f]
 * and
 * @f[
 * \forall 1 <= i < n, rp_i = \frac{p_i}{1-(p_1+...+p_{i-1})}
 * @f]
 * where @f$p_i@f$ stands for the probability of model @f$i@f$.
 *
 *
 * If there are @f$n@f$ models, @f$\rho_i@f$ is the rate and @f$p_i@f$
 * is the probability of model i (@f$\sum_{i=1}^{n} p_i * \rho_i =
 * 1@f$), the rates are defined by relative rates parameters
 * @f$r_i@f$ (called "relratei") with:
 * @f[
 * 1 <= i < n, \rho_i = K.(1-r_1)*(1-r_2)...(1-r_{i-1})*r_{i}
 * @f]
 * @f[
 * \rho_n = K.(1-r_1)*(1-r_2)*...*(1-r_{n-1})
 * @f]
 *
 * with @f[ K = \frac{1}{\sum_{i=1}^n p_i.\rho_i} @f]
 *
 * And on the reverse:
 *
 * @f[
 * \forall 1 <= i < n, r_i = \frac{\rho_i}{1-(\rho_1+...+\rho_{i-1})} < 1.
 * @f]
 *
 * So the relative rates parameters are independent from the probabilities.
 *
 * For example:
 *
 * Mixture(model1=HKY85(kappa=3), model2=T92(theta=0.1),
 *         model2=L95(gamma=2), relrate1=0.2, relrate2=0.9,
 *         relproba1=0.1,relproba2=0.8)
 *
 * define a model as a mixture of 3 different models: HKY85 has
 * probability 0.1 and rate 0.36, T92 has probability 0.72 and rate 1.3,
 * and L95 has probability 0.18 and rate 0.14.
 *
 * The parameters are named \c "Mixture.relrate1", \c
 * "Mixture.relrate2", \c "Mixture.relproba1", \c
 * "Mixture.relproba2"... in addition to the parameters of the
 * submodels that are prefixed by "Mixture.i_", where i is the order
 * of the model.
 */
class MixtureOfTransitionModels :
  public AbstractMixedTransitionModel
{
public:
  /**
   * @brief Constructor of a MixtureOfTransitionModels, where all
   * the models have rate 1 and equal probability.
   *
   * @param alpha pointer to the Alphabet
   * @param vpModel vector of pointers to TransitionModels. All the
   *   TransitionModels are owned by the instance.
   * @warning providing a vpModel with size 0 will generate a segmentation fault!
   */
  MixtureOfTransitionModels(
    std::shared_ptr<const Alphabet> alpha,
    std::vector<std::unique_ptr<TransitionModelInterface>>& vpModel);

  /**
   * @brief Constructor of a MixtureOfTransitionModels.
   *
   * @param alpha pointer to the Alphabet
   * @param vpModel vector of pointers to TransitionModels. All the
   *   TransitionModels are owned by the instance.
   * @param vproba vector of the probabilities of the models
   * @param vrate vector of the rates of the models
   * @warning providing a vpModel with size 0 will generate a segmentation fault!
   *
   * See above the constraints on the rates and the probabilities of
   * the vectors.
   */
  MixtureOfTransitionModels(
    std::shared_ptr<const Alphabet> alpha,
    std::vector<std::unique_ptr<TransitionModelInterface>>& vpModel,
    Vdouble& vproba, Vdouble& vrate);

  MixtureOfTransitionModels(const MixtureOfTransitionModels&);

  MixtureOfTransitionModels& operator=(const MixtureOfTransitionModels&);

  virtual ~MixtureOfTransitionModels();

  MixtureOfTransitionModels* clone() const override {
    return new MixtureOfTransitionModels(*this);
  }

public:
  std::string getName() const override { return "Mixture"; }

  /**
   * @brief retrieve a pointer to the submodel with the given name.
   */
  const TransitionModelInterface& model(const std::string& name) const override;

  const TransitionModelInterface& model(size_t i) const
  {
    return AbstractMixedTransitionModel::nModel(i);
  }

  void updateMatrices_() override;

  /**
   * @brief Sets the rates of the submodels to follow the constraint
   * that the mean rate of the mixture equals rate_.
   * @param vd a vector of positive values such that the rates of
   * the respective submodels are in the same proportions (ie this
   * vector does not need to be normalized).
   */
  virtual void setVRates(const Vdouble& vd) override;

  /**
   * @brief Returns the vector of numbers of the submodels in the
   * mixture that match a description of the parameters numbers.
   *
   * @param desc is the description of the class indexes of the mixed
   * parameters. Syntax is like: kappa_1,gamma_3,delta_2
   */
  Vuint getSubmodelNumbers(const std::string& desc) const override;

  /**
   * @brief applies setFreq to all the models of the mixture and
   * recovers the parameters values.
   */
  void setFreq(std::map<int, double>&) override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_MIXTUREOFTRANSITIONMODELS_H
