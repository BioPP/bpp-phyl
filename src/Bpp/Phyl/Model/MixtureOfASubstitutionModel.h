//
// File: MixtureOfASubstitutionModel.h
// Authors:
//   Laurent Gueguen
//   Date: lundi 13 septembre 2010, Ã 21h 31
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

#ifndef BPP_PHYL_MODEL_MIXTUREOFASUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_MIXTUREOFASUBSTITUTIONMODEL_H

#include <Bpp/Numeric/VectorTools.h>
#include <cstring> // C lib for string copy
#include <map>
#include <string>
#include <vector>

#include "MixtureOfATransitionModel.h"

namespace bpp
{
class MixtureOfASubstitutionModel :
  public MixtureOfATransitionModel
{
public:
  /**
   * @brief Constructor of a MixtureOfASubstitutionModel, where all
   * the models have rate 1 and equal probability.
   *
   * @param alpha pointer to the Alphabet
   * @param model pointer to the SubstitutionModel that will be mixed
   * @param parametersDistributionsList list from parameters names towards discrete distributions to will define the mixtures.
   * @param ffrom   index of the starting codon that will be used to homogeneize the rates of the submodels
   * @param tto     index of the arriving codon that will be used to homogeneize the rates of the submodels
   *
   *   If ffrom and tto are not -1, for all submodels the transition
   *   rate ffrom->tto is the same. Otherwise, all submodels are
   *   normalized to have a substitution/time unit at equilibrium.
   */
  MixtureOfASubstitutionModel(
      std::shared_ptr<const Alphabet> alpha,
      std::unique_ptr<SubstitutionModelInterface> model,
      std::map<std::string, std::unique_ptr<DiscreteDistributionInterface>>& parametersDistributionsList,
      int ffrom = -1,
      int tto = -1) :
    AbstractParameterAliasable(model->getNamespace()),
    AbstractTransitionModel(alpha, model->getStateMap(), model->getNamespace()),
    MixtureOfATransitionModel(alpha, move(model), parametersDistributionsList, ffrom, tto)
  {}

  MixtureOfASubstitutionModel(const MixtureOfASubstitutionModel& model) :
    AbstractParameterAliasable(model),
    AbstractTransitionModel(model),
    MixtureOfATransitionModel(model)
  {}


  MixtureOfASubstitutionModel& operator=(const MixtureOfASubstitutionModel& model)
  {
    MixtureOfATransitionModel::operator=(model);
    return *this;
  }

  MixtureOfASubstitutionModel* clone() const override { return new MixtureOfASubstitutionModel(*this); }

protected:

  void updateMatrices_() override
  {
    MixtureOfATransitionModel::updateMatrices_();
    // setting the rates, if to_ & from_ are different from -1

    if (to_ >= 0 && from_ >= 0)
    {
      Vdouble vd;
      for (size_t j = 0; j < modelsContainer_.size(); ++j)
      {
        vd.push_back(1 / subNModel(j).Qij(static_cast<size_t>(from_), static_cast<size_t>(to_)));
      }
      setVRates(vd);
    }
  }

public:
  
  /**
   * @brief retrieve a pointer to the subsitution model with the given name.
   */
  const SubstitutionModelInterface& subModel(const std::string& name) const
  {
    return dynamic_cast<const SubstitutionModelInterface&>(model(name));
  }

  const SubstitutionModelInterface& subNModel(size_t i) const
  {
    return dynamic_cast<const SubstitutionModelInterface&>(nModel(i));
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_MIXTUREOFASUBSTITUTIONMODEL_H
