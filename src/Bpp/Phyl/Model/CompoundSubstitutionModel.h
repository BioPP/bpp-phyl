//
// File: CompoundSubstitutionModel.h
// Authors:
//   Anaïs Prud'homme
//   Date: vendredi 3 décembre à 11h30
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

#ifndef BPP_PHYL_MODEL_CompoundSubstitutionModel_H
#define BPP_PHYL_MODEL_CompoundSubstitutionModel_H

#include <Bpp/Numeric/VectorTools.h>
#include <cstring> // C lib for string copy
#include <map>
#include <string>
#include <vector>

#include "CompoundTransitionModel.h"

namespace bpp
{
class CompoundSubstitutionModel :
  public CompoundTransitionModel
{
public:
  /**
   * @brief Constructor of a CompoundSubstitutionModel, where all
   * the models have rate 1 and equal probability. => changer
   *
   * @param alpha pointer to the Alphabet
   * @param vpModel vector of pointers to SubstitutionModels. All the
   *   SubstitutionModels are owned by the instance.
   * @warning providing a vpModel with size 0 will generate a segmentation fault!
   */
  CompoundSubstitutionModel(const Alphabet* alpha, std::vector<std::shared_ptr<TransitionModel>> vpModel) :
    AbstractParameterAliasable("Mixture."),
    AbstractTransitionModel(alpha, vpModel.size() ? vpModel[0]->shareStateMap() : 0, "Mixture."),
    CompoundTransitionModel(alpha, vpModel)
  {
    /*
     * Check that all models are substitutionmodels
     */

    for (const auto& model:vpModel)
    {
      if (!dynamic_cast<const SubstitutionModel*>(model.get()))
        throw Exception("CompoundSubstitutionModel can only be built with SubstitutionModels, not " + model->getName());
    }
  }

  /**
   * @brief Constructor of a CompoundSubstitutionModel.
   *
   * @param alpha pointer to the Alphabet
   * @param vpModel vector of pointers to SubstitutionModels. All the
   *   SubstitutionModels are owned by the instance.
   * @param vproba vector of the probabilities of the models
   * @param vrate vector of the rates of the models
   * @warning providing a vpModel with size 0 will generate a segmentation fault!
   *
   * See above the constraints on the rates and the probabilities of
   * the vectors.
   */

  CompoundSubstitutionModel(
    const Alphabet* alpha,
    std::vector<std::shared_ptr<TransitionModel>> vpModel,
    Vdouble& vproba, Vdouble& vrate) :
    AbstractParameterAliasable("Mixture."),
    AbstractTransitionModel(alpha, vpModel.size() ? vpModel[0]->shareStateMap() : 0, "Mixture."),
    CompoundTransitionModel(alpha, vpModel, vproba, vrate)
  {
    /*
     * Check that all models are substitutionmodels
     */

    for (const auto& model:vpModel)
    {
      if (!dynamic_cast<const SubstitutionModel*>(model.get()))
        throw Exception("CompoundSubstitutionModel can only be built with SubstitutionModels, not " + model->getName());
    }
  }


  CompoundSubstitutionModel(const CompoundSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    AbstractTransitionModel(model),
    CompoundTransitionModel(model)
  {}


  CompoundSubstitutionModel& operator=(const CompoundSubstitutionModel& model)
  {
    CompoundTransitionModel::operator=(model);
    return *this;
  }

  CompoundSubstitutionModel* clone() const { return new CompoundSubstitutionModel(*this); }

public:
  /**
   * @brief retrieve a pointer to the subsitution model with the given name.
   *
   * Return Null if not found.
   *
   */
  const SubstitutionModel* getSubModel(const std::string& name) const
  {
    return dynamic_cast<const SubstitutionModel*>(getModel(name));
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CompoundSubstitutionModel_H
