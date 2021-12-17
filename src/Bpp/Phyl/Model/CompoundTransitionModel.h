//
// File: CompoundTransitionModel.h
// Authors:
//   Anaïs Prud'homme
//   Date: vendredi 3 décembre à 11h
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

#ifndef BPP_PHYL_MODEL_COMPOUNDTRANSITIONMODEL_H
#define BPP_PHYL_MODEL_COMPOUNDTRANSITIONMODEL_H

#include <Bpp/Numeric/VectorTools.h>
#include <cstring> // C lib for string copy
#include <map>
#include <string>
#include <vector>

#include "AbstractTransitionModel.h"

namespace bpp
{
/**
 * @brief 
 * @author Anaïs Prud'homme
 */

class CompoundTransitionModel :
  public AbstractTransitionModel
{
public:
  /**
   * @brief Constructor of a CompoundTransitionModel, where all
   * the models have rate 1 and equal probability.
   * @param alpha pointer to the Alphabet
   * @param vpModel vector of pointers to TransitionModels. All the
   *   TransitionModels are owned by the instance.
   * @warning providing a vpModel with size 0 will generate a segmentation fault!
   */
  CompoundTransitionModel(
    const Alphabet* alpha,
    std::vector<std::shared_ptr<TransitionModel>> vpModel);

  /**
   * @brief Constructor of a CompoundTransitionModel.
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

  CompoundTransitionModel(
    const Alphabet* alpha,
    std::vector<std::shared_ptr<TransitionModel>> vpModel,
    Vdouble& vproba, Vdouble& vrate);

  CompoundTransitionModel(const CompoundTransitionModel&);

  CompoundTransitionModel& operator=(const CompoundTransitionModel&);

  virtual ~CompoundTransitionModel();

  CompoundTransitionModel* clone() const { return new CompoundTransitionModel(*this); }

public:
  std::string getName() const { return "Mixture"; }

  /**
   * @brief retrieve a pointer to the submodel with the given name.
   *
   * Return Null if not found.
   *
   */

  const TransitionModel* getModel(const std::string& name) const;

  const TransitionModel* getModel(size_t i) const
  {
    return AbstractMixedTransitionModel::getNModel(i);
  }

  // TransitionModel* getModel(size_t i)
  // {
  //   return AbstractMixedTransitionModel::getModel(i);
  // }

  void updateMatrices();

  /**
   * @brief Sets the rates of the submodels to follow the constraint
   * that the mean rate of the mixture equals rate_.
   * @param vd a vector of positive values such that the rates of
   * the respective submodels are in the same proportions (ie this
   * vector does not need to be normalized).
   */
  virtual void setVRates(const Vdouble& vd);

  /**
   * @brief Returns the vector of numbers of the submodels in the
   * mixture that match a description of the parameters numbers.
   *
   * @param desc is the description of the class indexes of the mixed
   * parameters. Syntax is like: kappa_1,gamma_3,delta_2
   */
  Vuint getSubmodelNumbers(const std::string& desc) const;

  /**
   * @brief applies setFreq to all the models of the mixture and
   * recovers the parameters values.
   */
  void setFreq(std::map<int, double>&);
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CompoundTransitionModel_H
