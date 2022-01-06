//
// File: CompoundTransitionModel.h
// Authors:
//   Anaïs Prud'homme
// Date :
//   Vendredi 3 décembre 2021
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
 * @brief Transition models defined as a composition of nested
 * substitution models.
 * @author Anaïs Prud'homme
 *
 */
class CompoundTransitionModel :
  public AbstractMixedTransitionModel
{
private:
  std::map<std::string, DiscreteDistribution*> distributionMap_;

protected:
  int from_, to_;

public:
  CompoundTransitionModel(
    const Alphabet* alpha,
    TransitionModel* model,
    std::map<std::string, DiscreteDistribution*> parametersDistributionsList,
    int ffrom = -1, int tto = -1);

  CompoundTransitionModel(const CompoundTransitionModel&);

  CompoundTransitionModel& operator=(const CompoundTransitionModel&);

  virtual ~CompoundTransitionModel();

  CompoundTransitionModel* clone() const { return new CompoundTransitionModel(*this); }

public:
  std::string getName() const { return "CompoundModel"; }

  void updateMatrices();

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

  /*
   *@brief Returns the vector of numbers of the submodels in the
   * mixture that match a description of the parameters numbers.
   *
   **@param desc is the description of the class indexes of the mixed
   **parameters. Syntax is like: kappa_1,gamma_3,delta_2
   *
   */

  Vuint getSubmodelNumbers(const std::string& desc) const;

  /**
   * @brief sets the eq frequencies of the first nested model, and
   * adapts the parameters at best to it (surely there is a better way
   * to manage this).
   *
   */

  void setFreq(std::map<int, double>&);

  /**
   * @brief returns the DiscreteDistribution associated with a given
   * parameter name.
   * @param parName name of the parameter
   **/

  const DiscreteDistribution* getDistribution(std::string& parName) const;

  /**
   *@brief Numbers of the states between which the substitution rates
   * of all the submodels must be equal. If they are set to -1, this
   * constraint does not exist among the submodels.
   *
   */
  int from() const { return from_; }
  int to() const { return to_; }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CompoundTransitionModel_H
