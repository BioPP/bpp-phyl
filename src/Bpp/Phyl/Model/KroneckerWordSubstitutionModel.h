//
// File: KroneckerWordSubstitutionModel.h
// Authors:
//   Laurent Gueguen
// Created: mercredi 27 juillet 2016, ÃÂ  00h 12
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

#ifndef BPP_PHYL_MODEL_KRONECKERWORDSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_KRONECKERWORDSUBSTITUTIONMODEL_H


#include "AbstractKroneckerWordSubstitutionModel.h"

// From bpp-core
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/BppVector.h>

namespace bpp
{
/**
 * @brief Basal class for words of  substitution models with multiple
 * substitutions.
 * @author Laurent Guéguen
 *
 * The equilibrium frequency of each word is due to the combination of
 * the multiple substitutions.
 */
class KroneckerWordSubstitutionModel :
  public AbstractKroneckerWordSubstitutionModel
{
public:
  /**
   * @brief Build a new KroneckerWordSubstitutionModel object from a
   * Vector of pointers to SubstitutionModels.
   *
   * @param modelList the list of substitution models to use, in
   *   the order of the positions in the words from left to right. All
   *   the models must be different objects to avoid parameters
   *   redundancy, otherwise only the first model is used. The used models
   *   are owned by the instance.
   * @param prefix the Namespace.
   */
  KroneckerWordSubstitutionModel(ModelList& modelList, const std::string& prefix = "");

  /**
   * @brief Build a new KroneckerWordSubstitutionModel object from a
   * Vector of pointers to SubstitutionModels.
   *
   * @param modelList the list of substitution models to use, in
   *   the order of the positions in the words from left to right. All
   *   the models must be different objects to avoid parameters
   *   redundancy, otherwise only the first model is used. The used models
   *   are owned by the instance.
   * @param vPos a vector of sets of simultaneously changing
   *   positions.
   * @param prefix the Namespace.
   */
  KroneckerWordSubstitutionModel(ModelList& modelList,
                                 const std::vector<std::set< size_t> >& vPos,
                                 const std::string& prefix = "");

  /**
   * @brief Build a new KroneckerWordSubstitutionModel object from a
   * pointer to an SubstitutionModel and a number of
   * desired models.
   *
   * @param pmodel pointer to the substitution model to use in all the
   *  positions. It is owned by the instance.
   * @param num The number of models involved.
   * @param prefix the Namespace.
   */
  KroneckerWordSubstitutionModel(
      std::unique_ptr<SubstitutionModelInterface> pmodel,
      unsigned int num,
      const std::string& prefix = "");

  /**
   * @brief Build a new KroneckerWordSubstitutionModel object from a
   * pointer to an SubstitutionModel and a number of
   * desired models.
   *
   * @param pmodel pointer to the substitution model to use in all the
   *  positions. It is owned by the instance.
   * @param num The number of models involved.
   * @param vPos a vector of sets of simultaneously changing
   *   positions.
   * @param prefix the Namespace.
   */
  KroneckerWordSubstitutionModel(
      std::unique_ptr<SubstitutionModelInterface> pmodel,
      unsigned int num,
      const std::vector<std::set< size_t> >& vPos,
      const std::string& prefix = "");

  virtual ~KroneckerWordSubstitutionModel() {}

  KroneckerWordSubstitutionModel* clone() const override {
    return new KroneckerWordSubstitutionModel(*this);
  }

protected:

  /**
   * @brief Constructor for the derived classes only
   */
  KroneckerWordSubstitutionModel(
      std::shared_ptr<const Alphabet> alph,
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& prefix = "");

  void completeMatrices_() override {}

public:
  virtual std::string getName() const override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_KRONECKERWORDSUBSTITUTIONMODEL_H
