//
// File: AbstractKroneckerWordSubstitutionModel.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: lundi 25 juillet 2016, ÃÂ  17h 00
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

#ifndef BPP_PHYL_MODEL_ABSTRACTKRONECKERWORDSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_ABSTRACTKRONECKERWORDSUBSTITUTIONMODEL_H

#include <set>

#include "AbstractWordSubstitutionModel.h"

// From bpp-core:
#include <Bpp/Numeric/Prob/Simplex.h>

namespace bpp
{

/**
 * @brief Abstract Kronecker Word Model.
 *
 * Objects of this class are built from several substitution models.
 * Each model corresponds to a position in the word. No model is
 * directly accessible.
 *
 * A substitution rate between words is the product of the rates of
 * changing letters of position generators.
 *
 * Either all substitutions are allowed, or any subset of them is
 * allowed. For example, in case of 3 words model, if @f$Q@f$ is the
 * generator of the position model, either all substitutions are
 * allowed (ie generator is @f$ Q_1 \otimes Q_1 \otimes Q_1 @f$ with
 * @f$Q_1@f$ is @f$Q@f$ with diagonal of 1), or all and only three
 * letters substitutions are allowed (ie generator is @f$ Q_0
 * \otimes Q_0 \otimes Q_0 @f$ with @f$Q_0@f$ is @f$Q@f$ with
 * diagonal of 0), or all and only two letters substitutions are
 * allowed (generator is @f$ I \otimes Q_0 \otimes Q_0 + Q_0 \otimes
 * I \otimes Q_0 + Q_0 \otimes Q_0 \otimes I @f$), or all and only
 * two neighbors letters substitutions are allowed (generator is @f$
 * I \otimes Q_0 \otimes Q_0 + Q_0 \otimes Q_0 \otimes I @f$).
 *
 * Default construction is all substitutions, but a vector of
 * vectors of simultaneously changing positions can be given.
 *
 */
class AbstractKroneckerWordSubstitutionModel :
  public AbstractWordSubstitutionModel
{
private:
  /**
   * @brief vector of sets of simultaneously changing positions.
   *
   * Shortcut : if empty, any number of any positions can change
   * simultaneously
   */
  std::vector<std::set< size_t> > sChangingPos_;

  /**
   * @brief vector of generators for computation purposes
   */
  std::vector< RowMatrix<double> > vGenerators_;

protected:
  void initGenerators_();

  /**
   * @brief checks that the vector of changing positions is valid
   */
  bool checkChangingPositions_();

private:
  
  /**
   * @brief First fill of the generator, from the position model
   */
  void fillBasicGenerator_();

public:

  /**
   * @brief Build a new AbstractKroneckerWordSubstitutionModel
   * object from a vector of pointers to SubstitutionModels,
   * allowing for all substitutions.
   *
   * @param modelList the list of substitution models to use, in
   *   the order of the positions in the words from left to right. All
   *   the models must be different objects to avoid parameters
   *   redundancy, otherwise only the first model is used. The used models
   *   are owned by the instance.
   * @param prefix the Namespace.
   */
  AbstractKroneckerWordSubstitutionModel(
    ModelList& modelList,
    const std::string& prefix);

  /**
   * @brief Build a new AbstractKroneckerWordSubstitutionModel
   * object from a vector of pointers to SubstitutionModels,
   * allowing for a given vector of simultaneous substitutions.
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
  AbstractKroneckerWordSubstitutionModel(
    ModelList& modelList,
    const std::vector<std::set< size_t> >& vPos,
    const std::string& prefix);

  /**
   * @brief Build a new AbstractWordSubstitutionModel object from a
   * pointer to an SubstitutionModel and a number of desired models,
   * allowing for all substitutions.
   *
   * @param pmodel A pointer to the substitution model to use in all
   * the positions. It will be owned by the instance.
   * @param num The number of models involved.
   * @param prefix the Namespace.
   */
  AbstractKroneckerWordSubstitutionModel(
    std::unique_ptr<SubstitutionModelInterface>& pmodel,
    unsigned int num,
    const std::string& prefix);

  /**
   * @brief Build a new AbstractWordSubstitutionModel object from a
   * pointer to an SubstitutionModel and a number of desired models,
   * allowing for a given vector of simultaneous substitutions.
   *
   * @param pmodel A pointer to the substitution model to use in all
   * the positions. It will be owned by the instance.
   * @param num The number of models involved.
   * @param vPos a vector of vectors of simultaneously changing
   *   positions.
   * @param prefix the Namespace.
   */
  AbstractKroneckerWordSubstitutionModel(
    std::unique_ptr<SubstitutionModelInterface>& pmodel,
    unsigned int num,
    const std::vector<std::set< size_t> >& vPos,
    const std::string& prefix);


  virtual ~AbstractKroneckerWordSubstitutionModel() {}

  AbstractKroneckerWordSubstitutionModel(const AbstractKroneckerWordSubstitutionModel&);

  AbstractKroneckerWordSubstitutionModel& operator=(const AbstractKroneckerWordSubstitutionModel&);

protected:
  
  /**
   * @brief Constructor for the derived classes only
   */
  AbstractKroneckerWordSubstitutionModel(
      std::shared_ptr<const Alphabet> alph,
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& prefix);

  void setChangingPositions(const std::vector<std::set<size_t> >& vPos);
};

} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ABSTRACTKRONECKERWORDSUBSTITUTIONMODEL_H
