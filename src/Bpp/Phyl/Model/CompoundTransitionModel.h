//
// File: CompoundTransitionModel.h
// Authors:
//   Anaïs Prud'homme
// Date :
//   Vendredi 3 décembre 2021 à 11h30
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

#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <cstring> // C lib for string copy
#include <map>
#include <string>
#include <vector>

#include "AbstractSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Transition models defined as a composition of nested
 * substitution models.
 * @author Anaïs Prud'homme
 *
 */
class CompoundTransitionModel :
  public virtual AbstractTransitionModel
{

protected:
  int from_, to_;

   /**
   * @brief vector of pointers to TransitionModels.
   *
   * Beware: these TransitionModels are owned by the object, so
   * will be deleted at destruction
   */
  std::vector<std::shared_ptr<TransitionModel> > modelsContainer_;

  /**
   * @brief vector of the probabilities of the models
   */
  std::vector<double> vProbas_;

public:
  CompoundTransitionModel(
    const Alphabet* alpha,
    TransitionModel* model,
    int ffrom = -1, int tto = -1);

  CompoundTransitionModel(const CompoundTransitionModel&);

  CompoundTransitionModel& operator=(const CompoundTransitionModel&);

  virtual ~CompoundTransitionModel(){};

  CompoundTransitionModel* clone() const { return new CompoundTransitionModel(*this); }

public:
  std::string getName() const { return "CompoundModel"; }

  void updateMatrices();

public:
  /**
   * @brief returns the number of models in the composition
   */
  virtual size_t getNumberOfModels() const
  {
    return modelsContainer_.size();
  }

  /**
   * @brief retrieve a pointer to the submodel with the given name.
   *
   * Return Null if not found.
   *
   */

  const TransitionModel* getModel(const std::string& name) const;

  /**
   * @brief Returns a specific model from the composition
   */
  const TransitionModel* getModel(size_t i) const
  {
    return modelsContainer_[i].get();
  }

  TransitionModel* getModel(size_t i)
  {
    return modelsContainer_[i].get();
  }

  Vuint getSubmodelNumbers(const std::string& desc) const;

  /**
   * @brief Returns the probability of a specific model from the
   * mixture
   */
  virtual double getNProbability(size_t i) const
  {
    return vProbas_[i];
  }

  /**
   * @brief Returns the vector of probabilities
   *
   */
  virtual const std::vector<double>& getProbabilities() const
  {
    return vProbas_;
  }

  /**
   * @brief Sets the  probability of a specific model from the mixture
   */
  virtual void setNProbability(size_t i, double prob)
  {
    if (prob < 0)
      prob = 0;
    if (prob > 1)
      prob = 1;

    vProbas_[i] = prob;
  }

  /**
   * @brief From TransitionModel interface
   *
   */

  virtual size_t getNumberOfStates() const;

  virtual const Matrix<double>& getPij_t(double t) const;
  virtual const Matrix<double>& getdPij_dt(double t) const;
  virtual const Matrix<double>& getd2Pij_dt2(double t) const;

  /**
   * @return Says if equilibrium frequencies should be computed (all
   * models are likewise, may be refined)
   */
  bool computeFrequencies() const
  {
    return modelsContainer_[0]->computeFrequencies();
  }

  /**
   * @return Set if equilibrium frequencies should be computed
   */
  void computeFrequencies(bool yn)
  {
    for (auto& sm : modelsContainer_)
    {
      sm->computeFrequencies(yn);
    }
  }

  /**
   * @brief sets the eq frequencies of the first nested model, and
   * adapts the parameters at best to it (surely there is a better way
   * to manage this).
   *
   */

  void setFreq(std::map<int, double>&);


  void setFreqFromData(const SequencedValuesContainer& data, double pseudoCount)
  {
    std::map<int, double> freqs;
    SequenceContainerTools::getFrequencies(data, freqs, pseudoCount);
    setFreq(freqs);
  }

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
