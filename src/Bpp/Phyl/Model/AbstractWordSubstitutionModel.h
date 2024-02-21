//
// File: AbstractWordSubstitutionModel.h
// Authors:
//   Laurent Gueguen
// Created: 2009-01-08 00:00:00
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

#ifndef BPP_PHYL_MODEL_ABSTRACTWORDSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_ABSTRACTWORDSUBSTITUTIONMODEL_H


#include "AbstractSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/WordAlphabet.h>

// From the STL:
#include <vector>

namespace bpp
{
/**
 * @brief A list of models, for building a WordSubstitutionModel
 *
 * Instance of ModelList contains a vector of pointers toward a substitution model which they do not own. Several methods are provided to build and check Alphabet and StateMaps.
 */
class ModelList
{
protected:
  /**
   * @brief Position-specific models are stored as shared_ptr to allow several positions to share the same model. The constructors, however, take unique_ptr pointers to ensure that the pointers are not shared outside the model instance.
   */
  std::vector<std::shared_ptr<SubstitutionModelInterface>> models_;
  std::shared_ptr<WordAlphabet> wordAlphabet_;

public:
  /**
   * @brief Create a ModelList from one template substitution model.
   *
   * @param models A vector of pointers toward substitution model objects.
   *
   * !! All pointers of the vector will be emptied.
   */
  
  ModelList(std::vector<std::unique_ptr<SubstitutionModelInterface>>& models) :
    models_(models.size()), wordAlphabet_(nullptr)
  {
    std::vector<std::shared_ptr<const Alphabet>> alphabets(models.size());
    for (size_t i = 0; i < models.size(); ++i)
    {
      alphabets[i] = models[i]->getAlphabet();
      models_[i] = std::move(models[i]);
    }
    wordAlphabet_ = std::make_shared<WordAlphabet>(alphabets);
  }

 private:
  ModelList(const ModelList& ml) {}

  ModelList& operator=(const ModelList& ml) { return *this; }

public:
  size_t size() const { return models_.size(); }

  std::shared_ptr<SubstitutionModelInterface> getModel(size_t i)
  {
    return models_[i];
  }

  std::shared_ptr<const WordAlphabet> getWordAlphabet()
  {
    return wordAlphabet_;
  }
};


/**
 * @brief Abstract Basal class for words of substitution models.
 * @author Laurent Guéguen
 *
 * Objects of this class are built from several substitution models.
 * Each model corresponds to a position in the word. No model is
 * directly accessible.
 *
 * Only substitutions with one letter changed are accepted.
 *
 * There is one substitution per word per unit of time
 * on the equilibrium frequency, and each position has its specific rate.
 * For example, if there are @f$n@f$ models and \f$\rho_i\f$ is the rate of
 * model i (@f$\sum_{i=1}^{n} \rho_i = 1@f$):
 * @f{eqnarray*}
 * Q_{abc \rightarrow abd} &=& \rho_2 Q^{(2)}_{c \rightarrow d}\\
 * Q_{abc \rightarrow aed} &=& 0\\
 * @f}
 *
 * The parameters of this word model are the same as the ones of the
 * models used. Their names have a new prefix :
 *
 * If there is one model per position, "i_" where i stands for the
 * position in the word.
 *
 * If there is only one model, "123..._" where all positions are
 * enumerated.
 */
class AbstractWordSubstitutionModel :
  public AbstractSubstitutionModel
{
private:
  /**
   * @brief boolean flag to check if a specific WordAlphabet has been built
   */
  bool newAlphabet_;

protected:
  
  /**
   * Vector of shared_ptr, to allow multiple positions to share the same model.
   */
  std::vector<std::shared_ptr<SubstitutionModelInterface>> VSubMod_;
  
  std::vector<std::string> VnestedPrefix_;

  std::vector<double> Vrate_;

protected:
  void updateMatrices_();

  /**
   * @brief Called by updateMatrices to handle specific modifications
   * for inheriting classes
   */
  virtual void completeMatrices_() = 0;

  /**
   * @brief First fill of the generator, from the position model
   */
  virtual void fillBasicGenerator_();

public:
  /**
   * @brief Build a new AbstractWordSubstitutionModel object from a
   * vector of pointers to SubstitutionModels.
   *
   * @param modelList the list of substitution models to use, in
   *   the order of the positions in the words from left to right. All
   *   the models must be different objects to avoid parameters
   *   redundancy, otherwise only the first model is used. The used models
   *   are owned by the instance.
   * @param prefix the Namespace.
   */
  AbstractWordSubstitutionModel(
    ModelList& modelList,
    const std::string& prefix);

  /**
   * @brief Build a new AbstractWordSubstitutionModel object from a
   * pointer to an SubstitutionModel and a number of
   * desired models.
   *
   * @param pmodel A pointer to the substitution model to use in all
   * the positions. It will be owned by the instance.
   * @param num The number of models involved.
   * @param prefix the Namespace.
   */
  AbstractWordSubstitutionModel(
    std::unique_ptr<SubstitutionModelInterface> model,
    unsigned int num,
    const std::string& prefix);

  AbstractWordSubstitutionModel(const AbstractWordSubstitutionModel&);

  AbstractWordSubstitutionModel& operator=(const AbstractWordSubstitutionModel&);

  virtual ~AbstractWordSubstitutionModel() {}

  void setNamespace(const std::string& prefix);

protected:
  /**
   * @brief Constructor for the derived classes only
   */
  AbstractWordSubstitutionModel(
      std::shared_ptr<const Alphabet> alph,
      std::shared_ptr<const StateMapInterface> stateMap,
      const std::string& prefix);

public:

  /**
   * @brief returns the ith model, or throw an exception if i is not a valid number.
   */
  const SubstitutionModelInterface& nModel(size_t i) const
  {
    if (i < VSubMod_.size())
      return *VSubMod_[i];
    else
      throw NullPointerException("AbstractWordSubstitutionModel::nModel. Invalid model requested.");
  }

  size_t getNumberOfModels() const
  {
    return VSubMod_.size();
  }

  /**
   *@brief Estimation of the parameters of the models so that the
   * equilibrium frequencies match the given ones.
   *
   *@param freqs  map of the frequencies
   *
   * When there is one submodel for all the positions, the submodel
   * parameters are fit on the means of the frequencies on each
   * position. Otherwise, each model is fit on the frequencies on its
   * corresponding position in the word.
   */
  virtual void setFreq(std::map<int, double>& freqs);
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ABSTRACTWORDSUBSTITUTIONMODEL_H
