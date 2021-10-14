//
// File: SimpleSubstitutionProcess.h
// Authors:
//   Julien Dutheil
// Created: 2012-05-15 13:11:00
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

#ifndef BPP_PHYL_LIKELIHOOD_SIMPLESUBSTITUTIONPROCESS_H
#define BPP_PHYL_LIKELIHOOD_SIMPLESUBSTITUTIONPROCESS_H


#include "AbstractSubstitutionProcess.h"

// From the stl:
#include <memory>

namespace bpp
{
/**
 * @brief Space and time homogeneous substitution process, without mixture.
 */
class SimpleSubstitutionProcess :
  public AbstractSubstitutionProcess
{
protected:
  std::shared_ptr<BranchModel> model_;

private:

public:
  SimpleSubstitutionProcess(std::shared_ptr<BranchModel> model, ParametrizablePhyloTree* tree);

  SimpleSubstitutionProcess(const SimpleSubstitutionProcess& ssp);

  SimpleSubstitutionProcess& operator=(const SimpleSubstitutionProcess& ssp);

public:
  SimpleSubstitutionProcess* clone() const { return new SimpleSubstitutionProcess(*this); }

  size_t getNumberOfStates() const { return model_->getNumberOfStates(); }

  const StateMap& getStateMap() const
  {
    return model_->getStateMap();
  }

  size_t getNumberOfModels() const { return 1; }

  bool isCompatibleWith(const AlignedValuesContainer& data) const
  {
    return data.getAlphabet()->getAlphabetType() == model_->getAlphabet()->getAlphabetType();
  }

  std::vector<size_t> getModelNumbers() const
  {
    return std::vector<size_t>(1, 1);
  }

  const BranchModel* getModel(unsigned int nodeId, size_t classIndex) const
  {
    return model_.get();
  }

  const BranchModel* getModel(size_t n) const
  {
    return model_.get();
  }

  const std::vector<unsigned int> getNodesWithModel(size_t i) const
  {
    throw Exception("SimpleSubstitutionProcess::getNodesWithModel not finished. Ask developpers.");
    return Vuint(0);
  }

  size_t getModelNumberForNode(unsigned int nodeId) const
  {
    return 1;
  }

  const BranchModel* getModelForNode(unsigned int nodeId) const
  {
    return model_.get();
  }

  const DiscreteDistribution* getRateDistribution() const
  {
    return 0;
  }

  // const Matrix<double>& getGenerator(unsigned int nodeId, size_t classIndex) const {
  //   return model_->getGenerator();
  // }

  std::shared_ptr<const FrequencySet> getRootFrequencySet() const
  {
    return std::shared_ptr<const FrequencySet>(0);
  }

  ParameterList getSubstitutionModelParameters(bool independent) const
  {
    return independent ? model_->getIndependentParameters() : model_->getParameters();
  }

  ParameterList getRateDistributionParameters(bool independent) const
  {
    return ParameterList();
  }

  ParameterList getRootFrequenciesParameters(bool independent) const
  {
    return ParameterList();
  }

  ParameterList getBranchLengthParameters(bool independent) const
  {
    return getParametrizablePhyloTree().getParameters();
  }

  const std::vector<double>& getRootFrequencies() const
  {
    if (std::dynamic_pointer_cast<const TransitionModel>(model_))
      return std::dynamic_pointer_cast<const TransitionModel>(model_)->getFrequencies();
    else
      throw Exception("SimpleSubstitutionProcess::getRootFrequencies not possible with a non Transition Model.");
  }

  /**
   * @brief Set the modelPath, after checking  it is valid
   * (ie modelpath has only the model of the process).
   *
   */

  void setModelScenario(std::shared_ptr<ModelScenario> modelpath);

  double getInitValue(size_t i, int state) const
  {
    return model_->getInitValue(i, state);
  }

  double getProbabilityForModel(size_t classIndex) const
  {
    if (classIndex != 0)
      throw IndexOutOfBoundsException("SimpleSubstitutionProcess::getProbabilityForModel.", classIndex, 0, 1);
    return 1;
  }

  Vdouble getClassProbabilities() const
  {
    return Vdouble(1, 1);
  }

  double getRateForModel(size_t classIndex) const
  {
    if (classIndex != 0)
      throw IndexOutOfBoundsException("SimpleSubstitutionProcess::getRateForModel.", classIndex, 0, 1);
    return 1;
  }

  // bool transitionProbabilitiesHaveChanged() const { return true; }

protected:
  void fireParameterChanged(const ParameterList& pl); // Forward parameters and updates probabilities if needed.

  bool modelChangesWithParameter_(size_t i, const ParameterList& pl) const
  {
    if (pl.getCommonParametersWith(model_->getParameters()).size() > 0)
      return true;
    ParameterList plbl = pTree_->getBranchLengthParameters(i);
    if (plbl.getCommonParametersWith(pl).size() > 0)
      return true;
    return false;
  }
};
} // end namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_SIMPLESUBSTITUTIONPROCESS_H
