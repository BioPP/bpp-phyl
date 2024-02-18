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


#include "AbstractAutonomousSubstitutionProcess.h"

// From the stl:
#include <memory>

namespace bpp
{
/**
 * @brief Space and time homogeneous substitution process, without mixture.
 */
class SimpleSubstitutionProcess :
  public AbstractAutonomousSubstitutionProcess
{
protected:
  std::shared_ptr<BranchModelInterface> model_;

public:
  SimpleSubstitutionProcess(
      std::shared_ptr<BranchModelInterface> model,
      std::shared_ptr<const PhyloTree> tree = nullptr,
      std::shared_ptr<FrequencySetInterface> rootFrequencies = nullptr);

  SimpleSubstitutionProcess(
      std::shared_ptr<BranchModelInterface> model,
      std::shared_ptr<ParametrizablePhyloTree> tree,
      std::shared_ptr<FrequencySetInterface> rootFrequencies = nullptr);

  SimpleSubstitutionProcess(const SimpleSubstitutionProcess& ssp);

  SimpleSubstitutionProcess& operator=(const SimpleSubstitutionProcess& ssp);

public:
  SimpleSubstitutionProcess* clone() const override { return new SimpleSubstitutionProcess(*this); }

  const StateMapInterface& stateMap() const override
  {
    return model_->stateMap();
  }

  std::shared_ptr<const StateMapInterface> getStateMap() const override
  {
    return model_->getStateMap();
  }

  size_t getNumberOfModels() const override { return 1; }

  std::vector<size_t> getModelNumbers() const override
  {
    return std::vector<size_t>(1, 1);
  }

  const BranchModelInterface& model(unsigned int nodeId, size_t classIndex) const override
  {
    return *model_;
  }

  std::shared_ptr<const BranchModelInterface> getModel(unsigned int nodeId, size_t classIndex) const override
  {
    return model_;
  }

  const BranchModelInterface& model(size_t n) const override
  {
    return *model_;
  }

  std::shared_ptr<const BranchModelInterface> getModel(size_t n) const override
  {
    return model_;
  }

  const std::vector<unsigned int> getNodesWithModel(size_t i) const override
  {
    throw Exception("SimpleSubstitutionProcess::getNodesWithModel not finished. Ask developpers.");
    return Vuint(0);
  }

  size_t getModelNumberForNode(unsigned int nodeId) const override
  {
    return 1;
  }

  virtual std::shared_ptr<const BranchModelInterface> getModelForNode(unsigned int nodeId) const override
  {
    return model_;
  }

  std::shared_ptr<const DiscreteDistributionInterface> getRateDistribution() const override
  {
    return nullptr;
  }

  std::shared_ptr<DiscreteDistributionInterface> getRateDistribution() override
  {
    return nullptr;
  }

  const DiscreteDistributionInterface& rateDistribution() const override
  {
    throw NullPointerException("SimpleSubstitutionProcess::rateDistribution. This process does not model a rate distribution.");
  }

  DiscreteDistributionInterface& rateDistribution() override
  {
    throw NullPointerException("SimpleSubstitutionProcess::rateDistribution. This process does not model a rate distribution.");
  }

 ParameterList getSubstitutionModelParameters(bool independent) const override
  {
    return independent ? model_->getIndependentParameters() : model_->getParameters();
  }

  ParameterList getRateDistributionParameters(bool independent) const override
  {
    return ParameterList();
  }

  ParameterList getBranchLengthParameters(bool independent) const override
  {
    if (getParametrizablePhyloTree())
      return getParametrizablePhyloTree()->getParameters();
    else
      return ParameterList();
  }

  const std::vector<double>& getRootFrequencies() const override
  {
    if (!hasRootFrequencySet())
    {
      if (std::dynamic_pointer_cast<const TransitionModelInterface>(model_))
        return std::dynamic_pointer_cast<const TransitionModelInterface>(model_)->getFrequencies();
      else
        throw Exception("SimpleSubstitutionProcess::getRootFrequencies not possible with a non Transition Model.");
    }
    else
      return getRootFrequencySet()->getFrequencies();
  }

  /**
   * @brief Set the modelPath, after checking  it is valid
   * (ie modelpath has only the model of the process).
   *
   */

  void setModelScenario(std::shared_ptr<ModelScenario> modelpath) override;

  double getProbabilityForModel(size_t classIndex) const override
  {
    if (classIndex != 0)
      throw IndexOutOfBoundsException("SimpleSubstitutionProcess::getProbabilityForModel.", classIndex, 0, 1);
    return 1;
  }

  Vdouble getClassProbabilities() const override
  {
    return Vdouble(1, 1);
  }

  double getRateForModel(size_t classIndex) const override
  {
    if (classIndex != 0)
      throw IndexOutOfBoundsException("SimpleSubstitutionProcess::getRateForModel.", classIndex, 0, 1);
    return 1;
  }

protected:
  void fireParameterChanged(const ParameterList& pl) override; // Forward parameters and updates probabilities if needed.

};
} // end namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_SIMPLESUBSTITUTIONPROCESS_H
