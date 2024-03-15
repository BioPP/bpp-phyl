// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_RATEACROSSSITESSUBSTITUTIONPROCESS_H
#define BPP_PHYL_LIKELIHOOD_RATEACROSSSITESSUBSTITUTIONPROCESS_H


#include "AbstractAutonomousSubstitutionProcess.h"

// From bpp-core:
#include <Bpp/Numeric/AbstractParameterAliasable.h>

// From the stl:
#include <memory>

namespace bpp
{
class RateAcrossSitesSubstitutionProcess :
  public AbstractAutonomousSubstitutionProcess
{
private:
  std::shared_ptr<BranchModelInterface> model_;
  std::shared_ptr<DiscreteDistributionInterface> rDist_;

public:
  RateAcrossSitesSubstitutionProcess(
    std::shared_ptr<BranchModelInterface> model,
    std::shared_ptr<DiscreteDistributionInterface> rdist,
    std::shared_ptr<const PhyloTree> tree = nullptr,
    std::shared_ptr<FrequencySetInterface> rootFrequencies = nullptr);

  RateAcrossSitesSubstitutionProcess(
    std::shared_ptr<BranchModelInterface> model,
    std::shared_ptr<DiscreteDistributionInterface> rdist,
    std::shared_ptr<ParametrizablePhyloTree> tree,
    std::shared_ptr<FrequencySetInterface> rootFrequencies = nullptr);

  RateAcrossSitesSubstitutionProcess(const RateAcrossSitesSubstitutionProcess& rassp);

  RateAcrossSitesSubstitutionProcess& operator=(const RateAcrossSitesSubstitutionProcess& rassp);

public:
  RateAcrossSitesSubstitutionProcess* clone() const override
  { 
    return new RateAcrossSitesSubstitutionProcess(*this);
  }

  size_t getNumberOfModels() const override { return 1; }

  std::vector<size_t> getModelNumbers() const override
  {
    return std::vector<size_t>(1, 1);
  }

  const StateMapInterface& stateMap() const override
  {
    return model_->stateMap();
  }

  std::shared_ptr<const StateMapInterface> getStateMap() const override
  {
    return model_->getStateMap();
  }

  const BranchModelInterface& model(size_t n) const override
  {
    return *model_;
  }

  std::shared_ptr<const BranchModelInterface> getModel(size_t n) const override
  {
    return model_;
  }

  const BranchModelInterface& model(unsigned int nodeId, size_t classIndex) const override
  {
    return *model_;
  }

  std::shared_ptr<const BranchModelInterface> getModel(unsigned int nodeId, size_t classIndex) const override
  {
    return model_;
  }

  const std::vector<unsigned int> getNodesWithModel(size_t i) const override
  {
    throw Exception("RateAcrossSitesSubstitutionProcess::getNodesWithModel not finished. Ask developpers.");
    return Vuint(0);
  }

  size_t getModelNumberForNode(unsigned int nodeId) const override
  {
    return 1;
  }

  std::shared_ptr<const BranchModelInterface> getModelForNode(unsigned int nodeId) const override
  {
    return model_;
  }

  const DiscreteDistributionInterface& rateDistribution() const override
  {
    if (!rDist_)
      throw NullPointerException("RateAcrossSitesSubstitutionProcess::rateDistribution. No associated rate distribution.");
    return *rDist_;
  }

  DiscreteDistributionInterface& rateDistribution() override
  {
    if (!rDist_)
      throw NullPointerException("RateAcrossSitesSubstitutionProcess::rateDistribution. No associated rate distribution.");
    return *rDist_;
  }

  std::shared_ptr<const DiscreteDistributionInterface> getRateDistribution() const override
  {
    return rDist_;
  }

  std::shared_ptr<DiscreteDistributionInterface> getRateDistribution() override
  {
    return rDist_;
  }

  ParameterList getSubstitutionModelParameters(bool independent) const override
  {
    return independent ? model_->getIndependentParameters() : model_->getParameters();
  }

  ParameterList getRateDistributionParameters(bool independent) const override
  {
    return independent ? rDist_->getIndependentParameters() : rDist_->getParameters();
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
    if (classIndex >= rDist_->getNumberOfCategories())
      throw IndexOutOfBoundsException("RateAcrossSitesSubstitutionProcess::getProbabilityForModel.", classIndex, 0, rDist_->getNumberOfCategories());
    return rDist_->getProbability(classIndex);
  }

  Vdouble getClassProbabilities() const override
  {
    Vdouble vProb;

    for (size_t i = 0; i < rDist_->getNumberOfCategories(); ++i)
    {
      vProb.push_back(rDist_->getProbability(i));
    }

    return vProb;
  }

  double getRateForModel(size_t classIndex) const override
  {
    if (classIndex >= rDist_->getNumberOfCategories())
      throw IndexOutOfBoundsException("RateAcrossSitesSubstitutionProcess::getRateForModel.", classIndex, 0, rDist_->getNumberOfCategories());
    return rDist_->getCategory(classIndex);
  }

protected:
  void fireParameterChanged(const ParameterList& pl) override; // Forward parameters and updates probabilities if needed.
};
} // end namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_RATEACROSSSITESSUBSTITUTIONPROCESS_H
