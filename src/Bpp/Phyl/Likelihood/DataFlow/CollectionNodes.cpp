//
// File: CollectionNodes.cpp
// Authors:
//   Laurent Guéguen (2018)
// Created: mercredi 8 avril 2020, ÃÂ  00h 31
//

#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>

#include "Bpp/Phyl/Likelihood/DataFlow/ProcessTree.h"
#include "CollectionNodes.h"

using namespace std;
using namespace bpp;

CollectionNodes::CollectionNodes(
    Context& context,
    shared_ptr<const SubstitutionProcessCollection> collection) :
  AbstractParametrizable(""),
  collection_(collection), context_(context)
{
  // add Independent Parameters
  const auto& paramProc = collection_->getIndependentParameters();

  for (size_t i = 0; i < paramProc.size(); ++i)
  {
    shareParameter_(ConfiguredParameter::create(this->context(), paramProc[i]));
  }

  // Share dependencies with aliased parameters

  for (size_t i = 0; i < paramProc.size(); ++i)
  {
    auto vs = collection_->getAlias(paramProc[i].getName());
    auto dep = dynamic_cast<const ConfiguredParameter*>(&getParameter(paramProc[i].getName()))->dependency(0);
    for (const auto& s:vs)
    {
      auto newacp = ConfiguredParameter::create(this->context(), {dep}, collection_->getParameter(s));
      shareParameter_(newacp);
    }
  }


  // rates nodes
  auto rN = collection_->getRateDistributionNumbers();

  for (auto num : rN)
  {
    std::string suff = "_" + TextTools::toString(num);

    auto obj = collection_->getRateDistribution(num);

    if (!dynamic_pointer_cast<const ConstantRateDistribution>(obj))
      distColl_.addObject(ConfiguredParametrizable::createConfigured<DiscreteDistribution, ConfiguredDistribution>(this->context(), *obj, getParameters_(), suff), num);
  }

  // models nodes
  auto mN = collection_->getModelNumbers();

  for (auto num : mN)
  {
    std::string suff = "_" + TextTools::toString(num);

    const auto& obj = collection_->getModel(num);

    modelColl_.addObject(ConfiguredParametrizable::createConfigured<BranchModelInterface, ConfiguredModel>(context_, *obj, getParameters_(), suff), num);
  }

  // frequencies nodes
  auto fN = collection_->getFrequenciesNumbers();

  for (auto num : fN)
  {
    std::string suff = "_" + TextTools::toString(num);

    const auto& obj = collection_->frequencySet(num);

    freqColl_.addObject(ConfiguredParametrizable::createConfigured<FrequencySetInterface, ConfiguredFrequencySet>(context_, obj, getParameters_(), suff), num);
  }

  ///////
  // tree nodes
  auto tN = collection_->getTreeNumbers();

  for (auto num : tN)
  {
    std::string suff = "_" + TextTools::toString(num);

    auto obj = collection_->getTree(num);

    treeColl_.addObject(make_shared<ProcessTree>(context_, *obj, getParameters_(), suff), num);
  }
}

std::shared_ptr<ProcessTree> CollectionNodes::getProcessTree(size_t treeIndex)
{
  return std::dynamic_pointer_cast<ProcessTree>(treeColl_[treeIndex]);
}
