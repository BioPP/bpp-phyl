//
// File: CollectionNodes.cpp
// Authors:
//   Laurent GuÃ©guen (2018)
// Created: mercredi 8 avril 2020, Ã  00h 31
//

#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>

#include "Bpp/Phyl/NewLikelihood/DataFlow/ProcessTree.h"
#include "CollectionNodes.h"

using namespace std;
using namespace bpp;

CollectionNodes::CollectionNodes(Context& context,
                                 const SubstitutionProcessCollection& collection) :
  AbstractParametrizable(""),
  collection_(collection), context_(context)
{
  // add Independent Parameters
  const auto& paramProc = collection_.getIndependentParameters();

  for (size_t i = 0; i < paramProc.size(); i++)
  {
    shareParameter_(ConfiguredParameter::create(getContext(), paramProc[i]));
  }

  // Share dependencies with aliased parameters

  for (size_t i = 0; i < paramProc.size(); i++)
  {
    auto vs = collection_.getAlias(paramProc[i].getName());
    auto dep = dynamic_cast<const ConfiguredParameter*>(&getParameter(paramProc[i].getName()))->dependency(0);
    for (const auto& s:vs)
    {
      auto newacp = ConfiguredParameter::create(getContext(), {dep}, collection_.getParameter(s));
      shareParameter_(newacp);
    }
  }


  // rates nodes
  auto rN = collection_.getRateDistributionNumbers();

  for (auto num : rN)
  {
    std::string suff = "_" + TextTools::toString(num);

    const auto& obj = collection_.getRateDistribution(num);

    if (dynamic_cast<const ConstantRateDistribution*>(&obj) == nullptr)
      distColl_.addObject(ConfiguredParametrizable::createConfigured<DiscreteDistribution, ConfiguredDistribution>(getContext(), obj, getParameters_(), suff), num);
  }

  // models nodes
  auto mN = collection_.getModelNumbers();

  for (auto num : mN)
  {
    std::string suff = "_" + TextTools::toString(num);

    const auto& obj = collection_.getModel(num);

    modelColl_.addObject(ConfiguredParametrizable::createConfigured<BranchModel, ConfiguredModel>(context_, *obj, getParameters_(), suff), num);
  }

  // frequencies nodes
  auto fN = collection_.getFrequenciesNumbers();

  for (auto num : fN)
  {
    std::string suff = "_" + TextTools::toString(num);

    const auto& obj = collection_.getFrequencies(num);

    freqColl_.addObject(ConfiguredParametrizable::createConfigured<FrequencySet, ConfiguredFrequencySet>(context_, obj, getParameters_(), suff), num);
  }

  ///////
  // tree nodes
  auto tN = collection_.getTreeNumbers();

  for (auto num : tN)
  {
    std::string suff = "_" + TextTools::toString(num);

    const auto& obj = collection_.getTree(num);

    treeColl_.addObject(make_shared<ProcessTree>(context_, obj, getParameters_(), suff), num);
  }
}

std::shared_ptr<ProcessTree> CollectionNodes::getProcessTree(size_t treeIndex)
{
  return std::dynamic_pointer_cast<ProcessTree>(treeColl_[treeIndex]);
}
