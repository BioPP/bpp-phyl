// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_COLLECTIONNODES_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_COLLECTIONNODES_H

#include <Bpp/Phyl/Likelihood/DataFlow/DataFlow.h>

#include "Bpp/Phyl/Likelihood/DataFlow/DiscreteDistribution.h"
#include "Bpp/Phyl/Likelihood/DataFlow/FrequencySet.h"
#include "Bpp/Phyl/Likelihood/DataFlow/Model.h"
#include "Bpp/Phyl/Likelihood/DataFlow/ProcessTree.h"
#include "Bpp/Phyl/Likelihood/SubstitutionProcessCollection.h"

namespace bpp
{
/**
 * @biref Construction of all the DataFlow objects linked with objects in
 * a SubstitutionProcessCollection.
 */
class CollectionNodes :
  public AbstractParametrizable
{
private:
  std::shared_ptr<const SubstitutionProcessCollection> collection_;

  Context& context_;

  /**
   * @brief A collection of Branch Models
   */
  ParametrizableCollection<ConfiguredModel> modelColl_;

  /**
   * @brief A collection of Frequencies Sets
   */
  ParametrizableCollection<ConfiguredFrequencySet> freqColl_;

  /**
   * @brief A collection of DiscreteDistributions
   */
  ParametrizableCollection<ConfiguredDistribution> distColl_;

  /**
   * @brief A collection of trees
   */
  ParametrizableCollection<ProcessTree> treeColl_;

public:
  CollectionNodes(
      Context& context,
      std::shared_ptr<const SubstitutionProcessCollection> collection);

  CollectionNodes* clone() const
  {
    throw Exception("CollectionNodes::clone should not be called.");
  }

  Context& context()
  {
    return context_;
  }

  const SubstitutionProcessCollection& collection() const
  {
    return *collection_;
  }

  std::shared_ptr<const SubstitutionProcessCollection> getCollection() const
  {
    return collection_;
  }

  ConfiguredModel& model(size_t modelIndex)
  {
    return dynamic_cast<ConfiguredModel&>(*modelColl_[modelIndex]);
  }

  std::shared_ptr<ConfiguredModel> getModel(size_t modelIndex)
  {
    return std::dynamic_pointer_cast<ConfiguredModel>(modelColl_[modelIndex]);
  }

  ParametrizableCollection<ConfiguredModel>& getModelCollection()
  {
    return modelColl_;
  }

  std::shared_ptr<ConfiguredFrequencySet> getFrequencies(size_t freqIndex)
  {
    return std::dynamic_pointer_cast<ConfiguredFrequencySet>(freqColl_[freqIndex]);
  }

  std::shared_ptr<ConfiguredDistribution> getRateDistribution(size_t distIndex)
  {
    return std::dynamic_pointer_cast<ConfiguredDistribution>(distColl_[distIndex]);
  }

  std::shared_ptr<ProcessTree> getProcessTree(size_t treeIndex);
};
} // namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_COLLECTIONNODES_H
