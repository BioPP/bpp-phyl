//
// File: CollectionNodes.h
// Authors: Laurent Guéguen
// Created: mardi 7 avril 2020, à 23h 52
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef COLLECTION_NODES_H
#define COLLECTION_NODES_H

#include "Bpp/Phyl/NewLikelihood/DataFlow/Model.h"
#include "Bpp/Phyl/NewLikelihood/DataFlow/DiscreteDistribution.h"
#include "Bpp/Phyl/NewLikelihood/DataFlow/FrequencySet.h"
#include "Bpp/Phyl/NewLikelihood/DataFlow/ProcessTree.h"
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlow.h>

#include "Bpp/Phyl/NewLikelihood/SubstitutionProcessCollection.h"

namespace bpp
{
/** Construction of all the DataFlow objects linked with objects in
 * a SubstitutionProcessCollection.
 *
 */

class CollectionNodes :
  public AbstractParametrizable
{
private:
  const SubstitutionProcessCollection& collection_;

  Context& context_;

  /**
   * A collection of Branch Models
   */

  ParametrizableCollection<ConfiguredModel> modelColl_;

  /*
   * A collection of Frequencies Sets
   */

  ParametrizableCollection<ConfiguredFrequencySet> freqColl_;

  /*
   * A collection of DiscreteDistributions
   */

  ParametrizableCollection<ConfiguredDistribution> distColl_;

  /*
   * A collection of trees
   *
   */

  ParametrizableCollection<ProcessTree> treeColl_;

public:
  CollectionNodes(Context& context,
                  const SubstitutionProcessCollection& collection);

  CollectionNodes* clone() const
  {
    throw Exception("CollectionNodes::clone should not be called.");
  }

  Context& getContext()
  {
    return context_;
  }

  const SubstitutionProcessCollection& getCollection() const
  {return collection_;}

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

#endif// LIKELIHOOD_CALCULATION_SINGLE_PROCESS_H
